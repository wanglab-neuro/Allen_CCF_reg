classdef cellcount
    properties
        ccf = [];
        metadata = [];
    end
    
    methods
        function cc = cellcount()
            disp('Done creating obj')
        end
        
        function [cc, ref] = prep_ccf(cc, tv, scl)
            if nargin < 3 || isempty(scl)
                scl = 0.25;
            end
            [nn, nh, nw] = size(tv);
            
            tvt = griddedInterpolant(single(tv), 'linear');
            sclr = round(1 / scl);
            [X, Y, Z] = ndgrid(1 + sclr / 2: sclr: nh, 1 + sclr / 2: sclr: nw, 1: nn);
            ref = double(tvt(Z, X, Y));
        end
        
    end
       
    methods (Static)
        %% main procedure %%
        function [tv, av, st] = load_ccf()
            if ~exist('tv', 'var')
                m = matfile('C:\Jinghao\research_temp\allenccf\AllenAtlas.mat');
                tv = m.tv;
                av = m.av;
                st = m.st;
            end
        end
    
        function cell_count(cc)
            %%% parameters %%%
            lens_flag = true;
            [file_name, path_name] = cc.data_info();
            cc.metadata.pname = path_name;

            % %%%% if you want a pop up window, use this: %%%%
            % nlist = cc.section_list(length(file_name));

            %%%% if not, put numbers in here (and don't run the above line): %%%%
            % nlist = [570: 8: 632];
            nlist = [678: 8: 742];

            %%% initialize %%%
            [tv, av, st] = cc.load_ccf();
            [cc, ref] = cc.prep_ccf(tv);

            %%% main loop per slice %%%
            cc = cc.load(cc);

            if isempty(cc)
                %%% load slices %%%
                cc = cc.load_slice(cc, file_name);

                %%% pyramid size selection %%%
                sl = cc.metadata.slice;
                ref_flag = 1;
                [gpy, loc, rsclp] = cc.find_loc_pyramid(cc, ref, sl, ref_flag, nlist);

                %%% coarse selection %%%
                idf = cc.find_loc_coarse(cc, ref, sl, ref_flag, nlist, gpy, loc, rsclp);

                %%% fine/rotate selection %%%
                rtl = size(tv, 2) / size(ref, 1);
                cct = cellcount();
                [theta, gamma, ap, rloc, imref] = cc.find_loc_fine(cct, ref, sl, ref_flag, idf, gpy, loc, rtl, rsclp);

                %%% nonrigid register %%%
                angles = struct('theta', theta, 'gamma', gamma, 'ap', ap + idf);
                [img, D, rlocn, imr, coordf] = cc.regist(cc, tv, sl, angles, rloc, imref);

                %%% Detect labeled neurons %%%
                opt.vis = false;
                % opt.vis = true;
                pl = cc.cell_detect(cc, sl, opt);

                %%% transform detected neurons %%%
                opt.stv = cellfun(@size, img, 'uniformoutput', false);
                opt.ssl = cellfun(@size, sl, 'uniformoutput', false);
                pln = cc.cell_transform(D, pl, rlocn, opt);

                %%% get 3d final point list %%%
                pl3 = cc.cell3d(pln, coordf);

                %%% metadata compile %%%
                metadata.slice = sl;
                metadata.pyramid_size = gpy;
                metadata.angles = angles;
                metadata.point_lists = pl;
                metadata.point_lists_new = pln;
                metadata.point_lists_3d = pl3;
                metadata.section_location = rlocn;
                metadata.reg_image = img;
                metadata.displacement = D;
                metadata.reg_ref_image = imr;
                metadata.coordf = coordf;
                metadata.fnames = file_name;

                %%% pass metadata %%%
                cc.metadata = metadata;
                cc.metadata.pname = path_name;

                %%% save data %%%
                cc.save(cc);
            end

            if lens_flag
                %%% find lens lesion %%%
                img = cc.metadata.reg_image;
                imr = cc.metadata.reg_ref_image;
                rlocn = cc.metadata.section_location;
                coordf = cc.metadata.coordf;
                [mask, pll] = cc.lens_loc(cc, img, imr, rlocn);
                pll3 = cc.cell3d(pll, coordf);
            else
                pll = NaN;
                pll3 = NaN;
            end

            %%% lens data compile %%%
            lensdata.point_lists = pll;
            lensdata.point_lists_3d = pll3;
            lensdata.mask = mask;
            cc.save(cc, lensdata)
        end
        
        %% aux functions: main procedure %%
        function [file_name, path_name] = data_info()
            [file_name_tmp, path_name] = uigetfile('*', 'Select imaging file', 'MultiSelect', 'on');
            if ~iscell(file_name_tmp)
                file_name{1} = file_name_tmp;
            else
                file_name = file_name_tmp;
            end
        end
        
        function nlist = section_list(k)
            prompt = ['put in the rough # section for each section, total ', num2str(k), ' sections'];
            nlist = inputdlg(prompt);
            nlist = nlist{1};
        end
        
        function img = load_slice_unit(fname)
            img = imread(fname);
        end
        
        function imgo = prep_slice(img)
            a = normalize(double(img));
            b = normalize(a(:, :, 3));
%             b1 = imgaussfilt(b, floor(min(size(b)) / 100));
%             b1t = normalize(b1) > 0.1;
            b1t = normalize(b) > 0.1;
            b1t = imfill(b1t, 'holes');
            [l, n] = bwlabeln(b1t);
            s = zeros(1, n);
            for i = 1: n
                tmp = l == i;
                s(i) = sum(tmp(:));
            end
            [~, id] = max(s);
            b1 = l == id;
%             b2 = sum(b1 > 0.1, 1);
%             b3 = sum(b1 > 0.1, 2);
            b2 = sum(b1, 1);
            b3 = sum(b1, 2);
            bb = a(find(b3, 1): find(b3, 1, 'last'), find(b2, 1): find(b2, 1, 'last'), :);
            imgo = normalize(bb);
        end
        
        function imgs = select_chan(img, nchan)
            if nargin < 2 || isempty(nchan)
                nchan = 3;
            end
            
            if iscell(img)
                if size(img{1}, 3) > 1
                    imgs = cellfun(@(x) x(:, :, nchan), img, 'uniformoutput', false);
                end
            else
                imgs = {img(:, :, nchan)};
            end
        end
        
        function [gpy, loc, rscl] = find_loc_pyramid(cc, ref, img, ref_flag, n)
            %%% parameters %%%
            if nargin < 4 || isempty(ref_flag)
                ref_flag = 1;
            end
            
            img = cc.select_chan(img);
            
            rt = linspace(0.3, 1, 40);
            istep = 0.5; %%% mm %%%
            res = 0.01; %% mm %%
            istepf = istep / res;
            [nh, nw, nn] = size(ref);
            srg = 5;
            if nargin < 5 || isempty(n)
                rg = repmat(2: istepf: nn - 1, length(img), 1);
            else
                rg = NaN(length(n), 2 * srg + 1);
                for i = 1: length(n)
                    idts = max(2, n(i) - srg) - n(i) + srg;
                    idtt = max(2, n(i) - srg): min(nn - 1, n(i) + srg);
                    idt = idts + 1: idts + length(idtt);
                    rg(i, idt) = idtt;
                end
            end
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            rscl = opt.rscl;
            opt.cfflag = 1; %%% 1 for pyramid, 2 for coarse and 3 for fine %%%
                        
            imgp = img;
            for i = 1: length(imgp)
                %%% slice prep %%%
                reft = ref(:, :, 1);
                reft = imresize(reft, opt.rscl);
                img1 = imresize(imgp{i}, size(reft));
                se = 5;
                imgt = cc.anidenoise(cc, img1, se);
                imgt = imgt(se + 1: end - se, se + 1: end - se);
                imgt = imgaussfilt(TVL1denoise(normalize(imgt), 0.4, 20) .^ 1, max(1, floor(min(size(reft)) / 100)));
                imgt = imresize(imgt, size(reft));
                imgp{i} = imgt;
            end
            
            %%% calculate measures %%%
            if length(imgp) == 1
                cs = zeros(length(rg), length(rt));
                ss = cs;
                mis = cs;
                ks = cs;
                locs = cell(length(rt), 1);
                for i = 1: length(rg)
                    for j = 1: length(rt)
                        bf = imresize(imgp{1}, rt(j));
                        opt.bf = bf;
                        [cst, sst, mit, locst, kst] = cc.mea_cal_unit(cc, ref(:, :, rg(i)), opt);
                        cs(i, j) = cst;
                        ss(i, j) = sst;
                        mis(i, j) = mit;
                        ks(i, j) = kst;
                        locs{j}(i, :) = locst;
                    end
                end
            else
                cs = NaN(size(rg, 2), length(rt), size(rg, 1));
                ss = cs;
                mis = cs;
                ks = cs;
                locs = cell(length(rt), size(rg, 1));
                locs = cellfun(@(x) NaN(size(rg, 2), 4), locs, 'uniformoutput', false);
                
                %%% random arrangement %%%
                idt = randsample(numel(cs), round(numel(cs) / 5));
                [i1, i2, i3] = ind2sub([size(rg, 1), size(rg, 2), length(rt)], idt);
                count = 1;
                denom = 10;
                for i = 1: length(i1)
                    bf = imresize(imgp{i1(i)}, rt(i3(i)));
                    opt.bf = bf;
                    [cst, sst, mit, locst, kst] = cc.mea_cal_unit(cc, ref(:, :, rg(i1(i), i2(i))), opt);
                    cs(i2(i), i3(i), i1(i)) = cst;
                    ss(i2(i), i3(i), i1(i)) = sst;
                    mis(i2(i), i3(i), i1(i)) = mit;
                    ks(i2(i), i3(i), i1(i)) = kst;
                    locs{i3(i), i1(i)}(i2(i), :) = locst;
                    if i > count * length(i1) / denom
                        disp([num2str(i), '/', num2str(length(i1)), ' done'])
                        count = count + 1;
                    end
                end
            end
            
            tmp = ss;
            tmp = nanmax(nanmean(tmp, 3), [], 1);
            [~, i] = max(tmp);
            gpy = rt(i);
            loc = round(nanmean(cell2mat(locs(i, :)'), 1));
        end
        
        function idf = find_loc_coarse(cc, ref, img, ref_flag, n, gpy, loc, rsclp)
            %%% parameters %%%
            if nargin < 4 || isempty(ref_flag)
                ref_flag = 1;
            end
            
            img = cc.select_chan(img);
            
            istep = 0.5; %%% mm %%%
            res = 0.01; %% mm %%
            istepf = istep / res;
            [nh, nw, nn] = size(ref);
            srg = 10;
            if nargin < 5 || isempty(n)
                rg = repmat(2: istepf: nn - 1, length(img), 1);
            else
                rg = NaN(length(n), 2 * srg + 1);
                for i = 1: length(n)
                    idts = max(2, n(i) - srg) - n(i) + srg;
                    idtt = max(2, n(i) - srg): min(nn - 1, n(i) + srg);
                    idt = idts + 1: idts + length(idtt);
                    rg(i, idt) = idtt;
                end
            end
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            lrt = opt.rscl / rsclp;
            opt.cfflag = 2; %%% 1 for pyramid, 2 for coarse and 3 for fine %%%
            opt.loc = loc * lrt;
            
            imgp = img;
            for i = 1: length(imgp)
                %%% slice prep %%%
                reft = ref(:, :, 1);
%                 reft = imresize(reft, opt.rscl);
                img1 = imresize(imgp{i}, size(reft));
                se = 5;
                imgt = cc.anidenoise(cc, img1, se);
                imgt = imgt(se + 1: end - se, se + 1: end - se);
                imgt = imgaussfilt(normalize(TVL1denoise(imgt, 0.4, 20)) .^ 2, max(1, floor(min(size(reft)) / 100)));
                imgt = imresize(imgt, size(reft));
                imgp{i} = imgt;
            end
                        
            %%% calculate measures %%%
            cs = zeros(size(rg));
            ss = cs;
            mis = cs;
            ks = cs;
            for k = 1: length(imgp)
                imgt = imgp{k};
                bf = imresize(imgt, gpy);
                opt.bf = bf;
                for i = 1: length(rg)
                    [cst, sst, mit, ~, kst] = cc.mea_cal_unit(cc, ref(:, :, rg(k, i)), opt);
                    cs(k, i) = cst;
                    ss(k, i) = sst;
                    mis(k, i) = mit;
                    ks(k, i) = kst;
                end
                disp([num2str(k), '/', num2str(size(rg, 1)), ' done'])
            end
            
            tmp = cs + ss;
            [~, i] = max(tmp, [], 2);
            idf = zeros(size(rg, 1), 1);
            for j = 1: size(rg, 1)
                idf(j) = rg(j, i(j));
            end
        end
        
        function [theta, gamma, ap, rloc, imgo] = find_loc_fine(cc, ref, img, ref_flag, n, gpy, loc, rtl, rsclp)
            img = cc.select_chan(img);
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            opt.cfflag = 3;
            lrt = opt.rscl / rsclp;
            opt.loc = lrt * loc;
            if isempty(gcp('nocreate'))
                parpool(feature('numCores'));
            end

            %%% generate interpolation grid %%%
            [nh, nw, nf] = size(ref);
            mx = 2;
            astep = 0.5;
            apmx = 2;
            mxtheta = atand(min(nf - n(end) - apmx, n(1) - 1 - apmx) / (nw / 2));
            mxgamma = atand(min(nf - n(end) - apmx, n(1) - 1 - apmx) / (nh / 2));
            idtheta = max(-mx, -mxtheta): astep: min(mx, mxtheta); %%% left right %%%
            idgamma = max(-mx, -mxgamma): astep: min(mx, mxgamma); %%% up down %%%
            idap = -apmx: apmx;
            nth = length(idtheta);
            nga = length(idgamma);
            nap = length(idap);
            reft = permute(ref, [3, 1, 2]);
            reft = griddedInterpolant(double(reft));
            
            %%% calculate measure %%%
            cs = NaN(nth * nga * nap, length(img));
            ss = cs;
            mis = cs;
            ks = cs;
            ns = nth * nga * nap;
            rlocs = NaN(ns, 4, length(img));
            imgp = img;
            
            %%%% set up arrangement %%%%
            nt = 2;
            id = [repmat((1: size(cs, 1))', nt, 1), randsample(length(img), nt * size(cs, 1), true)];
            coord = cc.prep_coord(nh, nw);
                
            h = tic;
            for ii = 1: length(imgp)
                %%% slice prep %%%
                img1 = imresize(imgp{ii}, [nh, nw]);
                imgt = imgaussfilt(normalize(TVL1denoise(img1, 0.5, 40)) .^ 1, max(1, floor(min(size(reft.Values)) / 100)));
                bf = imresize(imgt, gpy);
                opt.bf = bf;
                
                idc = unique(id(id(:, 2) == ii, 1));
                [i1, i2, i3] = ind2sub([nap, nga, nth], idc);
                imgs = zeros(nh, nw, length(i1));
                for i = 1: length(i1)
                    theta = idtheta(i3(i));
                    gamma = idgamma(i2(i));
                    ap = idap(i1(i)) + n(ii);
                    angles.theta = theta;
                    angles.gamma = gamma;
                    angles.ap = ap;
                    
                    coordf = cc.slice_coords(coord, angles, [nh, nw]);
                    imgt = reft(rtl * (coordf(:, :, 3) - ap) + ap, coordf(:, :, 2), coordf(:, :, 1));
                    imgs(:, :, i) = imgt;
                end
                            
                %%% calculate %%%
                cstt = zeros(length(i1), 1);
                sstt = cstt;
                mistt = cstt;
                kstt = cstt;
                rlocstt = zeros(length(i1), 4);
                parfor i = 1: length(i1)
                    [cst, sst, mit, rloc, kst] = cc.mea_cal_unit(cc, imgs(:, :, i), opt);
                    cstt(i) = cst;
                    sstt(i) = sst;
                    mistt(i) = mit;
                    rlocstt(i, :) = rloc;
                    kstt(i) = kst;
                end
                tt = toc(h);
                disp(['done slice #', num2str(ii), '/', num2str(length(imgp)), ', ', num2str(tt), ' seconds'])
                
                %%% collect %%%
                cs(idc, ii) = cstt;
                ss(idc, ii) = sstt;
                mis(idc, ii) = mistt;
                ks(idc, ii) = kstt;
                rlocs(idc, :, ii) = rlocstt;
            end
                
            tmp = ks .* mis;
            [~, ii] = nanmax(nanmax(tmp, [], 2));
            iu = floor((ii - 1) / (nga * nap));
            ju = floor((ii - iu * nga * nap - 1) / nap);
            ku = mod(ii - iu * nga * nap - ju * nap, nap);
            iu = iu + 1;
            ju = ju + 1;
            if ku == 0
                ku = nap;
            end
                
            theta = idtheta(iu);
            gamma = idgamma(ju);
            ap = idap(ku);
            rloc = round(nanmean(rlocs(ii, :, :), 3));
                
            angles.theta = theta;
            angles.gamma = gamma;
            imgo = imgp;
            for i = 1: length(imgp)
                angles.ap = ap + n(i);
                coordf = cc.slice_coords(coord, angles, [nh, nw]);
                imgo{i} = reft(rtl * (coordf(:, :, 3) - angles.ap) + angles.ap, coordf(:, :, 2), coordf(:, :, 1));
            end
        end
        
        function coordf = slice_coords(coord, angles, sz)
            theta = angles.theta;
            gamma = angles.gamma;
            ap = angles.ap;
            nh = sz(1);
            nw = sz(2);
            
            rx = rotx(gamma);
            ry = roty(theta);
            coordf = ry * rx * coord';
            coordf = reshape(coordf', nh, nw, []);
            coordf(:, :, 1) = coordf(:, :, 1) + (1 + nw) / 2;
            coordf(:, :, 2) = coordf(:, :, 2) + (1 + nh) / 2;
            coordf(:, :, 3) = coordf(:, :, 3) + ap;
        end
        
        function [imgg, D, rloct, imrefs, coordfts] = regist(cc, tv, sl, angles, rloc, imro)
            %%% estimate tv slices to be used %%%
            rt = 0.00;
            [nf, nh, nw] = size(tv);
            coord = cc.prep_coord(nh, nw);
            imgg = imro;
            D = imro;
            coordfts = imro;
            imrefs = imro;
            
            for ii = 1: length(sl)
                anglesuse = angles;
                anglesuse.ap = angles.ap(ii);
                coordft = cc.slice_coords(coord, anglesuse, [nh, nw]);
                coordfts{ii} = coordft;
                rg = [floor(min(min(coordft(:, :, 3)))), ceil(max(max(coordft(:, :, 3))))];
                reft = griddedInterpolant(double(tv(rg(1): rg(2), :, :)));
                coordf = coordft;
                coordf(:, :, 3) = coordft(:, :, 3) - rg(1);
                rloct = round(reshape(rloc, 2, 2)' * [(1 + rt), -rt; -rt, (1 + rt)])';
                rloct = round(nh / size(imro{ii}, 1)) * rloct(:);
                rloct(1) = max(1, rloct(1));
                rloct(3) = max(1, rloct(3));
                rloct(2) = min(nh, rloct(2));
                rloct(4) = min(nw, rloct(4));
                
                %%% imcur %%%
                sls = repmat(zeros(rloct(2) - rloct(1) + 1, rloct(4) - rloct(3) + 1), 1, 1, size(sl{ii}, 3));
                for i = 1: size(sl{ii}, 3)
                    sls(:, :, i) = imresize(sl{ii}(:, :, i), size(sls(:, :, 1)));
                end
                se = 5;
                imcur = cc.anidenoise(cc, sls(:, :, 3), se);
                imcur = imcur(se + 1: end - se, se + 1: end - se);
                imcurt = normalize(imgaussfilt(TVL1denoise(normalize(imcur), 0.4, 20) .^ 1, max(1, floor(min(size(imcur)) / 100))));
                
                %%% imref %%%
                imref = normalize(reft(coordf(:, :, 3), coordf(:, :, 2), coordf(:, :, 1)));
                imref = cc.anidenoise(cc, imref, se);
                imref = imref(se + 1: end - se, se + 1: end - se);
                imrefs{ii} = imref;
                mask = normalize(imref) > 0.05;
                mask = imfill(mask, 'holes');
                imreft = imref .* mask;
                %%%% test exponent %%%%
                imreft = normalize(imreft(rloct(1): rloct(2), rloct(3): rloct(4)));
                %%%% apply exponent %%%%
                imreft = normalize(imgaussfilt(TVL1denoise(normalize(imreft .^ 1), 0.4, 40), max(1, floor(min(size(imref)) / 400))));
                
                imcurt = cc.opt_exp(cc, imcurt, imreft);
                imreft = cc.opt_exp(cc, imreft, imcurt);
                
                %%% sharpen images %%%
                imreft = imsharpen(imreft, 'radius', 2, 'amount', 2);
                imcurt = imsharpen(imcurt, 'radius', 2, 'amount', 2);
                
                imreft = max(0, imreft);
                imcurt = max(0, imcurt);
                
                disp('doing registration')
                try
                    [Dt, img]=imregdemons(gpuArray(imcurt), gpuArray(imreft), 400, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);
                    D = gather(D);
                catch
                    [Dt, img]=imregdemons(imcurt, imreft, 400, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);
                end
                disp('done registration')
                
                %%% collect %%%
                D{ii} = Dt;
                imggt = sls;
                for i = 1: size(sl{ii}, 3)
                    imggt(:, :, i) = cc.iminterpolate(sls(:, :, i), Dt(:, :, 2), Dt(:, :, 1));
                end
                imgg{ii} = imggt;
            end
        end
        
        function pl = cell_detect(cc, a, opt)
            if nargin < 3 || isempty(opt)
                vis = false;
            else
                vis = opt.vis;
            end

            a = cc.select_chan(a, 2);

            pl = cell(length(a), 1);
            for ii = 1: length(a)
                %%% cell enhancing %%%
                aa0 = normalize(double(a{ii}));
                a0 = normalize(imgaussfilt(aa0, 0.5) .^ 2);
                
                nl = 2;
                t = a0;
                tt = imgaussfilt(t, 19);
                for i = 1: nl
                    t = abs(t - tt);
                    tt = imgaussfilt(t, 49);
                end
                
                t = normalize(t);
                nl = 1;
                for i = 1: nl
                    t = normalize(imgaussfilt(double(cc.multi_otsu(cc, t)), 2));
                end
                a1 = t;
                
                nl = 2;
                t = a1;
                tt = imgaussfilt(t, 19);
                for i = 1: nl
                    t = abs(t - tt);
                    tt = imgaussfilt(t, 49);
                end
                [gx, gy] = gradient(t);
                a2 = normalize(imgaussfilt(gx .^ 2 + gy .^ 2, 2));%imgaussfilt(a1, 199) .^ 0.2);
                a3 = sigmoid(a2, 1000, cc.intense_filter(a2, 1));
                
                t = normalize(a3);
                nl = 1;
                for i = 1: nl
                    t = normalize(imgaussfilt(double(cc.multi_otsu(cc, t)), 2));
                end
                t1 = imgaussfilt(a3, 199);
                a4 = normalize(t .* a1 .* a0 .* t1);
                af = normalize(imgaussfilt(a4 .* (a4 > cc.intense_filter(a4, 0.2)), 2));
                
                %%% find cells %%%
                mxs = imregionalmax(af);
                [l, n] = bwlabeln(mxs);
                [y, x] = find(mxs > 0);
                pl{ii} = [y, x];
            end
            
            %%% visualize %%%
            if vis
                figure(gcf)
                nw = 4;
                nh = ceil(length(sl) / nw);
                count = 1;
                for ii = 1: length(sl)
                    subplot(nw, nh, count)
                    imagesc(a)
                    hold on
                    for i = 1: length(y)
                        plot(x(i), y(i), '.r')
                    end
                    hold off
                end
            end
        end
        
        function plns = cell_transform(D, pl, rloc, opt)
            plns = pl;
            for ii = 1: length(pl)
                [nh, nw, nd] = size(D{ii});
                rh = nh / opt.ssl{ii}(1);
                wh = nw / opt.ssl{ii}(2);
                pltmp = round(pl{ii} .* [rh, wh]);
                [x, y] = ndgrid(0: nh - 1, 0: nw - 1); % coordinate image
                x_prime = x + D{ii}(:, :, 2); % updated x values (1st dim, rows)
                y_prime = y + D{ii}(:, :, 1); % updated y values (2nd dim, cols)
                
                %%% second logdemons tranform %%%
                pln = pltmp;
                tmpx = round(x_prime(:));
                tmpy = round(y_prime(:));
                id1 = tmpx > 0 & tmpx <= nh;
                id2 = tmpy > 0 & tmpy <= nw;
                tmpx = tmpx(id1 & id2);
                tmpy = tmpy(id1 & id2);
                idt = sub2ind([nh, nw], tmpx, tmpy);
                lidt = find(id1 & id2);
                for i = 1: size(pltmp, 1)
                    xt = pltmp(i, 1);
                    yt = pltmp(i, 2);
                    idc = sub2ind([nh, nw], xt, yt);
                    [~, idu] = min(abs(idt - idc));
                    idu = idu(1);
                    [x, y] = ind2sub([nh,nw], lidt(idu));
                    pln(i, :) = [x + rloc(1), y + rloc(3)];
                end
                plns{ii} = pln;
            end
        end
        
        function pl3 = cell3d(pl, coordf)
            pl3 = pl;
            for ii = 1: length(pl)
                tmp1 = griddedInterpolant(coordf{ii}(:, :, 1));
                tmp2 = griddedInterpolant(coordf{ii}(:, :, 2));
                tmp3 = griddedInterpolant(coordf{ii}(:, :, 3));
                pl3t = zeros(size(pl{ii}, 1), 3);
                for i = 1: size(pl{ii}, 1)
                    pl3t(i, :) = [tmp1(pl{ii}(i, :)), tmp2(pl{ii}(i, :)), tmp3(pl{ii}(i, :))];
                end
                pl3{ii} = pl3t;
            end
        end
        
        function save(cc, varargin)
            fn = [cc.metadata.pname, filesep, 'cell_count_results.mat'];
            if nargin < 2 || isempty(varargin)
                save(fn, 'cc')
            else
                eval([inputname(2), ' = varargin{1}'])
                save(fn, inputname(2), '-append')
            end
        end
        
        function cc = load(cc)
            fn = [cc.metadata.pname, filesep, 'cell_count_results.mat'];
            if exist(fn, 'file') == 2
                load(fn, 'cc')
            end
        end
        
        function cc = load_slice(cc, file_name)
            for i = 1: length(file_name)
                fname = file_name{i};
                slo = cc.load_slice_unit([cc.metadata.pname, fname]);
                sl = cc.prep_slice(slo);
                cc.metadata.slice{i} = sl;
            end
        end
        
        %% aux functions: lens locate %%
        function [masks, plls] = lens_loc(cc, img, imr, rlocn)
            nthres = 50;
            ps = 5;
            fthres = 10;
            stthres = 0.5;
            
            masks = img;
            plls = img;
            for ii = 1: length(img)
                sthres = size(img{ii}, 1) / 2;
                t = imr{ii}(rlocn(1): rlocn(2), rlocn(3): rlocn(4));
                %             t1 = cc.opt_exp(cc, t, img(:, :, 3));
                %             t2 = cc.opt_exp(cc, img(:, :, 3), t);
                t2 = normalize(imgaussfilt(normalize(TVL1denoise(normalize(img{ii}(:, :, 3)), 0.2, 40)), 5));
                
                %%%%%%% need to reduce nonrelevant regions %%%%%%%%%%
                m1 = normalize(t) > 0.1;
                m2 = t2 > 0.1;
                r1 = imfill(m1 > 0.01, 'holes');
                r2 = imfill(m2 > 0.01, 'holes');
                r2r = padarray(r2, [ps, ps], 0, 'both');
                r1t = padarray(r1, [ps, ps], 0, 'both');
                r2t = padarray(t2, [ps, ps], 0, 'both');
                r2t = imgaussfilt(r2t, 3);
                rf1 = imerode(r1t, strel('disk', 5));
                [l, n] = bwlabeln(rf1);
                s = zeros(n, 1);
                for i = 1: n
                    tmp = l == i;
                    s(i) = sum(tmp(:));
                end
                [~, id] = max(s);
                rf1 = l == id;
                bw = activecontour(r2t, rf1, 500, 'Chan-Vese', 'smoothfactor', 50);
                bw = imerode(bw, strel('disk', 5));
                bd1 = bwboundaries(bw);
                id = bd1{1}(:, 1) < size(img{ii}, 1) / 2 + ps;
                pt = bd1{1}(id, :);
                rf1t = griddedInterpolant(double(r2r));
                b1 = rf1t(pt);
                s = sum(~b1);
                pt = pt - [ps, ps];
                mask = false(size(m2));
            
                if s > fthres
                    %%% for open lesion %%%
                    h1 = ~m2 & bw(ps + 1: end - ps, ps + 1: end - ps);
                    %                 h1 = ~m2 & rf1(ps + 1: end - ps, ps + 1: end - ps);
                    tmpt = false(size(h1));
                    tmpt(sub2ind(size(tmpt), pt(:, 1), pt(:, 2))) = true;
                    hh = tmpt + h1;
                    [l, n] = bwlabeln(h1);
                    h2 = false(size(tmpt));
                    for i = 1: n
                        tmp = l == i;
                        tt = max(hh(tmp));
                        if tt > 1
                            h2 = h2 | tmp;
                        end
                    end
                    
                    [l, n] = bwlabeln(h2);
                    s = zeros(n, 2);
                    for i = 1: n
                        tmp = l == i;
                        [y, x] = find(tmp);
                        s(i, :) = [mean(y), length(y)];
                    end
                    id1 = s(:, 2) > nthres;
                    id2 = s(:, 1) < sthres;
                    id1t = find(id1 & id2);
                    [~, id] = max(s(id1t, 2));
                    id = id1t(id);
                    if ~isempty(id)
                        mask = l == id;
                    else
                        mask = isnan(l);
                    end
                end
                
                %%% for incontinuous/inside opening component %%%
                h1 = imerode(~m1 & rf1(ps + 1: end - ps, ps + 1: end - ps), strel('disk', 10));
                h1 = imdilate(h1, strel('disk', 10));
                %                 h2 = ~m2 & (r1 | r2);
                h2 = ~m2 & r2;
                hh = h1 + h2;
                [l, n] = bwlabeln(h2);
                s = zeros(n, 2);
                for i = 1: n
                    tmpt = l == i;
                    tmp = hh(tmpt);
                    nn = max(tmp(:));
                    [y, x] = find(tmpt);
                    s(i, 2: 3) = [mean(y), length(y)];
                    if nn == 1
                        s(i, 1) = 1;
                    else
                        st = sum(sum(tmp > 1)) / s(i, 3);
                        if st < stthres
                            s(i, 1) = 1;
                        end
                    end
                end
                id1 = s(:, 3) > nthres;
                id2 = s(:, 2) < sthres;
                id1t = find(s(:, 1) & id1 & id2);
                [~, id] = max(s(id1t, 2));
                id = id1t(id);
                if ~isempty(id)
                    mask2 = l == id;
                else
                    mask2 = isnan(l);
                end
                
                %%% refine mask %%%
                maskn = cc.mask_refine(mask2);
                m = max(maskn(:));
                s = zeros(m, 1);
                for i = 1: m
                    tmpt = maskn == i;
                    tmp = hh(tmpt);
                    nn = max(tmp(:));
                    if nn == 1
                        s(i) = 1;
                    end
                end
                idt = find(s > 0);
                s = zeros(length(idt), 1);
                for i = 1: length(idt)
                    tmp = maskn == idt(i);
                    s(i) = sum(tmp(:));
                end
                [~, idtt] = max(s);
                mask2 = maskn == idt(idtt);
                mask = mask | mask2;
                
                %%% final output data %%%
                [h, w] = find(mask);
                pll = [h + rlocn(1), w + rlocn(3)];
                masks{ii} = mask;
                plls{ii} = pll;
            end
        end
        
        function [model, pts, pti] = lens_recon(cc, lensdata)
            %%% has to have a series of slices with lens lesion %%%
            % get the major blob %
            img = lensdata.mask;
            idp = cellfun(@(x) find(x == true), img, 'uniformoutput', false);
            img = reshape(cell2mat(img), size(img{1}, 1), size(img{1}, 2), []);
            [l, n] = bwlabeln(img);
            s = zeros(n, 1);
            for i = 1: n
                tmp = l == i;
                s(i) = sum(tmp(:));
            end
            [~, idt] = max(s);
            imgot = l == idt;
            imgo = reshape(imgot, size(imgot, 1), []);
            imgo = mat2cell(imgo, size(imgo, 1), repmat(size(imgot, 2), 1, size(imgot, 3)));
            idpo = cellfun(@(x) find(x == true), imgo, 'uniformoutput', false);
            idf = cell(1, length(idpo));
            for i = 1: length(idpo)
                tmpa = idp{i};
                tmpb = idpo{i};
                [~, tmp] = intersect(tmpa, tmpb, 'stable');
                idf{i} = tmp;
            end
            
            pll3 = lensdata.point_lists_3d;
            pll3f = cellfun(@(x, y) x(y, :), pll3, idf, 'uniformoutput', false);
            
            t = cell2mat(pll3f');
            tt = t;
            tt(:, 3) = t(:, 3) + mean(diff(cc.metadata.angles.ap)) * (rand(size(t(:, 3))) - 0.5);
            maxDistance = 1;
            pt = pointCloud(tt);
            [model, inlierIndices] = pcfitcylinder(pt, maxDistance);
            pti = select(pt, inlierIndices);
            pts = pt;
        end
        
        function maskn = mask_refine(mask)
            t = bwdist(~mask);
            tt = imhmin(-t, prctile(t(t > 0), 50));
            maskn = watershed(tt);
            maskn(~mask) = 0;
        end
        
        %% aux functions: rendering %%
        function [bp, mask] = prep_renderer(cc, tv)
            %%% get 3d brain for rendering %%%
            if exist('ref', 'var')
                bp = ref;
            else
                [~, bp] = prep_ccf(cc, tv, 0.25);
            end
            bp = single(bp);
            bp = normalize(bp);
            mask = bp > 0.1;
            for i = 1: size(mask, 3)
                mask(:, :, i) = imfill(mask(:, :, i), 'holes');
            end
            bp(~mask) = 0;
        end
            
        function [bp, mx] = merge_cells(cc, bp, mask)
            %%% add points &/or lens %%%
            pl3 = cc.metadata.point_lists_3d;
            scl = size(bp, 1) / 800;
            [h, w, d] = size(bp);
            stp = 0;
            rds = 16;
            mx = 5;
            for i = 1: length(pl3)
                pt = pl3{i};
                ptr = pt;
                ptr(:, 1: 2) = round(pt(:, 1: 2) * scl);
                ptr(:, 3) = round(pt(:, 3)) + round(rand(size(ptr(:, 3))) * rds - rds / 2);
                for j = 1: size(ptr, 1)
                    tmp = mask(ptr(j, 1), ptr(j, 2), ptr(j, 3));
                    if tmp
                        bp = cc.fill_3d_neighbor(bp, ptr(j, 2), ptr(j, 1), ptr(j, 3), mx);
%                         bp(max(1, ptr(j, 2) - stp): min(h, ptr(j, 2) + stp), max(1, ptr(j, 1) - stp): min(w, ptr(j, 1) + stp), max(1, ptr(j, 3) - stp): min(d, ptr(j, 3) + stp)) = 10;
                    end
                end
            end
            bp = normalize(bp);
        end
        
        function bp = fill_3d_neighbor(bp, x, y, z, vl)
            bp(x, y, z) = vl;
%             bp(x - 1, y, z) = vl;
%             bp(x + 1, y, z) = vl;
%             bp(x, y - 1, z) = vl;
%             bp(x, y + 1, z) = vl;
            bp(x, y, z - 1) = vl;
            bp(x, y, z + 1) = vl;
            bp(x, y, z - 2) = vl;
            bp(x, y, z + 2) = vl;
        end
        
        function bpl = create_cyl_label(bp, model, scl)
            [nh, nw, nn] = size(bp);
%             sclt = [scl, scl, 1];
            A = model.Orientation;
            B = model.Parameters(1: 3);
            C = model.Parameters(4: 6);
            ctr = model.Center;
            bmx = max(model.Height / 2, model.Radius);
%             [X1, X2, X3] = ndgrid(1: round(nw / scl), 1: round(nh / scl), 1: nn);
            [X1, X2, X3] = ndgrid(max(1, round(ctr(1) - bmx)): min(round(nw / scl), round(ctr(1) + bmx)), max(1, round(ctr(2) - bmx)): min(round(nh / scl), round(ctr(2) + bmx)), max(1, round(ctr(3) - 1.5 * bmx)): min(nn, round(ctr(3) + bmx)));
            X1 = single(X1);
            X2 = single(X2);
            X3 = single(X3);
            
            %%% coarse selection %%%
            tp = [X1(:), X2(:), X3(:)];
            id2 = (tp - 1.5 * B) * A' >= 0;
            id3 = (tp - C) * (-A') >= 0;
            tp = tp(id2 & id3, :);
            
            %%% exact selection %%%
            %%%% point-line distance %%%%
            id1 = vecnorm(cross(tp - B, repmat(A, size(tp, 1), 1), 2), 2, 2) / norm(A) <= model.Radius;
            tpf = tp(id1, :);
            tpf(:, 1: 2) = round(tpf(:, 1: 2) * scl);
            tpf = unique(tpf, 'rows');
            
            %%% create label %%%
            bpl = zeros(size(bp));
            for i = 1: size(tpf, 1)
                bpl(tpf(i, 2), tpf(i, 1), tpf(i, 3)) = 1;
            end
        end
        
        %% aux functions: utility %%
        function [cs, ss, mui, rloc, ks] = mea_cal_unit(cc, ref, opt)
            bf = normalize(opt.bf);
            rscl = opt.rscl;
            ref_flag = opt.ref_flag;
            [hh, ww] = size(bf);
            r = ref;
            r = normalize(imresize(r, rscl));
            if opt.cfflag == 1
                rf = normalize(imgaussfilt(TVL1denoise(normalize(r .^ 0.2), 0.8, 40), max(0.5, floor(min(size(r)) / 200))));
            elseif opt.cfflag == 2
                rf = normalize(imgaussfilt(TVL1denoise(normalize(r .^ 0.2), 0.8, 40), max(0.5, floor(min(size(r)) / 200))));
            else
                rf = normalize(imgaussfilt(TVL1denoise(normalize(r .^ 0.2), 0.8, 40), max(0.5, floor(min(size(r)) / 200)))); 
            end
            
            if isfield(opt, 'loc')
                idxo = opt.loc(3) + ww - 1;
                idyo = opt.loc(1) + hh - 1;
            else
                t = normxcorr2(bf, rf);
                [idxo, idyo, scr] = cc.ref_locate(t, ref_flag);
            end
                        
            if ~isempty(idxo)
                hmin = idyo - hh;
                hmax = hmin + hh - 1;
                wmin = idxo - ww;
                wmax = wmin + ww - 1;
                try
                    rr = normalize(rf(hmin: hmax, wmin: wmax));
                    if opt.cfflag == 2
                        [optimizer, metric] = imregconfig('multimodal');
                        optimizer.InitialRadius = 0.0001;
                        optimizer.Epsilon = 1.5e-4;
                        optimizer.GrowthFactor = 1.02;
                        optimizer.MaximumIterations =100;
                        imgt = imregister(bf, rr, 'affine', optimizer, metric);
                    elseif opt.cfflag == 1
                        imgt = normalize(bf);
                        scrt = 100;
                    elseif opt.cfflag == 3
                        rt = normalize(r(hmin: hmax, wmin: wmax));
                        rr = cc.opt_exp(cc, rt, bf, 0.8);
            
                        [optimizer, metric] = imregconfig('multimodal');
                        optimizer.InitialRadius = 0.0001;
                        optimizer.Epsilon = 1.5e-4;
                        optimizer.GrowthFactor = 1.02;
                        optimizer.MaximumIterations =100;
                        imgt = imregister(bf, rr, 'affine', optimizer, metric);
                    else
                        [scrt, imgt] = cc.klt_ref_track(cc, bf, rr);
                    end
                    
                    imgt = normalize(imgt);
                    if opt.cfflag > 1
                        ks = cc.hog_sim(imgt, rr);
                    else
                        ks = NaN;
                    end
                    
                    stmp = multissim(imgt, rr);
                    mtmp = mi(imgt, rr);
                    cs = cc.cor_sim(imgt, rr);
                    ss = stmp;
                    mui = mtmp;
                    rloc = [hmin, hmax, wmin, wmax];
                catch
                    cs = 0;
                    ss = 0;
                    mui = 0;
                    rloc = NaN(1, 4);
                    ks = -200;
                end
            else
                cs = 0;
                ss = 0;
                mui = 0;
                rloc = NaN(1, 4);
                ks = -200;
            end
        end
        
        function ks = hog_sim(imgt, rr)
            t1 = extractHOGFeatures(imgt, 'cellsize', [2, 2]);
            t2 = extractHOGFeatures(rr, 'cellsize', [2, 2]);
            ks = t1 * t2' / (norm(t1) * norm(t2));
        end
        
        function cs = cor_sim(imgt, rr)
            cs = corrcoef(imgt(:), rr(:));
            cs = cs(1, 2);
        end
        
        function [idxo, idyo, scr] = ref_locate(cor, flag)
            [yt, xt] = size(cor);
            tt = imregionalmax(cor);
            st = cor(tt);
            id = find(tt(:));
            [idy, idx] = ind2sub(size(cor), id);
            [sts, idt] = sort(st, 'descend');
            idy = idy(idt);
            idx = idx(idt);
            switch flag
                case 1 %%% symmetric to y axis %%%
                    idd = find(idx > round(xt / 2) - 5 & idx < round(xt / 2) + 5, 1);
                    if ~isempty(idd)
                        idxo = idx(idd);
                        idyo = idy(idd);
                        scr = sts(idd);
                    else
                        idxo = [];
                        idyo = [];
                        scr = [];
                    end
                case 2 %%% partial left/right slice %%%
                    idd = find(idx < round(xt / 4) & idx > round(xt / 4), 1);
                    if ~isempty(idd)
                        idxo = idx(idd);
                        idyo = idy(idd);
                        scr = sts(idd);
                    else
                        idxo = [];
                        idyo = [];
                        scr = [];
                    end
                case 3
                    idd = idx;
                    if ~isempty(idd)
                        idxo = idx(idd);
                        idyo = idy(idd);
                        scr = sts(idd);
                    else
                        idxo = [];
                        idyo = [];
                        scr = [];
                    end
            end
        end
        
        function coord = prep_coord(nh, nw)
            [X, Y] = meshgrid(1: nw, 1: nh);
            Xt = X - (1 + nw) / 2;
            Yt = Y - (1 + nh) / 2;
            coord = [Xt(:), Yt(:), zeros(size(Xt(:)))];
        end
        
        function imgo = opt_exp(cc, imcur, imreft, lamb, iter, sz)
            if nargin < 4 || isempty(lamb)
                lamb = 0.2;
            end
            
            if nargin < 5 || isempty(iter)
                iter = 40;
            end
            
            if nargin < 6 || isempty(sz)
                sz = [200, 285];
            end
            
            stp = 5;
            ss = zeros(stp, 1);
            cs = ss;
            imcurtmp = imresize(imcur, sz);
            imrefttmp = imresize(imreft, sz);
            for i = 1: stp
                imcurt = normalize(imgaussfilt(normalize(TVL1denoise(normalize(imcurtmp), lamb, iter)) .^ (1 + 0.2 * i), max(1, floor(min(size(imreft)) / 200))));
                ss(i) = ssim(imcurt, imrefttmp);
                cs(i) = cc.cor_sim(imcurt, imrefttmp);
            end
            
            tmp = ss + cs;
            [~, id] = max(tmp);
            imgo = normalize(imgaussfilt(normalize(TVL1denoise(normalize(imcur), lamb, iter)) .^ (1 + 0.2 * id), max(1, floor(min(size(imreft)) / 200))));
        end
        
        function [img, xform] = klt2_reg(cc, imref, imcur, flag, maskc, pps)
            if nargin < 4
                flag = 1;
            end
            
            if nargin < 5 || isempty(maskc)
                maskc = true(size(imref));
            end
            
            if nargin < 6 || isempty(pps)
                [~, pcur, pold, ~] = cc.klt2(imref, imcur, [], [], [], [], maskc);
            else
                pcur = pps.pcur;
                pold = pps.pold;
            end
            
            [~, ~, xform] = cc.klt_geo(pold, pcur, [], flag); %%% 1: regular RANSAC trans; else: fitgeotrans %%%
            if ~isempty(xform)
                img = cc.klt_warp(imcur, xform);
            else
                img = imcur;
                xform = affine2d(diag(ones(1, 3)));
            end
        end
        
        function [inliercur, inlierold, xform] = klt_geo(pold, pcur, hthres, flag)
            if nargin < 3 || isempty(hthres)
                hthres = 1.5; %%% matlab default %%%
            end

            if nargin < 4
                flag = 1;
            end

            try
                if flag == 1
                    [xform, inliercur, inlierold] = ...
                        estimateGeometricTransform(pcur, pold, ...
                        'affine', 'MaxDistance', hthres);
                else
                    xform = fitgeotrans(pcur, pold, 'affine');
                    inliercur = pcur;
                    inlierold = pold;
                end
            catch
                xform = [];
                inliercur = [];
                inlierold = [];
            end
        end

        function [d, TPoints, SInliers, isFound] = klt2(S, T, biderr, mq, winsize, oldpoints, maskc)
            if nargin < 3 || isempty(biderr)
                biderr = 2;
            end
            if nargin < 4 || isempty(mq)
                mq = 0.02;
            end
            if nargin < 5 || isempty(winsize)
                winsize = 2 * [round(min(size(S)) / 4), round(min(size(S)) / 4)] + 1;
            end
            if nargin < 6 || isempty(oldpoints)
                oldpoints = detectMinEigenFeatures(S, 'MinQuality', mq);
                oldpoints = oldpoints.Location;
            end
            flagmask = true;
            if nargin < 7 || isempty(maskc)
                flagmask = false;
            end
            pointTracker = vision.PointTracker('MaxBidirectionalError', biderr, 'NumPyramidLevels', 4, 'BlockSize', winsize);
            initialize(pointTracker, oldpoints, S);
            [points, isFound] = step(pointTracker, T);
            TPoints = points(isFound, :);
            SInliers = oldpoints(isFound, :);
            if flagmask
                idS = [];
                for i = 1: size(SInliers, 1)
                    if maskc(round(SInliers(i, 2)), round(SInliers(i, 1)))
                        idS = [idS, i];
                    end
                end
                TPoints = TPoints(idS, :);
                SInliers = SInliers(idS, :);
            end
            d = TPoints - SInliers;
        end
        
        function img_warp = klt_warp(img, xform)
            img_warp = imwarp(img, xform, 'OutputView', imref2d(size(img)));
        end
        
        function [scrout, imgout, xformout] = klt_ref_track(cc, Y, imref, nmaxloop, maskc)
            if nargin < 4 || isempty(nmaxloop)
                nmaxloop = 10;
            end
            
            if nargin < 5 || isempty(maskc)
                maskc = true(size(imref));
            end
            
            %%% KLT tracker %%
            %%% simple version: no loop just once %%%
            nf = size(Y, 3);
            scrout = zeros(1, nf);
            xformout = cell(1, nf);
            imgout = zeros(size(Y));
            
            %%% loop through all neighboring frame pairs %%%
            for i = 1: nf
                imcur = Y(:, :, i);
                
                xforms = cell(size(imref, 3), nmaxloop + 1);
                imgs = cell(size(imref, 3), nmaxloop + 1);
                scrs = zeros(size(imref, 3), nmaxloop + 2); %%% additional score for untransformed %%%
                
                %%% nmaxloop times realizations %%%
                [~, pcur, pold, ~] = cc.klt2(imref, imcur, [], [], [], [], maskc);
                pps.pcur = pcur;
                pps.pold = pold;
                if nmaxloop == 1
                    gflag = 2;
                else
                    gflag = 1;
                end
                
                for ii = 1: nmaxloop
                    [img, xf] = cc.klt2_reg(cc, imref, imcur, gflag, maskc, pps);
                    scr = cc.get_trans_score_ref(cc, img, imref, maskc);
                    xforms{ii} = xf;
                    imgs{ii} = img;
                    scrs(ii) = scr;
                end
                
                %%% gather additional fitgeotrans and scores of unregistered frame pairs %%%
                [img, xf] = klt2_reg(imref, imcur, 2, maskc);
                scr = cc.get_trans_score_ref(cc, img, imref, maskc);
                xforms{nmaxloop + 1} = xf;
                imgs{nmaxloop + 1} = img;
                scrs(nmaxloop + 1) = scr;
                scr = cc.get_trans_score_ref(cc, img, imref, maskc);
                scrs(nmaxloop + 2) = scr;
                
                %%% find the transform to use %%%
                scr = min(scrs(:));
                [x, y] = find(scrs == scr);
                if y(1) < nmaxloop + 2
                    img = imgs{x(1), y(1)};
                    xform = xforms{x(1), y(1)};
                else
                    img = imcur;
                    xform = affine2d(diag(ones(1, 3)));
                end
                
                scrout(i) = scr;
                imgout(:, :, i) = img;
                xformout{i} = xform;
            end
        end
        
        function acorr = get_trans_score_ref(cc, Y, imref, maskc)
            if nargin < 4 || isempty(maskc)
                maskc = true(size(imref));
            end
            
            mq = 0.01;
            biderr = 2;
            [~, ~, nframes] = size(Y);
            acorr = zeros(1, nframes);
            for i = 1: nframes
                img_old = imref;
                img = Y(:, :, i);
                d = cc.klt2(normalize(img_old), normalize(img), biderr, mq, [], [], maskc);
                if ~isempty(d)
                    temp = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
                    acorr(i) = temp;
                else
                    acorr(i) = 100; %%% a large score %%%
                end
            end
        end
        
        function mxpsig = feature2_comp(cc, mxall, flag, mag, denom, sz)
            % compute 2nd features of the images
            %   Jinghao Lu, 01/26/2018
            
            if nargin < 3 || isempty(flag)
                flag = 0;
            end
            
            if nargin < 4 || isempty(mag)
                mag = 40;
            end
            
            if nargin < 5 || isempty(denom)
                denom = 10;
            end
            
            if nargin < 6 || isempty(sz)
                sz = 9;
            end
            
            mxp = mxall;
            n = size(mxall, 3);
            
            %%% if do anidenoise %%%
            if flag == 1
                parfor i = 1: n
                    tmp = anidenoise(imerode(mxp(:, :, i), strel('disk', 1)), sz);
                    mxp(:, :, i) = tmp(sz + 1: end - sz, sz + 1: end - sz);
                end
            end
            
            %%% 2nd feature, demons_prep %%%
            mxpsig = mxp; %%% mxpsig: just slightly enhance %%%
            for i = 1: n
                basecur = mxpsig(:, :, i);
                mxpsig(:, :, i) = cc.demons_prep(basecur, mag, denom);
            end
            mxpsig = bsxfun(@minus, mxpsig, min(min(mxpsig, [], 1), [], 2));
        end
        
        function imout = demons_prep(img, mag, denom)
            % create feature maps good for demons
            %   Jinghao Lu, 05/11/2017
            
            if nargin < 2
                mag = 20;
            end
            if nargin < 3
                denom = 3;
            end
            ithres = (max(img(:)) + min(img(:))) ./ denom;
            tmp = 1./(1 + exp(-mag * (img - ithres)));
            tmp = tmp .* (tmp >= median(tmp(:)));
            tmp(tmp == 0) = min(tmp(tmp > 0));
            imout = normalize(tmp);
        end
        
        function I = iminterpolate(I, sx, sy)
            % image interpolate for LogDemons
            %   Jinghao Lu, 06/11/2016
            
            %%% Find update points on moving image %%%
            [x, y] = ndgrid(0: (size(I, 1) - 1), 0: (size(I, 2) - 1)); % coordinate image
            x_prime = x + sx; % updated x values (1st dim, rows)
            y_prime = y + sy; % updated y values (2nd dim, cols)
            
            %%% Interpolate updated image %%%
            I = interpn(x, y, I, x_prime, y_prime, 'linear', 0); % moving image intensities at updated points
        end
        
        function img = multi_otsu(cc, ain)
            nl = 20;
            s = zeros(1, nl);
            level = multithresh(ain, nl);
            for i = 1: nl
                tmp = ain > level(i);
                [l, n] = bwlabeln(tmp);
                s(i) = n;
            end
            
            if (max(s) - min(s)) / max(s) > 0.1
                lvl = interp1(1: length(level), level, 1: 0.02: length(level));
                ss = interp1(1: length(s), s, 1: 0.02: length(s));
                ss = sgolayfilt(ss, 1, 2 * round(0.1 * length(ss)) + 1);
                
                thres = cc.intense_filter(ss(:), 1);
                id = find(ss < thres, 1);
                img = ain > lvl(id);
            else
                img = ain;
            end
        end
        
        function ithres = intense_filter(maxall, scl)
            if nargin < 2 || isempty(scl)
                scl = 0.2;
            end
            tmp1 = maxall;
            tmp1 = tmp1(tmp1 > 1/2^8); %%% above 0 in terms of single precision %%%
            tmp11 = tmp1;
            if ~isempty(tmp1)
                if ~isvector(maxall)
                    tmp11 = sort(tmp11);
                end
                tmp2 = linspace(tmp11(1), 1 * tmp11(end), length(tmp1));
                tmp = tmp11(:) - scl * tmp2(:);
                if sum(tmp) > 0
                    tmp11 = max(tmp11) - tmp11;
                end
                tmp2 = linspace(tmp11(1), 1 * tmp11(end), length(tmp1));
                tmp = tmp11(:) - scl * tmp2(:);
                x = 1: length(tmp);
                sker = 2 * round(length(tmp) / 100) + 1;
                xq = [1 - sker: 0, x, length(tmp) + 1: length(tmp) + sker];
                tmpt = interp1(x, tmp, xq, 'linear', 'extrap');
                tmpg = smooth(diff(smooth(tmpt, sker)), sker);
                tmpg = tmpg(x + sker + 1);
                idthres = find(tmpg >= 0, 1);
                ithres = min(prctile(tmp1, 90), tmp1(idthres));
            else
                ithres = 0;
            end
        end
        
        function YDeN = anidenoise(cc, Y, sz, ispara, iter, dt, kappa, opt)
            % Yblur = anidenoise(Y) denoise using anisotropic diffusion
            %   denoise image/sequence
            %   Y is input image/sequence
            %   Yblur is output image/sequence
            %   Jinghao Lu 05/17/2016
            
            %%% initialization %%%
            %%% initialize parameters %%%
            if nargin < 3 || isempty(sz)
                sz = 5;
            end
            
            if nargin < 4 || isempty(ispara)
                ispara = 0;
            end
            
            if nargin < 5 || isempty(iter)
                iter = 4;
            end
            
            if nargin < 6 || isempty(dt)
                dt = 0.1429; %%% maximum for numerical stability, and reduce iterations %%%
            end
            
            if nargin < 7 || isempty(kappa)
                kappa = 0.5; %%% any value above 0.1 for normalized image %%%
            end
            
            if nargin < 8 || isempty(opt)
                opt = 1;
            end
            
            %%% prepare data %%%
            Y = padarray(Y, [sz, sz], 'replicate');
            nframes = size(Y, 3);
            YDeN = zeros(size(Y), class(Y));
            
            %%% begin anisotropic diffusion %%%
            if nframes == 1
                YDeN = cc.anisodiff_unit(Y, iter, dt, kappa, opt);
            else
                if ispara
                    if isempty(gcp('nocreate'))
                        parpool(feature('numCores'));
                    end
                    parfor i = 1: nframes
                        I = Y(:, :, i);
                        YDeN(:, :, i) = anisodiff_unit(I, iter, dt, kappa, opt);
                        if mod(i, 100) == 0
                            disp(['done frame #', num2str(i), '/', num2str(nframes)])
                        end
                    end
                else
                    for i = 1: nframes
                        I = Y(:, :, i);
                        YDeN(:, :, i) = anisodiff_unit(I, iter, dt, kappa, opt);
                        if mod(i, 100) == 0
                            disp(['done frame #', num2str(i), '/', num2str(nframes)])
                        end
                    end
                end
            end
        end
        
        function diff_im = anisodiff_unit(im, num_iter, delta_t, kappa, option)
            % diff_im = anisodiff_unit(im, num_iter, delta_t, kappa, option)
            %   conventional anisotropic diffusion (Perona & Malik) on a gray scale
            %   image. A 2D network structure of 8 neighboring nodes is considered for
            %   diffusion conduction.
            %
            %       ARGUMENT DESCRIPTION:
            %              IM       - gray scale image (MxN).
            %              NUM_ITER - number of iterations.
            %              DELTA_T  - integration constant (0 <= delta_t <= 1/7).
            %                         Usually, due to numerical stability this
            %                         parameter is set to its maximum value.
            %              KAPPA    - gradient modulus threshold that controls the conduction.
            %              OPTION   - conduction coefficient functions proposed by Perona & Malik:
            %                         1 - c(x,y,t) = exp(-(nablaI/kappa).^2),
            %                             privileges high-contrast edges over low-contrast ones.
            %                         2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2),
            %                             privileges wide regions over smaller ones.
            %
            %       OUTPUT DESCRIPTION:
            %               DIFF_IM - (diffused) image with the largest scale-space parameter.
            %
            %   Example
            %   -------------
            %   s = phantom(512) + randn(512);
            %   num_iter = 15;
            %   delta_t = 1/7;
            %   kappa = 30;
            %   option = 2;
            %   ad = anisodiff2D(s,num_iter,delta_t,kappa,option);
            %   figure, subplot 121, imshow(s,[]), subplot 122, imshow(ad,[])
            %
            % References:
            %   P. Perona and J. Malik.
            %   Scale-Space and Edge Detection Using Anisotropic Diffusion.
            %   IEEE Transactions on Pattern Analysis and Machine Intelligence,
            %   12(7):629-639, July 1990.
            %
            %   G. Grieg, O. Kubler, R. Kikinis, and F. A. Jolesz.
            %   Nonlinear Anisotropic Filtering of MRI Data.
            %   IEEE Transactions on Medical Imaging,
            %   11(2):221-232, June 1992.
            %
            %   Codes adapted from Daniel Simoes Lopes
            %   ICIST
            %   Instituto Superior Tecnico - Universidade Tecnica de Lisboa
            %   danlopes (at) civil ist utl pt
            %   http://www.civil.ist.utl.pt/~danlopes
            %
            %   Jinghao Lu 01/06/2016
            
            % Convert input image to double.
            im = double(im);
            
            % PDE (partial differential equation) initial condition.
            diff_im = im;
            
            % Center pixel distances.
            dx = 1;
            dy = 1;
            dd = sqrt(2);
            
            % 2D convolution masks - finite differences.
            hN = [0 1 0; 0 -1 0; 0 0 0];
            hS = [0 0 0; 0 -1 0; 0 1 0];
            hE = [0 0 0; 0 -1 1; 0 0 0];
            hW = [0 0 0; 1 -1 0; 0 0 0];
            hNE = [0 0 1; 0 -1 0; 0 0 0];
            hSE = [0 0 0; 0 -1 0; 0 0 1];
            hSW = [0 0 0; 0 -1 0; 1 0 0];
            hNW = [1 0 0; 0 -1 0; 0 0 0];
            
            % Anisotropic diffusion.
            for t = 1:num_iter
                
                % Finite differences. [imfilter(.,.,'conv') can be replaced by conv2(.,.,'same')]
                %         nablaN = imfilter(diff_im,hN,'conv');
                %         nablaS = imfilter(diff_im,hS,'conv');
                %         nablaW = imfilter(diff_im,hW,'conv');
                %         nablaE = imfilter(diff_im,hE,'conv');
                %         nablaNE = imfilter(diff_im,hNE,'conv');
                %         nablaSE = imfilter(diff_im,hSE,'conv');
                %         nablaSW = imfilter(diff_im,hSW,'conv');
                %         nablaNW = imfilter(diff_im,hNW,'conv');
                nablaN = conv2(diff_im, hN, 'same');
                nablaS = conv2(diff_im, hS, 'same');
                nablaW = conv2(diff_im, hW, 'same');
                nablaE = conv2(diff_im, hE, 'same');
                nablaNE = conv2(diff_im, hNE, 'same');
                nablaSE = conv2(diff_im, hSE, 'same');
                nablaSW = conv2(diff_im, hSW,'same');
                nablaNW = conv2(diff_im, hNW, 'same');
                
                % Diffusion function.
                if option == 1
                    cN = exp(-(nablaN/kappa).^2);
                    cS = exp(-(nablaS/kappa).^2);
                    cW = exp(-(nablaW/kappa).^2);
                    cE = exp(-(nablaE/kappa).^2);
                    cNE = exp(-(nablaNE/kappa).^2);
                    cSE = exp(-(nablaSE/kappa).^2);
                    cSW = exp(-(nablaSW/kappa).^2);
                    cNW = exp(-(nablaNW/kappa).^2);
                elseif option == 2
                    cN = 1./(1 + (nablaN/kappa).^2);
                    cS = 1./(1 + (nablaS/kappa).^2);
                    cW = 1./(1 + (nablaW/kappa).^2);
                    cE = 1./(1 + (nablaE/kappa).^2);
                    cNE = 1./(1 + (nablaNE/kappa).^2);
                    cSE = 1./(1 + (nablaSE/kappa).^2);
                    cSW = 1./(1 + (nablaSW/kappa).^2);
                    cNW = 1./(1 + (nablaNW/kappa).^2);
                end
                
                % Discrete PDE solution.
                diff_im = diff_im + ...
                    delta_t*(...
                    (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
                    (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE + ...
                    (1/(dd^2))*cNE.*nablaNE + (1/(dd^2))*cSE.*nablaSE + ...
                    (1/(dd^2))*cSW.*nablaSW + (1/(dd^2))*cNW.*nablaNW );
                %                   delta_t*(...
                %                   (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
                %                   (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE);
                
                % Iteration warning.
                %         fprintf('\rIteration %d\n',t);
            end
        end
        
        %% aux functions: visualization %%
        function tp = image_fuse(imref, img, rloc)
            if nargin < 3 || isempty(rloc)
                rloc = [1; size(img, 1); 1; size(img, 2)];
            end
            tp = zeros([size(imref), size(img, 3)]);
            tp(rloc(1): rloc(2), rloc(3): rloc(4), :) = img;
            tp =  imfuse(imref * 4, tp * 4, 'blend', 'Scaling', 'none');
        end
        
        function cell_overlay(img, pl)
%             figure(gcf)
%             clf
            imagesc(img)
            hold on
            for i = 1: size(pl, 1)
                plot(pl(i, 2), pl(i, 1), '.r')
            end
            hold off
        end

        function tp = lens_overlay(img, mask, rloc)
            if nargin < 3 || isempty(rloc)
                rloc = [1; size(img, 1); 1; size(img, 2)];
            end
            tpt = zeros([size(img(:, :, 1))]);
            tpt(rloc(1): rloc(2), rloc(3): rloc(4), :) = mask;
            tp =  imfuse(tpt, img, 'blend', 'Scaling', 'none');
        end
        
        function vol = plot_3d_brain(bpm, mx)
            ns = 256;
            pmt = [1, 2, 3];
            intensity = [0, mx / 200, 1];
            alphat = [0.4, mx / 8, 0.8];
            colort = [1, 1, 1; 0.8, 0.8, 0.8; 0, 1, 0];
            qp = linspace(0, 1, ns);
            alpham = interp1(intensity, alphat, qp)';
            colorm = interp1(intensity, colort, qp);
            set(gcf, 'renderer', 'opengl')
            vol = volshow(permute(bpm, pmt), Colormap = colorm, Alphamap = alpham);
            vol.Parent.BackgroundColor = 'k';
            vol.Parent.BackgroundGradient = 'off';
            tmp = [4, 4, 1];
            vol.Transformation.A(1: 3, 1: 3) = diag(tmp);
            vol.Parent.CameraPosition = [1100, 400, 670];
            vol.Parent.CameraUpVector = [0, -1, 0];
            vol.RenderingStyle = 'GradientOpacity';
        end
        
        function plot_3d_brain_lens(cc, bp, bpl, mx, vol)
            if nargin < 5 || isempty(vol)
                vv = cc.plot_3d_brain(bp, mx);
            else
                vv = vol;
            end
            vv.OverlayData = bpl;
            vv.OverlayAlphamap = 1;
        end

        function brain_animation(vol, fn, pos, upv, zoom)
            if nargin < 3 || isempty(pos)
                t = linspace(0, 2 * pi, 200)';
                pos = [cos(-t) * 1100, zeros(200, 1), sin(-t) * 1100] + [425, 400, 670];
            end
            if nargin < 4 || isempty(upv)
                upv = repmat([0, -1, 0], 1000, 1);
            end
            if nargin < 5 || isempty(zoom)
                zoom = ones(1, 1000);
            end

            id = min([size(pos, 1), size(upv, 1), length(zoom)]);
            pos = pos(1: id, :);
            upv = upv(1: id, :);
            zoom = zoom(1: id);

%             v = VideoWriter(fn, 'MPEG-4');
%             v = VideoWriter(fn, 'uncompressed avi');
%             v.Quality = 100;
%             v.FrameRate = 30;
%             v.LosslessCompression = true;
%             open(v)
            for i = 1: id
                vol.Parent.CameraPosition = pos(i, :);
                vol.Parent.CameraUpVector = upv(i, :);
                vol.Parent.CameraZoom = zoom(i);

                I = getframe(vol.Parent.Parent);
%                 writeVideo(v, I);
                [indI, cm] = rgb2ind(I.cdata, 256);
                if i == 1
                    imwrite(indI, cm, fn, 'gif', Loopcount = inf, DelayTime = 0.025)
                else
                    imwrite(indI, cm, fn, 'gif', WriteMode = 'append', DelayTime = 0.025)
                end
            end
%             close(v)
        end
        
        function pc_lens_model(pts, pti, model)
            pcshow(pts.Location, 'w')
            hold on
            pcshow(pti.Location, 'r', 'markersize', 10)
            plot(model)
        end
    end
end




