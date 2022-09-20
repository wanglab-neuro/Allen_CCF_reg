classdef cellcount
    properties
        ccf = [];
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
            cc.pname = path_name;
            
            %%%% if you want a pop up window, use this: %%%%
            nlist = cc.section_list(length(file_name));
            
            %%%% if not, put numbers in here (and don't run the above line): %%%%
            nlist = [];
            
            %%% initialize %%%
            cc = cellcount();
            [tv, av, st] = cc.load_ccf();
            [cc, ref] = cc.prep_ccf(tv);
            
            %%% main loop per slice %%%
            for i = 1: length(file_name)
                %%% load slice %%%
                fname = file_name{i};
                slo = cc.load_slice(fname);
                
                %%% pyramid size selection %%%
                sl = cc.prep_slice(slo);
                ref_flag = 1;
                n = nlist(i);
                [gpy, loc, rsclp] = cc.find_loc_pyramid(cc, ref, sl(:, :, 3), ref_flag, n);
                
                %%% coarse selection %%%
                idf = cc.find_loc_coarse(cc, ref, sl(:, :, 3), ref_flag, n, gpy, loc, rsclp);
                
                %%% fine/rotate selection %%%
                rtl = size(tv, 2) / size(ref, 1);
                [theta, gamma, ap, rloc, imref] = cc.find_loc_fine(cc, ref, sl(:, :, 3), ref_flag, idf, gpy, loc, rtl, rsclp);
                
                %%% nonrigid register %%%
                angles = struct('theta', theta, 'gamma', gamma, 'ap', ap + idf);
                [img, D, rlocn, imr, coordf] = cc.regist(cc, tv, sl, angles, rloc, imref);
                
                %%% Detect labeled neurons %%%
                opt.vis = false;
                %                 opt.vis = true;
                pl = cc.cell_detect(cc, sl(:, :, 2), opt);
                
                %%% transform detected neurons %%%
                opt.stv = size(img);
                opt.ssl = size(sl);
                pln = cc.cell_transform(D, pl, rlocn, opt);
                
                %%% get 3d final point list %%%
                pl3 = cc.cell3d(pln, coordf);
                
                if lens_flag
                    %%% find lens lesion %%%
                    [mask, pll] = cc.lens_loc(cc, img, imr, rlocn);
                    pll3 = cc.cell3d(pll, coordf);
                else
                    pll = NaN;
                    pll3 = NaN;
                end
                
                %%% metadata compile %%%
                cc.pyramid_size{i} = gpy;
                cc.angles{i} = angles;
                cc.point_lists{i} = pl;
                cc.point_lists_new{i} = pln;
                cc.point_lists_3d{i} = pl3;
                cc.lens_point_lists{i} = pll;
                cc.lens_point_lists_3d{i} = pll3;
                cc.section_location{i} = rlocn;
                cc.reg_image{i} = img;
                cc.displacement{i} = D;
                cc.reg_ref_image{i} = imr;
                cc.coordf{i} = coordf;
                cc.fnames{i} = fname;
                
                %%% save data %%%
                cc.save(cc);
            end

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
        
        function img = load_slice(fname)
            img = imread(fname);
        end
        
        function imgo = prep_slice(img)
            a = normalize(double(img));
            b = normalize(a(:, :, 3));
            b1 = imgaussfilt(b, floor(min(size(b)) / 100));
            b2 = sum(b1 > 0.1, 1);
            b3 = sum(b1 > 0.1, 2);
            bb = a(find(b3, 1): find(b3, 1, 'last'), find(b2, 1): find(b2, 1, 'last'), :);
            imgo = normalize(bb);
        end
        
        function [gpy, loc, rscl] = find_loc_pyramid(cc, ref, img, ref_flag, n)
            %%% parameters %%%
            if nargin < 4 || isempty(ref_flag)
                ref_flag = 1;
            end
            rt = linspace(0.3, 1, 40);
            istep = 0.5; %%% mm %%%
            res = 0.01; %% mm %%
            istepf = istep / res;
            [nh, nw, nn] = size(ref);
            srg = 5;
            if nargin < 5 || isempty(n)
                rg = 2: istepf: nn - 1;
            else
                rg = max(2, n - srg): min(nn - 1, n + srg);
            end
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            rscl = opt.rscl;
            opt.cfflag = 1; %%% 1 for pyramid, 2 for coarse and 3 for fine %%%
            
            %%% slice prep %%%
            cs = zeros(length(rg), length(rt));
            ss = cs;
            mis = cs;
            ks = cs;
            locs = cell(length(rt), 1);
            reft = ref(:, :, 1);
            reft = imresize(reft, opt.rscl);
            img = imresize(img, size(reft));
            se = 5;
            imgt = cc.anidenoise(cc, img, se);
            imgt = imgt(se + 1: end - se, se + 1: end - se);
            imgt = imgaussfilt(TVL1denoise(normalize(imgt), 0.4, 20) .^ 1, max(1, floor(min(size(reft)) / 100)));
            imgt = imresize(imgt, size(reft));
            
            %%% calculate measures %%%
            for i = 1: length(rg)
                for j = 1: length(rt)
                    bf = imresize(imgt, rt(j));
                    opt.bf = bf;
                    [cst, sst, mit, locst, kst] = cc.mea_cal_unit(cc, ref(:, :, rg(i)), opt);
                    cs(i, j) = cst;
                    ss(i, j) = sst;
                    mis(i, j) = mit;
                    ks(i, j) = kst;
                    locs{j}(i, :) = locst;
                end
            end
            
            tmp = ss;
            tmp = max(tmp, [], 1);
            [~, i] = max(tmp);
            gpy = rt(i);
            loc = round(mean(locs{i}, 1));
        end
        
        function idf = find_loc_coarse(cc, ref, img, ref_flag, n, gpy, loc, rsclp)
            %%% parameters %%%
            if nargin < 4 || isempty(ref_flag)
                ref_flag = 1;
            end
            istep = 0.5; %%% mm %%%
            res = 0.01; %% mm %%
            istepf = istep / res;
            [nh, nw, nn] = size(ref);
            if nargin < 5 || isempty(n)
                rg = 2: istepf: nn - 1;
            else
                rg = max(2, n - 10): min(nn - 1, n + 10);
            end
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            lrt = opt.rscl / rsclp;
            opt.cfflag = 2; %%% 1 for pyramid, 2 for coarse and 3 for fine %%%
            
            %%% slice prep %%%
            cs = zeros(length(rg), 1);
            ss = cs;
            mis = cs;
            ks = cs;
            reft = ref(:, :, 1);
            img = imresize(img, size(reft));
            se = 5;
            imgt = cc.anidenoise(cc, img, se);
            imgt = imgt(se + 1: end - se, se + 1: end - se);
            imgt = imgaussfilt(normalize(TVL1denoise(imgt, 0.4, 20)) .^ 2, floor(min(size(reft)) / 100));
            imgt = imresize(imgt, size(reft));
            bf = imresize(imgt, gpy);
            opt.bf = bf;
            opt.loc = loc * lrt;
            
            %%% calculate measures %%%
            for i = 1: length(rg)
                [cst, sst, mit, ~, kst] = cc.mea_cal_unit(cc, ref(:, :, rg(i)), opt);
                cs(i) = cst;
                ss(i) = sst;
                mis(i) = mit;
                ks(i) = kst;
            end
            
            tmp = cs + ss;
            [~, i] = max(tmp);
            idf = rg(i);
        end
        
        function [theta, gamma, ap, rloc, imgo] = find_loc_fine(cc, ref, img, ref_flag, n, gpy, loc, rtl, rsclp)
            %%% generate interpolation grid %%%
            [nh, nw, nf] = size(ref);
            mx = 2;
            astep = 0.5;
            apmx = 2;
            mxtheta = atand(min(nf - n - apmx, n - 1 - apmx) / (nw / 2));
            mxgamma = atand(min(nf - n - apmx, n - 1 - apmx) / (nh / 2));
            idtheta = max(-mx, -mxtheta): astep: min(mx, mxtheta); %%% left right %%%
            idgamma = max(-mx, -mxgamma): astep: min(mx, mxgamma); %%% up down %%%
            idap = -apmx: apmx;
            nth = length(idtheta);
            nga = length(idgamma);
            nap = length(idap);
            reft = permute(ref, [3, 1, 2]);
            reft = griddedInterpolant(double(reft));
            opt.ref_flag = ref_flag;
            opt.rscl = 1;
            opt.cfflag = 3;
            lrt = opt.rscl / rsclp;
            opt.loc = lrt * loc;
            
            %%% slice prep %%%
            img = imresize(img, [nh, nw]);
            imgt = imgaussfilt(normalize(TVL1denoise(img, 0.5, 40)) .^ 1, max(1, floor(min(size(reft.Values)) / 100)));
            bf = imresize(imgt, gpy);
            opt.bf = bf;
            cs = zeros(nth * nga * nap, 1);
            ss = cs;
            mis = cs;
            ks = cs;

            coord = cc.prep_coord(nh, nw);
            imgs = zeros(nh, nw, length(idap) * length(idtheta) * length(idgamma));
            for i = 1: nth
                for j = 1: nga
                    for k = 1: nap
                        theta = idtheta(i);
                        gamma = idgamma(j);
                        ap = idap(k) + n;
                        angles.theta = theta;
                        angles.gamma = gamma;
                        angles.ap = ap;
                        
                        coordf = cc.slice_coords(coord, angles, [nh, nw]);
                        imgt = reft(rtl * (coordf(:, :, 3) - ap) + ap, coordf(:, :, 2), coordf(:, :, 1));
                        count = (i - 1) * nga * nap + (j - 1) * nap + k;
                        imgs(:, :, count) = imgt;
                    end
                end
            end
            
            tic
            ns = size(imgs, 3);
            rlocs = zeros(ns, 4);
            if isempty(gcp('nocreate'))
                parpool(feature('numCores'));
            end
            parfor i = 1: ns
                [cst, sst, mit, rloc, kst] = cc.mea_cal_unit(cc, imgs(:, :, i), opt);
                cs(i) = cst;
                ss(i) = sst;
                mis(i) = mit;
                rlocs(i, :) = rloc;
                ks(i) = kst;
                if mod(i, floor(ns / 10)) == 0
                    disp(['done ', num2str(i), '/', num2str(ns)])
                end
            end
            toc
            
            tmp = ks .* mis;
            [~, ii] = max(tmp(:));
            iu = floor(ii / (nga * nap));
            ju = floor((ii - iu * nga * nap) / nap);
            ku = mod(ii - iu * nga * nap - ju * nap, nap);
            iu = iu + 1;
            ju = ju + 1;
            if ku == 0
                ku = nap;
            end
            
            theta = idtheta(iu);
            gamma = idgamma(ju);
            ap = idap(ku);
            rloc = rlocs(ii, :);
            imgo = imgs(:, :, ii);
        end
        
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
        
        function [imgg, D, rloct, imref, coordft] = regist(cc, tv, sl, angles, rloc, imro)
            %%% estimate tv slices to be used %%%
            rt = 0.00;
            [nf, nh, nw] = size(tv);
            coord = cc.prep_coord(nh, nw);
            coordft = cc.slice_coords(coord, angles, [nh, nw]);
            rg = [floor(min(min(coordft(:, :, 3)))), ceil(max(max(coordft(:, :, 3))))];
            reft = griddedInterpolant(double(tv(rg(1): rg(2), :, :)));
            coordf = coordft;
            coordf(:, :, 3) = coordft(:, :, 3) - rg(1);
            rloct = round(reshape(rloc, 2, 2)' * [(1 + rt), -rt; -rt, (1 + rt)])';
            rloct = round(nh / size(imro, 1)) * rloct(:);
            rloct(1) = max(1, rloct(1));
            rloct(3) = max(1, rloct(3));
            rloct(2) = min(nh, rloct(2));
            rloct(4) = min(nw, rloct(4));
            
            %%% imcur %%%
            sls = repmat(zeros(rloct(2) - rloct(1) + 1, rloct(4) - rloct(3) + 1), 1, 1, size(sl, 3));
            for i = 1: size(sl, 3)
                sls(:, :, i) = imresize(sl(:, :, i), size(sls(:, :, 1)));
            end
            se = 5;
            imcur = cc.anidenoise(cc, sls(:, :, 3), se);
            imcur = imcur(se + 1: end - se, se + 1: end - se);
            imcurt = normalize(imgaussfilt(TVL1denoise(normalize(imcur), 0.4, 20) .^ 1, max(1, floor(min(size(imcur)) / 100))));
            
            %%% imref %%%
            imref = normalize(reft(coordf(:, :, 3), coordf(:, :, 2), coordf(:, :, 1)));
            imref = cc.anidenoise(cc, imref, se);
            imref = imref(se + 1: end - se, se + 1: end - se);
            mask = normalize(imref) > 0.05;
            mask = imfill(mask, 'holes');
            imreft = imref .* mask;
            %%%% test exponent %%%%
            imreft = normalize(imreft(rloct(1): rloct(2), rloct(3): rloct(4)));
            %%%% apply exponent %%%%
            imreft = normalize(imgaussfilt(TVL1denoise(normalize(imreft .^ 1), 0.4, 40), max(1, floor(min(size(imref)) / 400))));
            
% %             [optimizer, metric] = imregconfig('multimodal');
% %             optimizer.InitialRadius = 0.009;
% %             optimizer.Epsilon = 1.5e-4;
% %             optimizer.GrowthFactor = 1.02;
% %             optimizer.MaximumIterations =100;
% % %             imgt = imregister(imcurt, imreft, 'affine', optimizer, metric);
% %             tf = imregtform(imcurt, imreft, 'affine', optimizer, metric);
% %             
% %             for i = 1: size(sls, 3)
% %                 sls(:, :, i) = imwarp(sls(:, :, i), tf, 'outputview', imref2d(size(imreft)));
% %             end
% % %             imcurtt = imwarp(imcurt, tf, 'outputview', imref2d(size(imreft)));
% %             imcurt = imwarp(imcur, tf, 'outputview', imref2d(size(imreft)));
            imcurt = cc.opt_exp(cc, imcurt, imreft);
            imreft = cc.opt_exp(cc, imreft, imcurt);
            
            %%% sharpen images %%%
            imreft = imsharpen(imreft, 'radius', 2, 'amount', 2);
            imcurt = imsharpen(imcurt, 'radius', 2, 'amount', 2);
            
            imreft = max(0, imreft);
            imcurt = max(0, imcurt);
                        
            disp('doing registration')
            try
                [D, img]=imregdemons(gpuArray(imcurt), gpuArray(imreft), 400, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);
                D = gather(D);
            catch
                [D, img]=imregdemons(imcurt, imreft, 400, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);
            end
            disp('done registration')
            imgg = sls;
            for i = 1: size(sl, 3)
                imgg(:, :, i) = cc.iminterpolate(sls(:, :, i), D(:, :, 2), D(:, :, 1));
            end
        end
        
        function pl = cell_detect(cc, a, opt)
            if nargin < 3 || isempty(opt)
                vis = false;
            else
                vis = opt.vis;
            end
            
            %%% cell enhancing %%%
            a = normalize(double(a));
            a0 = normalize(imgaussfilt(a, 0.5) .^ 2);
            
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
            pl = {[y, x]};
            
            %%% visualize %%%
            if vis
                figure(gcf)
                imagesc(a)
                hold on
                for i = 1: length(y)
                    plot(x(i), y(i), '.r')
                end
                hold off
            end
        end
        
        function pln = cell_transform(D, pl, rloc, opt)
            [nh, nw, nd] = size(D);
            rh = nh / opt.ssl(1);
            wh = nw / opt.ssl(2);
            pltmp = round(pl{1} .* [rh, wh]);
            [x, y] = ndgrid(0: nh - 1, 0: nw - 1); % coordinate image
            x_prime = x + D(:, :, 2); % updated x values (1st dim, rows)
            y_prime = y + D(:, :, 1); % updated y values (2nd dim, cols)
                        
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
        
        function pl3 = cell3d(pl, coordf)
            tmp1 = griddedInterpolant(coordf(:, :, 1));
            tmp2 = griddedInterpolant(coordf(:, :, 2));
            tmp3 = griddedInterpolant(coordf(:, :, 3));
            pl3 = zeros(size(pl, 1), 3);
            for i = 1: size(pl, 1)
                pl3(i, :) = [tmp1(pl(i, :)), tmp2(pl(i, :)), tmp3(pl(i, :))];
            end
        end
        
        function save(cc)
            fn = [cc.pname, filesep, 'cell_count_results.mat'];
            save(fn, 'cc')
        end
        
        function [mask, pll] = lens_loc(cc, img, imr, rlocn)
            t = imr(rlocn(1): rlocn(2), rlocn(3): rlocn(4));
            t1 = cc.opt_exp(cc, t, img(:, :, 3));
            t2 = cc.opt_exp(cc, img(:, :, 3), t1);
            m1 = imfill(t1 > 0.01, 'holes');
            m2 = imfill(t2 > 0.01, 'holes');
            mask1 = xor(m1, m2);
            mask2 = imopen(mask1, strel('disk', 5));
            [l, n] = bwlabeln(mask2);
            s = zeros(n, 1);
            for i = 1: n
                tmp = l == i;
                s(i) = sum(tmp(:));
            end
            [~, id] = max(s);
            mask = l == id;
            [h, w] = find(mask);
            pll = [h + rlocn(1), w + rlocn(3)];
        end

        %% aux functions: utility %%
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
            tp = zeros([size(imref), size(img, 3)]);
            tp(rloc(1): rloc(2), rloc(3): rloc(4), :) = img;
            tp =  imfuse(imref * 4, tp * 4, 'blend', 'Scaling', 'none');
        end
        
        function cell_overlay(img, pl)
            figure(gcf)
            clf
            imagesc(img)
            hold on
            for i = 1: size(pl, 1)
                plot(pl(i, 2), pl(i, 1), '.r')
            end
            hold off
        end
    end
end




