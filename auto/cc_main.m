%% main %%
%%% parameters %%%
cc = cellcount();
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

%% analyzing lens location %%
m = matfile('C:\Jinghao\research_temp\allenccf\Lens_track\G014\cell_count_results.mat');
cc = m.cc;
lensdata = m.lensdata;
[model, pts, pti] = cc.lens_recon(cc, lensdata);

%% 3d rendering %%
%%% try 3d plot %%%
[bp, mask] = cc.prep_renderer(cc, tv);
[bp, mx] = cc.merge_cells(cc, bp, mask);
bpl = cc.create_cyl_label(bp, model, scl);

vol = cc.plot_3d_brain(bp, mx);
cc.plot_3d_brain_lens(cc, bp, bpl, mx, vol)

%%% plot av %%%
scl = 0.25;
[~, refa] = cc.prep_ccf(av, scl);
[bpa, mask] = cc.prep_renderer(cc, tv);
[bpa, mx] = cc.merge_cells(cc, bpa, mask);
bpal = cc.create_cyl_label(bpa, model, scl);

vola = cc.plot_3d_brain(bpa, mx);
cc.plot_3d_brain_lens(cc, bpa, bpal, mx, vola)

%%% write gif %%%
fn = 'C:\Jinghao\research_temp\allenccf\Lens_track\G014\brain_cell.gif';
cc.brain_animation(vol, fn)

%% visualization for lab meeting %%
%%% prepare proof of principle visualization %%%
figure(1)
clf
figure(2)
clf
nn = length(cc.metadata.slice);
for i = 1: nn
    figure(1)
    subplot(ceil(sqrt(nn)), ceil(sqrt(nn)), i)
    sz = cc.metadata.section_location;
    tp = cc.image_fuse(cc.metadata.reg_ref_image{i}, cc.metadata.reg_image{i}, sz);
    imshow(tp)
    figure(2)
    subplot(ceil(sqrt(nn)), ceil(sqrt(nn)), i)
    cc.cell_overlay(tp, cc.metadata.point_lists_new{i});
end

%%% cell overlay %%%
figure(1)
clf
nn = length(cc.metadata.slice);
for i = 1: nn
    subplot(ceil(sqrt(nn)), ceil(sqrt(nn)), i)
    sz = cc.metadata.section_location;
end

%%% plot old raw images %%%
[file_name, path_name] = cc.data_info();
figure(1)
clf
nn = length(cc.metadata.slice);
for i = 1: nn
    subplot(ceil(sqrt(nn)), ceil(sqrt(nn)), i)
    tp1t = normalize(single(cc.load_slice_unit([cc.metadata.pname, file_name{i}])));
    tp2 = cc.metadata.reg_ref_image{i};
    tp1 = zeros(size(tp2, 1), size(tp2, 2), size(tp1t, 3));
    for j = 1: size(tp1t, 3)
        tp1(:, :, j) = imresize(tp1t(:, :, j), size(tp2));
    end
    tp = cc.image_fuse(tp2, tp1);
    imshow(tp)
end

%%% lens track 2d images %%%
figure(1)
clf
nn = length(cc.metadata.slice);
for i = 1: nn
    subplot(ceil(sqrt(nn)), ceil(sqrt(nn)), i)
    sz = cc.metadata.section_location;
    tp = cc.image_fuse(cc.metadata.reg_ref_image{i}, cc.metadata.reg_image{i}, sz);
    tp = cc.lens_overlay(tp, lensdata.mask{i}, sz);
    imshow(tp)
end

%%% plot lens cylinder model %%%
figure(1)
clf
cc.pc_lens_model(pts, pti, model)


%% pool %%



