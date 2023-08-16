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
%%% plot av %%%
scl = 0.25;
[~, refa] = cc.prep_ccf(av, scl);
[bpa, mask] = cc.prep_renderer(cc, tv);
[bpa, mx] = cc.merge_cells(cc, bpa, mask);
bpal = cc.create_cyl_label(bpa, model, scl);

vola = cc.plot_3d_brain(bpa, mx);
cc.plot_3d_brain_lens(cc, bpa, bpal, mx, vola)

%% get imaging target regions %%
scl = 0.25;
bpal = cc.target_regions(bpa, model, scl);
[nf, nh, nw] = size(av);
t = imresize3(bpal, [nh, nw, nf]);
avt = permute(av, [2, 3, 1]);
id = avt(t > 0);
id = unique(id);
tnames = st.name(id);

%% visualizing utilities %%
%%% write gif %%%
fn = 'C:\Jinghao\research_temp\allenccf\Lens_track\G014\brain_cell.gif';
cc.brain_animation(vol, fn)

