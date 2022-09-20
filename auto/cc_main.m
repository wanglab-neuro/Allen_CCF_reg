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

