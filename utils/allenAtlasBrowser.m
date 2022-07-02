function f = allenAtlasBrowser(f, templateVolume, annotationVolume, structureTree, slice_figure, fname, curr_slice, curr_angle, cell_detect_mode)
    % create figure and adjust to user's screen size
    % f = figure('Name','Atlas Viewer'); 
    figure(f);
    try
        screen_size = get(0, 'ScreenSize'); 
        screen_size = screen_size(1,3:4)./[2560 1440];
    catch 
        screen_size = [1900 1080]./[2560 1440];
    end 
    set(f,'Position', [1050*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
    movegui(f,'onscreen')

    % initialize user data variables held by the figure
    Ft = griddedInterpolant(single(templateVolume));
    Fa = griddedInterpolant(single(annotationVolume));
    ud.Ft = Ft;
    ud.Fa = Fa;
    
    ud.bregma = allenCCFbregma; 
    if isempty(curr_slice)
        ud.currentSlice = ud.bregma(1);
    else
        ud.currentSlice = curr_slice;
    end
    
    if isempty(curr_angle)
        ud.currentAngle = zeros(2,1);
    else
        ud.currentAngle = curr_angle;
    end
    
    if nargin < 9 || isempty(cell_detect_mode)
        ud.cell_detect_mode = 1; % 1 for auto, others for manual %
    else
        ud.cell_detect_mode = cell_detect_mode;
    end
    
    ud.scrollMode = 0;

%     ud.transform_type = 'projective'; %can change to 'affine' or 'pwl'
    ud.transform_type = 'pwl'; %can change to 'affine' or 'pwl'

    ud.oldContour = [];
    ud.showContour = false;
    ud.showOverlay = false; ud.overlayAx = [];
    ud.pointList = cell(1,3); ud.pointList{1} = zeros(0,3); 
    ud.pointHands = cell(1,3);
    ud.probe_view_mode = false;
    ud.currentProbe = 0; ud.ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .5 0 .6; 1 .6 0]; 
    ud.ProbeColor =  {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','red','purple','orange'};
    ud.getPoint_for_transform =false; ud.pointList_for_transform = zeros(0,2); ud.pointHands_for_transform = [];
    ud.current_pointList_for_transform = zeros(0,2); ud.curr_slice_num = 1;
    ud.clicked = false;
    ud.showAtlas = false;
    ud.viewColorAtlas = false;
    ud.histology_overlay = 0; 
    ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);
    ud.transform = [];
    ud.transformed_slice_figure = [];
    ud.slice_shift = 0;
    ud.loaded_slice = 0;
    ud.slice_at_shift_start = 1;
    ud.text = [];

    reference_image = squeeze(templateVolume(ud.currentSlice,:,:));
    ud.im = plotTVslice(reference_image);
    ud.ref_size = size(reference_image);
    ud.ref = squeeze(templateVolume(ud.currentSlice,:,:));
    ud.curr_im = squeeze(templateVolume(ud.currentSlice,:,:));
    ud.curr_slice_trans = squeeze(templateVolume(ud.currentSlice,:,:));
    ud.im_annotation = squeeze(annotationVolume(ud.currentSlice,:,:));
    ud.atlas_boundaries = zeros(ud.ref_size,'uint16');
    ud.offset_map = zeros(ud.ref_size);
    ud.loaded = 0;

    % create functions needed to interact with the figure
    ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
        'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');
    
    ud.angleText = annotation('textbox', [.7 0.95 0.4 0.05], ...
        'EdgeColor', 'none', 'Color', 'k');
    allData.tv = templateVolume;
    allData.av = annotationVolume;
    allData.st = structureTree;
    hold(ud.atlasAx, 'on');
    set(ud.atlasAx, 'HitTest', 'off');
%     set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, slice_figure, k, allData, save_location, save_suffix));
    set(f, 'UserData', ud);
    set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure));
    set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, slice_figure, k, allData, fname));
    set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData))
    set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData)); % Set the motion detector.
end

function atlasClickCallback(im, keydata, slice_figure)
    f = get(get(im, 'Parent'), 'Parent');
    ud = get(f, 'UserData');
    ud_slice = get(slice_figure, 'UserData');


    if ud.getPoint_for_transform
        clickX = round(keydata.IntersectionPoint(1));
        clickY = round(keydata.IntersectionPoint(2));
        if ud.showOverlay; clickY = size(ud.ref,1) - clickY; end
        
        if ud.curr_slice_num ~= ud.slice_at_shift_start+ud.slice_shift
            if ~ud.loaded
                ud.current_pointList_for_transform = zeros(0,2);
                disp('transforming new slice');
            end
        end
        ud.pointList_for_transform(end+1, :) = [clickX, clickY];
        ud.current_pointList_for_transform(end+1, :) = [clickX, clickY];
        set(ud.pointHands_for_transform(:), 'color', [.7 .3 .3]);
        ud.pointHands_for_transform(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro', 'color', [0 .9 0],'LineWidth',2,'markers',4);
        
        ud.slice_at_shift_start = ud_slice.slice_num;
        ud.slice_shift = 0;
        ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift;
        ud.loaded = 0;
        ud.clicked = true;
    end
    set(f, 'UserData', ud);
end

function hotkeyFcn(f, fs, keydata, allData, fname)
    % retrieve user data from figure
    ud = get(f, 'UserData');
    ud_slice = get(fs, 'UserData');
    
    key_letter = lower(keydata.Key);
    switch key_letter
        case 't'
            ud.getPoint_for_transform = ~ud.getPoint_for_transform;
            ud.loaded = false;
            
            if ud.getPoint_for_transform
                disp('transform point mode on');
                
                ud.currentProbe = 0;
                
                % launch transform point mode
                if ~size(ud.current_pointList_for_transform,1)% ) (ud.curr_slice_num ~= (ud.slice_at_shift_start+ud.slice_shift) ||  && ~ud.loaded
                    ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift; %ud_slice.slice_num;
                    ud.current_pointList_for_transform = zeros(0,2);
                    set(ud.pointHands_for_transform(:), 'Visible', 'off');
                    num_hist_points = size(ud_slice.pointList,1);
                    template_point = 1; template_points_shown = 0;
                    updateBoundaries(f,ud, allData); ud = get(f, 'UserData');
                end
            else
                disp('transform point mode OFF');
            end
        case 'd'
            if ud.getPoint_for_transform
                ud.current_pointList_for_transform = zeros(0,2); set(ud.pointHands_for_transform(:), 'Visible', 'off');
                ud.pointHands_for_transform = []; 
                ud_slice.pointList = []; 
                set(fs, 'UserData', ud_slice);
                disp('current transform erased');
            end
        case 'h'
%             trans = slice_trans(f, fs);
%             ud.trans = trans;
% %             trans = double(roipoly(uint8(ud.ref)));
            fcur = imagesc(uint8(ud.ref));
            temp = drawpolygon;
            if size(temp.Position, 1) > 1
                trans = double(roipoly(uint8(ud.ref), temp.Position(:, 1), temp.Position(:, 2)));
            else
                trans = true(size(ud.ref));
            end
            if isempty(trans)
                trans = true(size(ud.ref));
            end
            ud.trans = trans;
            [D, mask1, mask2] = overlay_slice(f, fs, trans);
            ud.D = D;
            ud.mask1 = mask1;
            ud.mask2 = mask2;
        case 'p'
            sliceBrowser(fs, ud_slice.processed_images_folder, f, size(allData.tv), fname, ud.cell_detect_mode)
        case 's'
            pl = overlay_cell(f, fs, ud);
            if ~isempty(pl)
                ud.pl = pl;
                save(fullfile(ud_slice.processed_images_folder, ['cell_points_', ud_slice.processed_image_names{1}(1: end - 4), '.mat']), 'pl')
                disp('cell points saved');
            else
                disp('cell points not saved');
            end
        case 'a'
            objectPoints = ud.pl;
            objects = 1: length(objectPoints);
            av = allData.av;
            st = allData.st;
            
            % initialize cell array containing info on each clicked point
            if length(objects) > 1
                roi_annotation = cell(length(objects),1);
                roi_location = cell(length(objects),1);
            end
            
            % generate needed values
            bregma = allenCCFbregma; % bregma position in reference data space
            atlas_resolution = 0.010; % mm
            
            fwireframe = plotBrainGrid([], [], [], true); hold on;
            fwireframe.InvertHardcopy = 'off';
            
            roi_table = cell(1, length(objects));
            roi_table_group = cell(1, length(objects));
            for object_num = objects
                
                selected_object = objects(object_num);
                
                % get the object points for the currently analyzed object
                curr_objectPoints = round(objectPoints{selected_object}(:, [3 2 1]));
                
                % plot points on the wire frame brain
                figure(fwireframe); hold on
                hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.', 'color', ud_slice.CellColors(object_num, :), 'markers', 8);
                
                % use the point's position in the atlas to get the AP, DV, and ML coordinates
                ap = -(curr_objectPoints(:,1)-bregma(1))*atlas_resolution;
                dv = (curr_objectPoints(:,2)-bregma(2))*atlas_resolution;
                ml = (curr_objectPoints(:,3)-bregma(3))*atlas_resolution;
                
                roi_location_curr = [ap dv ml];
                
                % initialize array of region annotations
                roi_annotation_curr = cell(size(curr_objectPoints,1),3);
                
                % loop through every point to get ROI locations and region annotations
                for point = 1:size(curr_objectPoints,1)
                    
                    % find the annotation, name, and acronym of the current ROI pixel
                    ann = av(curr_objectPoints(point,1),curr_objectPoints(point,2),curr_objectPoints(point,3));
                    name = st.safe_name{ann};
                    acr = st.acronym{ann};
                    
                    roi_annotation_curr{point,1} = ann;
                    roi_annotation_curr{point,2} = name;
                    roi_annotation_curr{point,3} = acr;
                    
                end
                
                % save results in cell array
                if length(objects) > 1
                    roi_annotation{object_num} = roi_annotation_curr;
                    roi_location{object_num} = roi_location_curr;
                else
                    roi_annotation = roi_annotation_curr;
                    roi_location = roi_location_curr;
                end
                
                % display results in a table
                disp(['Clicked points for object ' num2str(selected_object)])
                roi_table{object_num} = table(roi_annotation_curr(:,2),roi_annotation_curr(:,3), ...
                    roi_location_curr(:,1),roi_location_curr(:,2),roi_location_curr(:,3), roi_annotation_curr(:,1), ...
                    'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});
                disp(roi_table{object_num})
                
                [rgall, ida, idb] = unique(roi_annotation_curr(:, 2));
                ncell = zeros(length(rgall), 1);
                for j = 1: length(rgall)
                    ncell(j) = sum(idb == j);
                end
                roi_table_group{object_num} = table(rgall, ncell, 'VariableNames', {'name', 'cell_count'});
                disp(roi_table_group{object_num})
            end
            save(fullfile(ud_slice.processed_images_folder, ['table_summary_', ud_slice.processed_image_names{1}(1: end - 4), '.mat']), 'roi_table', 'roi_table_group')
            disp('cell points analyzed');
        case 'uparrow' % scroll angles along D/V axis
            ud.scrollMode = 1;
            disp('switch scroll mode -- tilt D/V')
        case 'rightarrow' % scroll angles along M/L axis
            ud.scrollMode = 2;
            disp('switch scroll mode -- tilt M/L')
        case 'downarrow' % scroll along A/P axis
            ud.scrollMode = 0;
            disp('switch scroll mode -- scroll along A/P axis')
    end

    set(f, 'UserData', ud);
end
    
function pl = overlay_cell(f, fs, udf)
    udfs = get(fs, 'UserData');
    if ~isfield(udfs,'D')
        disp('No slice registration found. Press "h" to select the outline you want the slice to be registered to in the atlas viewer')
        pl=[];
        return
    else
        D = udf.D;
    end
    thetay = udf.currentAngle(1);
    thetaz = udf.currentAngle(2);
    h1 = udf.mask1;
    h2 = udf.mask2;
%     [Ry, Rz] = rotation_matrix(udf.currentAngle(1), udf.currentAngle(2));
    pp = udfs.pointList;    
    pl = cell(1, length(pp));
    for ii = 1: length(pp)
        pr = zeros(size(pp{ii}, 1), 3);
        for j = 1: size(pp{ii}, 1)
            mtx = zeros(size(udfs.current_slice_image(:, :, 1)));
            idd = sub2ind(size(udfs.current_slice_image), round(pp{ii}(j, 1)), round(pp{ii}(j, 2)));
            mtx(idd) = 1;
            mtx = mtx(h1(1): h1(2), h1(3): h1(4));
            mtx = imresize(mtx, [h2(2) - h2(1) + 1, h2(4) - h2(3) + 1]);
            
            mtxt = zeros(size(D(:, :, 1)));
            mtxt(h2(1): h2(2), h2(3): h2(4)) = mtx;
            
            img = iminterpolate(mtxt, D(:, :, 2), D(:, :, 1));
            tmp = imgaussfilt(img, 1);
            tmp = double(normalize(tmp) > 0.3);
            [l, ~] = bwlabeln(tmp);
            ptmp = zeros(1, 2);
            tt = regionprops(l == 1, 'Centroid');
            ptmp(1, :) = fliplr(tt.Centroid);
            ptmp = ptmp + 0.1 * rand(size(ptmp));
        
            prrn = zeros(size(ptmp, 1), 3);
            prrn(:, 1) = udf.currentSlice;
            prrn(:, 2) = prrn(:, 2) + ptmp(:, 2);
            prrn(:, 3) = prrn(:, 3) + ptmp(:, 1);
            pr(j, :) = prrn(:, [2, 3, 1]);
        end
        
        figure(f)
        hold on
        for i = 1: size(pr, 1)
            plot(gca, pr(i, 1), pr(i, 2), '.r', 'color', udfs.CellColors(ii, :), 'markersize', 8)
        end
        hold off
        
        %%% transform the points to 3D coordinates %%%
        dp = pr(:, [3, 1, 2]);
        n1 = udf.ref_size(1);
        n2 = udf.ref_size(2);
%         dp(:, 3) = n1 - dp(:, 3);
        b = [dp(1, 1), n2 / 2, n1 / 2];
        dp = dp - b;
        [Ry, Rz] = rotation_matrix(thetay, thetaz);
        dpf = dp * Ry' * Rz';
        dpf = dpf + b;
%         dpf(:, 3) = n1 - dpf(:, 3);
        pr = dpf;
        
        pr = pr(:, [2, 3, 1]);
        pl{ii} = pr;
        
    end
end

function [D, mask1, mask2] = overlay_slice(f, fs, trans)
    udf = get(f, 'UserData');
    udfs = get(fs, 'UserData');
    % remove color atlas
    udf.viewColorAtlas = false;
    set(udf.im, 'CData', udf.ref)
    colormap(udf.atlasAx, 'gray'); caxis(udf.atlasAx, [0 400]);
    % remove overlayt
    udf.showOverlay = 0;
    delete(udf.overlayAx); udf.overlayAx = [];
    
    %%% fuse the images %%%
    imgref = normalize(double(udf.ref) .* trans);
    imgslice = double(udfs.current_slice_image);
    imgcur = normalize(imgslice(:, :, 3));    
    mask = imgcur > 0.1;
    hmin = find(sum(mask, 2), 1);
    hmax = find(sum(mask, 2), 1, 'last');
    wmin = find(sum(mask, 1), 1);
    wmax = find(sum(mask, 1), 1, 'last');
    imgcur = imgcur(hmin: hmax, wmin: wmax);
    imgslice = imgslice(hmin: hmax, wmin: wmax, :);  
    mask1 = [hmin, hmax, wmin, wmax];
    
    bw = imgref > 0.1;
    hmin = find(sum(bw, 2), 1);
    hmax = find(sum(bw, 2), 1, 'last');
    wmin = find(sum(bw, 1), 1);
    wmax = find(sum(bw, 1), 1, 'last');
    mask2 = [hmin, hmax, wmin, wmax];
    
    imgs = zeros([size(imgref), 3]);
    for i = 1: size(imgslice, 3)
        imgs(hmin: hmax, wmin: wmax, i) = normalize(imresize(imgslice(:, :, i), [hmax - hmin + 1, wmax - wmin + 1]));
    end
    
    mask = imgs(:, :, 3) > 0.1;
    for i = 1: size(imgs, 3)
        imgs(:, :, i) = imgs(:, :, i) .* mask;  
        imgs(imgs < 0.1) = 0;
    end
    imgcur = imgs(:, :, 3);
    
%     [~,imgcur,~]=klt_ref_track(imgcur,imgref);
%     imref = normalize(imgaussfilt(imgref, 19));
%     imcur = normalize(imgaussfilt(imgcur, 19));
%     imref = TVL1denoise(imref, 0.1);
%     imcur = TVL1denoise(imcur, 0.1);
    disp('prepare for registration')
% % %     imref = normalize(feature2_comp(imgref, [], [], 4));
% % %     imcur = normalize(feature2_comp(imgcur, [], [], 4));
% % % %     imref = imgaussfilt(imref, 13);
% % % %     imcur = imgaussfilt(imcur, 13);
% % %     imref = anidenoise(imref, 0, 0, 100, 0.2, 2);
% % %     imcur = anidenoise(imcur, 0, 0, 100, 0.2, 2);
    
%     denomref = imgaussfilt(imgref, 99);
%     denomcur = imgaussfilt(imgcur, 99);
%     imref = anidenoise(imgref, 0, 0, 10, 0.2, 2);
%     imcur = anidenoise(imgcur, 0, 0, 10, 0.2, 2);
    imref = normalize(imgref);
    imcur = normalize(imgcur);
% %     [l, n] = bwlabeln(imcur > 0.1);
% %     s = zeros(1, n);
% %     for k = 1: n
% %         t = l == k;
% %         s(k) = sum(t(:));
% %     end
% %     [~, id] = max(s);
% %     l = l == id;
% %     imcur = imcur .* l;
%     imref = imgaussfilt(imgref, 9);
%     imcur = imgaussfilt(imgcur, 9);
%     denomref = min(imref(:)) + 0.5 * (max(imref(:)) - min(imref(:))) * normalize(imgaussfilt(TVL1denoise(imref, 0.2, 10), 3));
%     denomcur = min(imcur(:)) + 0.5 * (max(imcur(:)) - min(imcur(:))) * normalize(imgaussfilt(TVL1denoise(imcur, 0.2, 10), 3));
    imref = imgaussfilt(TVL1denoise(imref, 0.2, 10), 3);
    imcur = imgaussfilt(TVL1denoise(imcur, 0.2, 10), 3);
    denomref = imgaussfilt(imref, 13);
    denomcur = imgaussfilt(imcur, 13);
    imref = normalize(feature2_comp(imref, [], [], denomref / 1.2));
    imcur = normalize(feature2_comp(imcur, [], [], denomcur / 1.2));
%         imref = imgaussfilt(imref, 5);
%         imcur = imgaussfilt(imcur, 5);
    imref = anidenoise(imref, 0, 0, 50, 0.2, 2);
    imcur = TVL1denoise(imcur, 0.1, 20);
    imcur = anidenoise(imcur, 0, 0, 50, 0.2, 2);
    disp('doing registration')
    try
        [D, img]=imregdemons(gpuArray(imcur), gpuArray(imref), 1000, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);
        D = gather(D);
    catch
        [D, img]=imregdemons(imcur, imref, 1000, 'AccumulatedFieldSmoothing', 3, 'displaywaitbar', false);        
    end
    disp('done registration')
%     img = imwarp(imgcur, trans, 'outputview', imref2d(size(imgref)));

    imgg = imgs;
    for i = 1: size(imgg, 3)
        imgg(:, :, i) = iminterpolate(normalize(imgg(:, :, i)), D(:, :, 2), D(:, :, 1)); 
    end
    imgref = normalize(double(udf.ref));
    imgref = repmat(imgref, 1, 1, 3);
    imgref(:, :, 3) = 0;
    tp = imgg;
    tp(:, :, 1) = 0;
    image_blend =  imfuse(imgref * 4, tp * 4, 'blend', 'Scaling', 'none');
    set(udf.im, 'CData', image_blend);
    imshow(image_blend)
    udf.D = D;
    set(f, 'UserData', udf);
end

function trans = slice_trans(f, fs)
    udf = get(f, 'UserData');
    udfs = get(fs, 'UserData');
    ttype = udf.transform_type;
    ptar = udfs.pointList{1};
    pref = udf.current_pointList_for_transform;
    
    %%% transform %%%
    if ttype == 'lwm'
        trans = fitgeotrans(ptar, pref, ttype, 6);
    else
        trans = fitgeotrans(ptar, pref, ttype);
    end
end

% -----------------------------------------
% Update slice (from scrolling or loading)
% -----------------------------------------
function updateSlice(f, evt, allData)

    ud = get(f, 'UserData');

    % scroll through slices
    if ud.scrollMode==0
        ud.currentSlice = ud.currentSlice+evt.VerticalScrollCount*2;

        if ud.currentSlice>size(allData.tv,1); ud.currentSlice = 1; end %wrap around
        if ud.currentSlice<1; ud.currentSlice = size(allData.tv,1); end %wrap around
    elseif ud.scrollMode == 1
        ud.currentAngle(1) = ud.currentAngle(1)+evt.VerticalScrollCount*2*0.1;
    elseif ud.scrollMode == 2
        ud.currentAngle(2) = ud.currentAngle(2)+evt.VerticalScrollCount*2*0.1;
    end  

    % update coordinates at the top
    pixel = getPixel(ud.atlasAx);
    updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText, ud.angleText, ud.currentSlice, ud.currentAngle(1), ud.currentAngle(2), ud.ref_size);
    
    % ----------------------------------------
    % if no angle, just change reference slice
    % ----------------------------------------
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        
        reference_slice = squeeze(allData.tv(ud.currentSlice,:,:));
        ud.im_annotation = squeeze(allData.av(ud.currentSlice,:,:));
        
        if ud.viewColorAtlas
            set(ud.im, 'CData', ud.im_annotation);
        else
            set(ud.im, 'CData', reference_slice);
        end
        
        
        % update title/overlay with brain region
        [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);
        updateTitle(ud.atlasAx, name, acr)
        if ud.showOverlay
            updateOverlay(f, allData, ann);
        end
        ud.ref = reference_slice;
        set(ud.pointHands_for_transform(:), 'Visible', 'off');
        ud.offset_map = zeros(ud.ref_size);
        set(f, 'UserData', ud);
        
    else
        imszt = size(allData.av);
        imsz = imszt(3: -1: 2);
        imszh = round(imsz / 2);
        [idsy, idsz] = meshgrid((1: imsz(1)) - imszh(1), (1: imsz(2)) - imszh(2));
        dp = [zeros(prod(imsz), 1), idsy(:), idsz(:)];
        thetay = ud.currentAngle(1);
        thetaz = ud.currentAngle(2);
        [Ry, Rz] = rotation_matrix(thetay, thetaz);
        dpf = Rz * Ry * dp';
        dpf = reshape(dpf', [size(idsy), 3]);
        dpf(:, :, 1) = dpf(:, :, 1) + ud.currentSlice;
        dpf(:, :, 2) = dpf(:, :, 2) + imszh(1);
        dpf(:, :, 3) = dpf(:, :, 3) + imszh(2);
        Ft = ud.Ft;
        angle_slice = Ft(dpf(:,:,1), dpf(:,:,3), dpf(:,:,2));
        Fa = ud.Fa;
        ud.im_annotation = Fa(dpf(:,:,1), dpf(:,:,3), dpf(:,:,2));
        
        %%% add mesh 3D true location %%%
        ud.dpf = dpf;

        if ud.viewColorAtlas
            set(ud.im, 'CData', ud.im_annotation);
        elseif ~ud.showAtlas
            set(ud.im, 'CData', angle_slice);
        end
        
        ud.ref = angle_slice;
        set(ud.pointHands_for_transform(:), 'Visible', 'off');
    end
    
    
    % in all cases. update histology overlay
    if ud.histology_overlay == 1 || ud.histology_overlay == 2
        updateHistology(f,ud); ud = get(f, 'UserData');
    else
        if ud.viewColorAtlas
            ud.curr_im = ud.im_annotation;
        else
            ud.curr_im = ud.ref;
        end
    end
    
    % then update boundary overlay
    if ud.showAtlas
        updateBoundaries(f,ud, allData)
    end
        
    set(f, 'UserData', ud);
end

function [Ry, Rz] = rotation_matrix(thetay, thetaz)
    Ry = [cosd(thetay), 0, sind(thetay); 0, 1, 0; -sind(thetay), 0, cosd(thetay)];
    Rz = [cosd(thetaz), -sind(thetaz), 0; sind(thetaz), cosd(thetaz), 0; 0, 0, 1];
end

% ---------------------------------------------------------------
% update the image shown if histology is currently being overlaid
% ---------------------------------------------------------------
function updateHistology(f, ud)
    if ud.histology_overlay == 2
        image_blend =  imfuse(ud.ref*.6, ud.curr_slice_trans(:,:,:),'blend','Scaling','none');
        set(ud.im, 'CData', image_blend);
        ud.curr_im = image_blend;
    elseif ud.histology_overlay == 1
        set(ud.im, 'CData', ud.curr_slice_trans);
        ud.curr_im = ud.curr_slice_trans;
    end
    set(f, 'UserData', ud);
end

    
% -------------------------------------------------    
% update the position of the region boundary image
% -------------------------------------------------
function updateBoundaries(f, ud, allData)
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        curr_annotation = squeeze(allData.av(ud.currentSlice,:,:));
    else
        curr_annotation = ud.im_annotation;
    end
    
    atlas_vert_1 = double(curr_annotation(1:end-2,:));
    atlas_vert_2 = double(curr_annotation(3:end,:));
    atlas_vert_offset = abs( atlas_vert_1 - atlas_vert_2 ) > 0;
    shifted_atlas_vert1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert1(3:end,:) = atlas_vert_offset;
    shifted_atlas_vert2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert2(1:end-2,:) = atlas_vert_offset;

    atlas_horz_1 = double(curr_annotation(:,1:end-2));
    atlas_horz_2 = double(curr_annotation(:,3:end));
    atlas_horz_offset = abs( atlas_horz_1 - atlas_horz_2 )>0;
    shifted_atlas_horz1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz1(:,3:end) = atlas_horz_offset;
    shifted_atlas_horz2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz2(:,1:end-2) = atlas_horz_offset;

    shifted_atlas = shifted_atlas_horz1 + shifted_atlas_horz2 + shifted_atlas_vert1 + shifted_atlas_vert2;

    atlas_boundaries = (shifted_atlas>0); ud.atlas_boundaries = atlas_boundaries;

    if ud.showAtlas
        image_blend =  imfuse(ud.curr_im, atlas_boundaries/3.5*(1+.35*isa(ud.curr_im,'uint16')),'blend','Scaling','none') * 2;
        set(ud.im, 'CData', image_blend); 
    end
    
    set(f, 'UserData', ud);
end
    
% ------------------------
% react to mouse hovering
% ------------------------
function fh_wbmfcn(f, allData)
    % WindowButtonMotionFcn for the figure.

    ud = get(f, 'UserData');
    ax = ud.atlasAx;
    pixel = getPixel(ax);

    %get offset due to angling
    if 0<pixel(1) && pixel(1)<=ud.ref_size(1) && 0<pixel(2) && pixel(2)<=ud.ref_size(2)
        offset = ud.offset_map(pixel(1),pixel(2));
    else; offset = 0;
    end

    % show bregma coords
    updateStereotaxCoords(ud.currentSlice + offset, pixel, ud.bregma, ud.bregmaText, ud.angleText, ud.currentSlice, ud.currentAngle(1), ud.currentAngle(2), ud.ref_size);

    % get annotation for this pixel
    [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice+offset);

    updateTitle(ax, name, acr);

    if ~isempty(name)
        if ud.showContour
            if ~isempty(ud.oldContour)
                delete(ud.oldContour);
            end
            [~,ch] = contour(squeeze(allData.av(ud.currentSlice,:,:)==ann), 1, 'r');
            ud.oldContour = ch;
            set(f, 'UserData', ud);
        end

        if ud.showOverlay
            updateOverlay(f, allData, ann)
        end    
    end
end

% ---------------------------------------------
% update the coordinates shown in the top left
% ---------------------------------------------
function updateStereotaxCoords(currentSlice, pixel, bregma, bregmaText, angleText, slice_num, ap_angle, ml_angle, ref_size)
    atlasRes = 0.010; % mm
    ap = -(currentSlice-bregma(1))*atlasRes;
    dv = (pixel(1)-bregma(2))*atlasRes;
    ml = (pixel(2)-bregma(3))*atlasRes;
    set(bregmaText, 'String', sprintf('%.2f AP, %.2f DV, %.2f ML', ap, dv, ml));
%     set(angleText, 'String', ['Slice ' num2str((bregma(1) - slice_num)/100) ' AP, DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, ML angle ' num2str(round(atand(ml_angle/570),1)) '^{\circ}']);
    set(angleText, 'String', ['Slice ' num2str((bregma(1) - slice_num)/100) ' AP, DV angle ' num2str(ap_angle) '^{\circ}, ML angle ' num2str(ml_angle) '^{\circ}']);
end

% ---------------------------------
% update the current mouse location
% ---------------------------------
function pixel = getPixel(ax)

    currPoint = get(ax,'currentpoint');  % The current point w.r.t the axis.

    Cx = currPoint(1,1); Cy = currPoint(1,2);
    pixel = round([Cy Cx]);
end

% ---------------------------------
% update the overlaid brain region
% ---------------------------------
function updateOverlay(f, allData, ann)
    ud = get(f, 'UserData');
    if isempty(ud.overlayAx) % first time
        if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
            avo = plotAVoverlay(fliplr(squeeze(allData.av(ud.currentSlice,:,:))'), ann, ud.atlasAx);
        else
            avo = plotAVoverlay(fliplr(ud.im_annotation)', ann, ud.atlasAx);
        end
        ud.overlayAx = avo;
        set(ud.overlayAx, 'HitTest', 'off');
        set(f, 'UserData', ud);
    else
        ovIm = get(ud.overlayAx, 'Children');
    %     set(ovIm, 'HitTest', 'off');

        if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
            thisSlice = squeeze(allData.av(ud.currentSlice,:,:));
        else
            thisSlice = ud.im_annotation;
        end

        set(ovIm, 'CData', flipud(thisSlice));    
        plotAVoverlay(fliplr(thisSlice'), ann, ud.atlasAx, ud.overlayAx);
        set(ovIm, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));

    end
end

% ---------------------------------
% find the region being hovered on
% ---------------------------------
function [name, acr, ann] = getPixelAnnotation(allData, pixel, currentSlice)
    if pixel(1)>0&&pixel(1)<size(allData.av,2) && pixel(2)>0&&pixel(2)<=size(allData.av,3)
        ann = allData.av(currentSlice,pixel(1),pixel(2));
        name = allData.st.safe_name(ann);
        acr = allData.st.acronym(ann);
    else
        ann = []; name = []; acr = [];
    end
end

% ---------------------------------
% update the title, showing region
% ---------------------------------
function updateTitle(ax, name, acr)
    if ~isempty(name)
        title(ax, [name{1} ' (' acr{1} ')']);
    else
        title(ax, 'not found');
    end
end

function [f, h] = plotBrainGrid(brainGridData, ax, brain_figure, black_brain)
    % function plotBrainGrid([brainGridData], [ax])
    % 
    % To plot the wire mesh data loaded from brainGridData.npy. 

    if nargin<1 || isempty(brainGridData)
        mf = mfilename('fullpath');
        brainGridData = load(fullfile(fileparts(mf), 'brainGridData.mat'));
        brainGridData = brainGridData.brainGridData;
    end

    bp = double(brainGridData); 
    bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were

    if nargin<2||isempty(ax)    
        if nargin<3||isempty(brain_figure)
            brain_figure = figure('Name','Brain View');
        end

        ax = axes('Parent', brain_figure);  
    end

    if nargin<4||~black_brain
        h = plot3(ax, bp(:,1), bp(:,2), bp(:,3), 'Color', [0 0 0 0.3]);
    else
        h = plot3(ax, bp(:,1), bp(:,2), bp(:,3), 'Color', [.7 .7 .7 0.3]);
        set(get(ax, 'Parent'),'color','k')
    end

    set(ax, 'ZDir', 'reverse')
    axis(ax, 'equal');
    axis(ax, 'vis3d');
    axis(ax, 'off');
    f = get(ax, 'Parent');
end

function im = plotTVslice(thisSlice)
    % function plotTVslice(thisSlice)
    % 
    % useful way to plot a slice you get from sliceByVector, if the slice was
    % made from the allen atlas ccf "template volume"

    im = imagesc(thisSlice);
    % set(gca, 'YDir','normal')
    axis image
    axis off
    colormap gray
    caxis([0 400]);
end