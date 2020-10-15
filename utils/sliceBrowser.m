function sliceBrowser(slice_figure, processed_images_folder, f, reference_size, processed_images, cell_detect_mode)

    if nargin < 6 || isempty(cell_detect_mode)
        cell_detect_mode = 1;
    end
    
    % initialize user data variables held by the figure
    ud_slice.processed_image_names = processed_images;
    total_num_files = size(processed_images,1);
    disp(['found ' num2str(total_num_files) ' processed slice images']);
    ud_slice.total_num_files = total_num_files;
    ud_slice.break = 0; 
    ud_slice.slice_num = 1;
    ud_slice.key = 1; 
    ud_slice.pointList = {[]}; 
    ud_slice.pointHands = {[]};
    ud_slice.getPoint = 0;
    ud_slice.ref_size = reference_size(2:3);
    ud_slice.CellColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .5 0 .6; 1 .6 0];
    ud_slice.CellGroupIdx = 1;

    figure(slice_figure)
    ud_slice.processed_images = processed_images;
    ud_slice.processed_images_folder = processed_images_folder;
    ud_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9], 'Ydir', 'reverse');
    hold(ud_slice.sliceAx, 'on');
%     set(ud_sclice.sliceAx, 'HitTest', 'off');
    ud_slice.im = plotTVslice(zeros(ud_slice.ref_size, 'uint8'));

    % create functions needed to interact with the figure
    set(ud_slice.im, 'ButtonDownFcn', @(slice_figure,k)sliceClickCallback(slice_figure, k));
    set(slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(slice_figure, keydata, f));
    set(slice_figure, 'UserData', ud_slice)

    % adjust figure to user's screen size
    try
        screen_size = get(0, 'ScreenSize'); 
        screen_size = screen_size(1,3:4)./[2560 1440];
    catch
        screen_size = [1900 1080]./[2560 1440];
    end
    set(slice_figure,'Position', [150*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
    movegui(slice_figure,'onscreen')

    % set up first slice image
    ud_slice = updateSliceImage(ud_slice, slice_figure);
end

% ------------------------------------------------    
% Clicking function to register transform points  
% ------------------------------------------------
function sliceClickCallback(im, keydata)
    f = get(get(im, 'Parent'), 'Parent');
    ud = get(f, 'UserData');


    if ud.getPoint
        clickX = round(keydata.IntersectionPoint(1));
        clickY = round(keydata.IntersectionPoint(2));

        ud.pointList{ud.CellGroupIdx}(end+1, :) = [clickX, clickY];
        ud.pointHands{ud.CellGroupIdx}(end+1) = plot(ud.sliceAx, clickX, clickY, 'r.', 'color', ud.CellColors(ud.CellGroupIdx, :), 'markersize',10);    

         if clickX < 100 && (ud.ref_size(1) - clickY) < 100 % if click in corner, break
            ud.pointList = []; 
            set(ud.pointHands(:), 'Visible', 'off');     
         end

    end
    set(f, 'UserData', ud);
end

% ------------------------
% react to keyboard press
% ------------------------
function SliceAtlasHotkeyFcn(fig, keydata, f)

    ud = get(fig, 'UserData');
    
    % left arrow -- go to previous slice
    if strcmp(keydata.Key,'leftarrow')
        if ud.slice_num > 1
            ud.slice_num = ud.slice_num - 1;
            ud = updateSliceImage(ud);
        end
        
        % right arrow -- go to next slice
    elseif strcmp(keydata.Key,'rightarrow')
        if ud.slice_num < ud.total_num_files
            ud.slice_num = ud.slice_num + 1;
            ud = updateSliceImage(ud);
        end
        % d -- delete current transform points
    elseif strcmp(keydata.Key,'d')
        disp('current transform points deleted')
        set(ud.pointHands(:), 'Visible', 'off');
        ud.pointList = [];
        % t -- transform point mode
    elseif strcmp(keydata.Key,'t')
        ud.getPoint = ~ud.getPoint;
        if ud.getPoint
            disp('transform point mode on!');
        else
            disp('transform point mode off!')
        end
        % c -- click cell mode
    elseif strcmp(keydata.Key,'c')
        ud.getPoint = ~ud.getPoint;
        if ud.getPoint
            disp('click cell mode on!');
        else
            disp('click cell mode off!')
        end
        % y -- auto cell detect
    elseif strcmp(keydata.Key,'y')
        ud.pointList = cell_detect(double(ud.current_slice_image(:, :, 2)));
        if cell_detect_mode == 2
            ud.pointList{2} = cell_detect(double(ud.current_slice_image(:, :, 1)));
        end
        disp('Done auto cell detect')
        % n -- new cell group in click cell mode
    elseif strcmp(keydata.Key,'n')
        ud.CellGroupIdx = ud.CellGroupIdx + 1;
        ud.pointList{ud.CellGroupIdx} = [];
        ud.pointHands{ud.CellGroupIdx} = [];
        if ud.getPoint
            disp(['Cell Group ', num2str(ud.CellGroupIdx)]);
        end
    else
    % otherwise -- call function to atlas browser       
        figure(f);
        fcn = get(f, 'KeyPressFcn'); 
        fcn(f, keydata);
    end
    set(fig, 'UserData', ud);
end


function ud = updateSliceImage(ud, slice_figure)

    title_ending = '';
    
    processed_image_name = ud.processed_image_names{ud.slice_num};
    current_slice_image = imread(fullfile(ud.processed_images_folder, processed_image_name));
%     if size(current_slice_image,1) > ud.ref_size(1)+2 || size(current_slice_image,2) > ud.ref_size(2)+2
%         disp(['shrinking image to reference size ' num2str(ud.ref_size(1)) ' x ' num2str(ud.ref_size(2)) ' pxl'])
%         current_slice_image = imresize(current_slice_image, ud.ref_size);
%     end  
    ud.current_slice_image = current_slice_image;
    set(slice_figure, 'UserData', ud)
    set(ud.im, 'CData', current_slice_image); 


    file_transformations = fullfile(ud.processed_images_folder, 'transformations\\' ,...
                            [processed_image_name(1:end-4) '_transform_data.mat']);

    set(ud.pointHands{1}(:), 'Visible', 'off'); 

    if exist(file_transformations,'file')
        % load transform data
        transform_data = load(file_transformations);
        transform_data = transform_data.save_transform;
        if ~isempty(transform_data.transform_points{2})
            ud.pointList = transform_data.transform_points{2};
            title_ending = ' (transform points loaded)';
        end
    else
        ud.pointList = [];
    end
    title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files) title_ending])    
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