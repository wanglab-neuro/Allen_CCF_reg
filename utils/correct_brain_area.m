[fname, pname] = uigetfile('*.mat', 'Select imaging file', 'MultiSelect', 'off');
load(fname)

if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    ppname = mfilename('fullpath');
    load(fullfile(ppname, 'AllenAtlas.mat'))
end
bregma = allenCCFbregma; % bregma position in reference data space
atlas_resolution = 0.010; % mm


n = length(roi_table);
for i = 1: n
    nn = size(roi_table_group{i}, 1);
    for j = 1: nn
        figure(1)
        clf
        name_cur = roi_table_group{i}.name{j};
        islice = strcmp(roi_table{i}.name, name_cur);
        slices = table2cell(roi_table{i});
        slices = cell2mat(slices(:, 3));
        islice = unique(-slices(islice) / atlas_resolution + bregma(1));
        imagesc(squeeze(mean(tv(islice, :, :), 1)))
        hold on
        for ii = 1: size(roi_table{i}, 1)
            if strcmp(roi_table{i}.name{ii}, name_cur)
                plot(roi_table{i}.ML_location(ii) / atlas_resolution + bregma(3), roi_table{i}.DV_location(ii) / atlas_resolution + bregma(2), '.r', 'markersize', 8)
            end
        end
        hold off
        title(['Cells labeled for ', roi_table_group{i}.name{j}])
        waitforbuttonpress;
    end
end


