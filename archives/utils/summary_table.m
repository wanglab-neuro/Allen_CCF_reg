function [roi_table, roi_table_group] = summary_table(pls, av, st)
    bregma = allenCCFbregma; % bregma position in reference data space
    atlas_resolution = 0.010; % mm
    roi_table = cell(1, length(pls));
    roi_table_group = cell(1, length(pls));
    for i = 1: length(pls)
        tmp = round(pls{i}(:, [3, 2, 1]));
        roi_annotation_curr = cell(size(tmp, 1),3);
        ap = -(tmp(:,1)-bregma(1))*atlas_resolution;
        dv = (tmp(:,2)-bregma(2))*atlas_resolution;
        ml = (tmp(:,3)-bregma(3))*atlas_resolution;
        roi_location_curr = [ap dv ml];

        for j = 1: size(tmp, 1)
            ann = av(tmp(j, 1), tmp(j, 2), tmp(j, 3));
            name = st.safe_name{ann};
            acr = st.acronym{ann};

            roi_annotation_curr{j,1} = ann;
            roi_annotation_curr{j,2} = name;
            roi_annotation_curr{j,3} = acr;
        end

        roi_table{i} = table(roi_annotation_curr(:,2),roi_annotation_curr(:,3), ...
            roi_location_curr(:,1),roi_location_curr(:,2),roi_location_curr(:,3), roi_annotation_curr(:,1), ...
            'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});

        [rgall, ida, idb] = unique(roi_annotation_curr(:, 2));
        ncell = zeros(length(rgall), 1);
        for j = 1: length(rgall)
            ncell(j) = sum(idb == j);
        end
        roi_table_group{i} = table(rgall, ncell, 'VariableNames', {'name', 'cell_count'});
    end
end
