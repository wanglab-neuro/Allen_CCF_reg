function analyze_cell_count(fname, pname, av, tv, st)
    filefmt = fname(end - 3: end);
    fnamest = dir([pname, ['*', filefmt]]);
    fnames = cell(1, length(fnamest));
    for i = 1: length(fnamest)
        fnames{i} = fnamest(i).name;
    end

    %% analysis & summary %%
    plt = cell(1, length(fnames));
    for i = 1: length(fnames)
        load([pname, 'cell_points_', fnames{i}(1: end - 4), '.mat'])
        plt{i} = pl;
    end

    ngcell = max(cellfun(@length, plt));
    pls = cell(ngcell, 1);
    for i = 1: length(pls)
        pls{i} = [];
    end

    for i = 1: length(plt)
        for j = 1: length(plt{i})
            pls{j} = [pls{j}; plt{i}{j}];
        end
    end

    %%% get separate roi table %%%
    for i = 1: length(pls)
        tmp = pls{i};
        id = tmp(:, 1) >= 0 & tmp(:, 1) <= 1140 & tmp(:, 2) >= 0 & tmp(:, 2) <= 800 & tmp(:, 3) >= 0 & tmp(:, 3) <= 1320;
        pls{i} = pls{i}(id, :);
    end
    
    [roi_table, roi_table_group] = summary_table(pls, av, st);
    roi_table{1}.ML_location = -1.*roi_table{1}.ML_location % added 20210330
    
    
    save([pname, 'table_summary', fnames{1}(1: end - 4), '.mat'], 'roi_table', 'roi_table_group')

    %%
%     %% reconstruction %% % remove this later
%     %%% get 3-perspective map %%%
%     cls = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .5 0 .6; 1 .6 0];
%     figure(1)
%     clf
%     %%% coronal %%%
%     subplot(1, 3, 1)
%     islice = unique(round(cell2mat(cellfun(@(x) squeeze(x(:, 3)), pls, 'uniformoutput', false))));
%     imagesc(squeeze(mean(tv(islice, :, :), 1)))
%     hold on
%     for ii = 1: length(pls)
%         plot(pls{ii}(:, 1), pls{ii}(:, 2), '.', 'color', cls(ii, :), 'markersize', 8)
%     end
%     hold off
% 
%     %%% tangential %%%
%     subplot(1, 3, 2)
%     islice = unique(round(cell2mat(cellfun(@(x) squeeze(x(:, 2)), pls, 'uniformoutput', false))));
%     imagesc(squeeze(mean(tv(:, islice, :), 2))')
%     hold on
%     for ii = 1: length(pls)
%         plot(pls{ii}(:, 3), pls{ii}(:, 1), '.', 'color', cls(ii, :), 'markersize', 8)
%     end
%     hold off
% 
%     %%% sagittal %%%
%     subplot(1, 3, 3)
%     islice = unique(round(cell2mat(cellfun(@(x) squeeze(x(:, 1)), pls, 'uniformoutput', false))));
%     imagesc(squeeze(mean(tv(:, :, islice), 3))')
%     hold on
%     for ii = 1: length(pls)
%         plot(pls{ii}(:, 3), pls{ii}(:, 2), '.', 'color', cls(ii, :), 'markersize', 8)
%     end
%     hold off
%     print([pname, 'three_view'], '-dtiff')
%     savefig([pname, 'three_view.fig'])
end