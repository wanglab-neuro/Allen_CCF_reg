%%% post-correction for the 9 mice %%%
% mice = {'genio_029_genio2_N2c', 'genio024_genio3', 'genio024_genio4', 'Mas018_Mas2_N2c', 'wp007', 'wp023_wp3_N2c'};
% mice = {'M1', 'M2', 'W1'};
mice = {'M1_Re'};
% angles = [-10, 2; -5, 3; -2, -2; -12, 2; -15, -3; -8, 0];
% angles = [-17, 0; 6, 5.8; 0, -1.2];
angles = [17, 0];
% pname = 'Y:\Jinghao\for_others\jun\For coordinate fix_3D rotation_20200623\';
pname = 'Y:\Jinghao\for_others\jaehong\ofofacial_CCF\redo\';
% pname = 'C:\Users\Jinghao\Documents\MATLAB\For coordinate fix_3D rotation_20200623\';
% allmice = {'genio_029_genio2_N2c', 'genio012_genio2', 'genio024_genio3', 'genio024_genio4', 'Mas016_Mas2', 'Mas018_Mas2_N2c', 'wp005', 'wp007', 'wp023_wp3_N2c'};
% allmice = {'M1', 'M2', 'W1'};
allmice = {'M1_Re'};
plall = cell(1, length(allmice));
count = 1;

for i = 1: length(allmice)
    if any(contains(mice, allmice(i)))
        fnames = dir([pname, allmice{i}, filesep]);
        pls = [];
        for j = 1: length(fnames)
            namet = fnames(j).name;
            if contains(namet, 'cell_points') && contains(namet, '.mat')
                load([pname, allmice{i}, filesep, namet])
                
                %%% correct the pl %%%
                ng = length(pl);
                for k = 1: ng
                    if ~isempty(pl{k})
                        angt = angles(count, :);
                        thetay = angt(1);
                        thetaz = angt(2);
                        dp = pl{k}(:, [3, 1, 2]);
                        dp(:, 3) = 800 - dp(:, 3);
                        b = [dp(1, 1), 570, 400];
                        dp = dp - b;
                        Ry = [cosd(thetay), 0, sind(thetay); 0, 1, 0; -sind(thetay), 0, cosd(thetay)];
                        Rz = [cosd(thetaz), -sind(thetaz), 0; sind(thetaz), cosd(thetaz), 0; 0, 0, 1];
                        dpf = dp * Ry' * Rz';
                        dpf = dpf + b;
                        dpf(:, 3) = 800 - dpf(:, 3);
                        pls = [pls; dpf];
                    end
                end
            end
        end
        count = count + 1;
    else
        fnames = dir([pname, allmice{i}, filesep]);
        pls = [];
        for j = 1: length(fnames)
            namet = fnames(j).name;
            if contains(namet, 'cell_points') && contains(namet, '.mat')
                load([pname, allmice{i}, filesep, namet])
                ng = length(pl);
                for k = 1: ng
                    if ~isempty(pl{k})
                        dpf = pl{k}(:, [3, 1, 2]);
                        pls = [pls; dpf];
                    end
                end
            end
        end
    end
    plall{i} = pls(:, [2, 3, 1]);
end

for i = 1: length(plall)
    tmp = plall{i};
    id = tmp(:, 1) >= 0 & tmp(:, 1) <= 1140 & tmp(:, 2) >= 0 & tmp(:, 2) <= 800 & tmp(:, 3) >= 0 & tmp(:, 3) <= 1320;
    plall{i} = plall{i}(id, :);
end

save([pname, 'all_points.mat'], 'plall', 'allmice')

%% get roi table %%
[roi_table, roi_table_group] = summary_table(plall, av, st);
save([pname, 'table_summary_all_mice.mat'], 'roi_table', 'roi_table_group')




