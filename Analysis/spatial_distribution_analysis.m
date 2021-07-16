%% file parse %%
pname = 'Y:\Jinghao\for_others\jun\roi_table\old_coordinates\';
namelist = {'Genioglossus', 'Masseter', 'whisker'};
folders = dir(pname);
fname = cell(size(namelist));
count = 1;
for i = 1: length(folders)
    if any(contains(namelist, folders(i).name))
        fnamet = [pname, folders(i).name, filesep];
        fnames = dir(fnamet);
        fname{count} = [];
        for j = 1: length(fnames)
            if contains(fnames(j).name, '.mat')
                fname{count} = [fname{count}, {[fnamet, fnames(j).name]}];
            end
        end
        count = count + 1;
    end
end

%% read mat file %%
data = fname;
for i = 1: length(fname)
    for j = 1: length(fname{i})
        load(fname{i}{j})
        data{i}{j} = table2array(roi_table{1}(:, 3: 5));
    end
end
% for i = 1: length(data)
%     for j = 1: length(data{i})
%         data{i}{j} = table2array(data{i}{j}(:, 3: 5));
%     end
% end

%% estimate 3d distribution %%
%%% prepare parameters %%%
rgs = [-7.8, 5.4; 0, 8; -5.7, 5.7];
res = 21;
[ap, dv, ml] = ndgrid(linspace(rgs(1, 1), rgs(1, 2), res), linspace(rgs(2, 1), rgs(2, 2), res), linspace(rgs(3, 1), rgs(3, 2), res));
pts = [ap(:), dv(:), ml(:)];

%%% kernel density estimation %%%
dataks = zeros(size(pts, 1), sum(cellfun(@(x) length(x), data)));
count = 1;
for i = 1: length(data)
    for j = 1: length(data{i})
        dataks(:, count) = mvksdensity(data{i}{j}, pts, 'bandwidth', 1);
        count = count + 1;
    end
end

%%% correlation matrix %%%
mtx = norm_inner(dataks', dataks);

t = data{2}{3};
tt = data{2}{4};





