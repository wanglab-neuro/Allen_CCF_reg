%% count cell %%
[fname, pname] = uigetfile('*', 'Select imaging file', 'MultiSelect', 'off');
filefmt = fname(end - 3: end);
fnames = {fname};

%%% load the reference brain and region annotations %%%
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    ppname = mfilename('fullpath');
    ppname = fileparts(ppname);
    load(fullfile(ppname, 'AllenAtlas.mat'))
end

%%% main code %%%
curr_slice = [];
curr_angle = [];
for i = 1: length(fnames)
    %%% generate viewers %%%
    close all
    f = figure('Name','Atlas Viewer');
    fs = figure('Name','Slice Viewer');
    sliceBrowser(fs, pname, f, size(tv), fnames(i));
    f = allenAtlasBrowser(f, tv, av, st, fs, fnames(i), curr_slice, curr_angle);
    waitfor(fs)
    ud = get(f, 'UserData');
    curr_slice = ud.currentSlice;
    curr_angle = ud.currentAngle;
    disp(num2str(curr_slice))
end

%% analysis & reconstruction %%
analyze_cell_count(fname, pname, av, tv, st);
