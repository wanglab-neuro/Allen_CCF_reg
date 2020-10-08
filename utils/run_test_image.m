fnames = dir('*.jpg');
for i = 1: length(fnames)
    fname = fnames(i).name;
    a = imread(fname);
    a = imresize(a, [1400, 2000]);
    pl = cell_detect(a(:, :, 2));
    figure(gcf)
    clf
    hold on
    imagesc(a)
    for j = 1: size(pl, 1)
        plot(pl(j, 2), pl(j, 1), '.r')
    end
    axis tight
    set(gca, 'YDir', 'reverse')
    saveas(gcf, [fname(1: end - 4), '_detected.png'])
end