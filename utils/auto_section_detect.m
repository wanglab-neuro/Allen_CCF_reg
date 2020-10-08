%%% initial align %%%
[nf, nh, nw] = size(tv);
a = imresize(a, [nh, nw]);
img = normalize(double(a(:, :, 3)));
l = img > 0.1;
hmin = find(max(l, [], 2), 1);
hmax = find(max(l, [], 2), 1, 'last');
wmin = find(max(l, [], 1), 1);
wmax = find(max(l, [], 1), 1, 'last');
imcur = img .* l;
imcur = imcur(hmin: hmax, wmin: wmax);
nh = 400;
nw = 570;
imcur =  imresize(imcur, [nh, nw]);
imcur = imgaussfilt(imcur, 5);

%%% hierarchical first section location detection %%%
h1stp = 10;
nhier = ceil(nf / h1stp);
ngl = 3;
stride = 10;
scr = cell(nhier, ngl);
% scr = zeros(nhier, ngl);
inds = [1: h1stp: nf, nf];
dats = zeros(nh, nw, nhier);
datsa = zeros(nh, nw, nhier);
for i = 1: nhier
    tmp = normalize(double(squeeze(mean(tv(inds(i): inds(i + 1), :, :), 1))));
    dats(:, :, i) =  imresize(tmp, [400, 570]);
    tmp = double(squeeze(av(inds(i), :, :)));
    datsa(:, :, i) =  imresize(tmp, [400, 570], 'nearest');
end

tic
parfor i = 1: length(inds) - 1
%     tmp = normalize(dats(:, :, i));
    imrefa = datsa(:, :, i);
    imref = dats(:, :, i);
%     imref = tmp .* double(tmp > 0.01);
    imref = imgaussfilt(imref, 5);
    tmp = imcur;
    for j = 1: ngl % scale %
        [nth, ntw] = size(tmp);
        hid = 1: stride: nh - nth + 1;
        wid = 1: stride: nw - ntw + 1;
        if hid(end) ~= nh - nth + 1
            hid = [hid, nh - nth + 1];
        end
        if wid(end) ~= nw - ntw + 1
            wid = [wid, nw - ntw + 1];
        end
        scr{i, j} = zeros(length(hid), length(wid));
        for k = hid
            for l = wid
                try
                    imgcur = normalize(tmp);
                    imgrefa = imrefa(k: k + nth - 1, l: l + ntw - 1);
                    imgref = imref(k: k + nth - 1, l: l + ntw - 1);
                    idt = unique(imgrefa(:) .* (imgcur(:) > 0.1));
                    temp = false(size(imgrefa));
                    for kk = 1: length(idt)
                        temp = temp | imgrefa == idt(kk);
                    end
                    imgref = imgref .* (temp > 0);
                    scr{i, j}(k, l) = norm_inner(imgref(:)', imgcur(:));
%                     scr{i, j}(k, l) = corr(imgref(:), imgcur(:));
                catch
                    scr{i, j}(k, l) = 0;
                end
            end
        end
        tmp = impyramid(tmp, 'reduce');
    end
    disp(num2str(i))
end
toc

%%% second finer location %%%
scrm = cellfun(@(x) max(x(:)), scr);
[y, x] = find(scrm == max(scrm(:)));
h2stp = 5;
idh = (y - h2stp) * h1stp + 1: (y + h2stp) * h1stp;
dats = zeros(nh, nw, length(idh));
for i = 1: length(idh)
    dats(:, :, i) = normalize(double(squeeze(mean(tv(inds(i): inds(i + 1), :, :), 1))));
end

parfor i = 1: length(idh)
    tmp = normalize(dats(:, :, i));
%     denomref = imgaussfilt(TVL1denoise(tmp, 0.2, 10), 3);
%     denomref = imgaussfilt(tmp, 9);
    imref = tmp .* double(tmp > 0.01);
    imref = imgaussfilt(imref, 5);
%     imref = normalize(feature2_comp(tmp, [], [], denomref / 1.2));
%     imref = anidenoise(imref, 0, 0, 50, 0.2, 2);
    tmp = imcur;
    for j = 1: ngl % scale %
        [nth, ntw] = size(tmp);
        hid = 1: stride: nh - nth + 1;
        wid = 1: stride: nw - ntw + 1;
        if hid(end) ~= nh - nth + 1
            hid = [hid, nh - nth + 1];
        end
        if wid(end) ~= nw - ntw + 1
            hid = [wid, nw - ntw + 1];
        end
        scr{i, j} = zeros(length(hid), length(wid));
        for k = hid
            for l = wid
                try
                    imgcur = tmp;
                    imgref = imref(k: k + nth - 1, l: l + ntw - 1);
%                     img = klt2_reg(imgref, imgcur);
                    scr{i, j}(k, l) = norm_inner(imgref(:)', imgcur(:));
                catch
                    scr{i, j}(k, l) = 0;
                end
            end
        end
%         t = conv2(tmp, imref, 'same');
%         scr(i, j) = max(max(t));
        tmp = impyramid(tmp, 'reduce');
    end
    disp(num2str(i))
end

