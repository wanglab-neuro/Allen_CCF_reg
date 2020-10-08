a = normalize(double(a));
b = imopen(a, strel('disk', 5));
a0 = a .* b;
a1 = imgaussfilt(a, 3);
a2 = medfilt2(a1, [29, 29]);
a3 = max(0, a1 - a2);
a4 = a3 .^ 2;
thres = sgolayfilt(sgolayfilt(double(a4), 1, 25, [], 1), 1, 25, [], 2);
a5 = sigmoid(a4, 1000, thres) .* a4;
a6 = imgaussfilt(a5, 3);
af = a6 > intensity_filter(a6);
af = a6 .* af;
af = imgaussfilt(af, 7);

%%% find local max %%%
mxs = imregionalmax(af);
[l, n] = bwlabeln(mxs);



