function pl = cell_detect(a)
    %%% cell enhancing %%%
    a = normalize(double(a));
    a0 = normalize(imgaussfilt(a, 0.5) .^ 2);
    
    nl = 2;
    t = a0;
    tt = imgaussfilt(t, 19);
    for i = 1: nl
        t = abs(t - tt);
        tt = imgaussfilt(t, 49);
    end
    
    t = normalize(t);
    nl = 1;
    for i = 1: nl
        t = normalize(imgaussfilt(double(multi_otsu(t)), 2));
    end
    a1 = t;
%     tt = max(0, a0 - t);
%     [gx, gy] = gradient(a0);
%     a1 = normalize(imgaussfilt(gx .^ 2 + gy .^ 2, 2));
%     a1 = normalize(imgaussfilt(double(multi_otsu(t)), 2));
%     a1t = sgolayfilt(sgolayfilt(double(multi_otsu(a0)), 1, 125, [], 1), 1, 125, [], 2);
%     a1 = normalize(imgaussfilt(a1t .* a0, 1));
    
    nl = 2;
    t = a1;
    tt = imgaussfilt(t, 19);
    for i = 1: nl
        t = abs(t - tt);
        tt = imgaussfilt(t, 49);
    end
    [gx, gy] = gradient(t);
    a2 = normalize(imgaussfilt(gx .^ 2 + gy .^ 2, 2));
    a3 = sigmoid(a2, 100, intense_filter(a2, 0.6));
    
    t = normalize(a3);
    nl = 1;
    for i = 1: nl
        t = normalize(imgaussfilt(double(multi_otsu(t)), 2));
    end
    a4 = normalize(t .* a1);
%     thres = sgolayfilt(sgolayfilt(double(a2), 1, 25, [], 1), 1, 25, [], 2);
%     a3 = normalize(max(0, a2 - thres));
        
%     threst = imgaussfilt(sgolayfilt(sgolayfilt(double(a3), 1, 225, [], 1), 1, 225, [], 2), 149);
%     thresmask = imgaussfilt(sigmoid(threst, 1000, intensity_filter(threst)), 19);
%     thresmask = normalize(threst) > 0.01;
%     a4 = normalize((a3 .* threst));
    
%     thres2 = ksdense_filter(a2);
%     thres2 = intensity_filter(a3);
%     thres2 = 0.04;
    
%     a4 = normalize(sigmoid(a3, 100, thres2));
% %     a4 = a3 .* a2;
%     a5 = a4 .* (a4 > 0.2);
%     
%     a6 = normalize(imgaussfilt(a5, 1));
%     a7 = a6 > intensity_filter(a6);
%     a8 = normalize(a7 .* a6);
% %     a8 = a7 > intensity_filter(a7);
% %     a9 = normalize(a7 .* a8);
%     
% %     thres3 = ksdense_filter(a9);
% %     thres3 = intensity_filter(a7);
% %     a8 = normalize(sigmoid(a7, 100, thres3)) > 0.1;
% %     a9 = a8 .* a7;
%     thres3 = ksdense_filter(a8);
% %     thres2 = intensity_filter(a3);
% %     thres2 = 0.04;
%     
%     a9 = normalize(sigmoid(a8, 100, thres3));
    af = normalize(imgaussfilt(a4 .* (a4 > intense_filter(a4, 0.2)), 1));

    %%% find cells %%%
    mxs = imregionalmax(af);
    [l, n] = bwlabeln(mxs);
    [y, x] = find(mxs > 0);
    pl = [y, x];
end

nr = 20;
a0 = normalize(imgaussfilt(double(a), 2));
s = zeros([size(a0), nr]);
stp = 0.02;
k = sgolayfilt(sgolayfilt(a0, 1, 125, [], 1), 1, 125, [], 2);
t = a0;
for i = 1: nr
%     s(:, :, i) = a0 .* (a0 > k + (i) * stp);
%     s(:, :, i) = a0 .* max(0, (a0 - (k + (i) * stp)));
%     k = sgolayfilt(sgolayfilt(t, 1, 125, [], 1), 1, 125, [], 2);    
%     x = [k(:), ones(size(k(:)))] \ t(:);
%     k = [k(:), ones(size(k(:)))] * x;
%     k = reshape(k, size(t));
    t = max(0, a0 - k - i * stp);
    s(:, :, i) = a0 .* t;
end

%%% otsu %%%
nl = 20;
a = normalize(double(a));
a0 = normalize(imgaussfilt(a, 0.5) .^ 2);
s = zeros(1, nl);
level = multithresh(a0, nl);
for i = 1: nl
%     s(i) = sum(sum(a0 .* (a0 > level(i))));
    tmp = a0 > level(i);
    [l, n] = bwlabeln(tmp);
%     ss = zeros(1, n);
%     for j = 1: n
%         ss(j) = length(find(l == j));
%     end
    s(i) = n;
end
