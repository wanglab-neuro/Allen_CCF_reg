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
    
    nl = 2;
    t = a1;
    tt = imgaussfilt(t, 19);
    for i = 1: nl
        t = abs(t - tt);
        tt = imgaussfilt(t, 49);
    end
    [gx, gy] = gradient(t);
%     t1 = normalize(sgolayfilt(sgolayfilt(a1, 1, 125, [], 1), 1, 125, [], 2));
%     a2 = normalize(imgaussfilt(gx .^ 2 + gy .^ 2, 2) .* t1 .^ 0.2);%imgaussfilt(a1, 199) .^ 0.2);
    a2 = normalize(imgaussfilt(gx .^ 2 + gy .^ 2, 2));%imgaussfilt(a1, 199) .^ 0.2);
    a3 = sigmoid(a2, 1000, intense_filter(a2, 1));
    
    t = normalize(a3);
    nl = 1;
    for i = 1: nl
        t = normalize(imgaussfilt(double(multi_otsu(t)), 2));
    end
    t1 = imgaussfilt(a3, 199);
    a4 = normalize(t .* a1 .* a0 .* t1);
    af = normalize(imgaussfilt(a4 .* (a4 > intense_filter(a4, 0.2)), 2));

    %%% find cells %%%
    mxs = imregionalmax(af);
    [l, n] = bwlabeln(mxs);
    [y, x] = find(mxs > 0);
    pl = {[y, x]};
    
    %%% visualize %%%
    figure(gcf)
    hold on
    for i = 1: length(y)
        plot(x(i), y(i), '.r')
    end
    hold off
end
