function img = multi_otsu(ain)
    nl = 20;
    s = zeros(1, nl);
    level = multithresh(ain, nl);
    for i = 1: nl
    %     s(i) = sum(sum(a0 .* (a0 > level(i))));
        tmp = ain > level(i);
        [l, n] = bwlabeln(tmp);
    %     ss = zeros(1, n);
    %     for j = 1: n
    %         ss(j) = length(find(l == j));
    %     end
        s(i) = n;
    end
    
    if (max(s) - min(s)) / max(s) > 0.1
        lvl = interp1(1: length(level), level, 1: 0.02: length(level));
        ss = interp1(1: length(s), s, 1: 0.02: length(s));
        ss = sgolayfilt(ss, 1, 2 * round(0.1 * length(ss)) + 1);
        
        thres = intense_filter(ss(:), 1);
        id = find(ss < thres, 1);
        img = ain > lvl(id);
    else
        img = ain;
    end
end