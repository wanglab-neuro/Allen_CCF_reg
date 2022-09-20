function thres = ksdense_filter(ain)
    t1 = normalize(ain);
    t1 = t1(t1 > 0.01);
    edges = linspace(min(t1), max(t1), 101);
    n = histcounts(t1, edges);
    n = n / sum(n);
%     [~, id] = max(n);
%     mx = edges(id);
%     x = edges(1: max(2 * id, 5));
    x = edges(1: end - 1) + diff(edges) / 2;
    f = fit(x(1: 20)', n(1: 20)', 'gauss1', 'lower', [0, -1, 0], 'upper', [100, 1, 1000]);
    nn = feval(f, x');
%     thres = x(find(nn < max(nn) * 0.5, 1));
    [~, id] = find(n - nn > max(nn) * 0.01 & x > x(2), 1);
    thres = x(id);
end
