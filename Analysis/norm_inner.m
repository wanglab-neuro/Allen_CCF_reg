function prod = norm_inner(a, b) %%% a: ncell X t; b: t X n %%%
    prod = a * b ./ (sqrt(sum(a .^ 2, 2)) * sqrt(sum(b .^ 2, 1)));
end
