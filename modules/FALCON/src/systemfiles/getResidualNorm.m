function norm_r = getResidualNorm(r)

if iscell(r)
    norm_r = [];
    for ii = 1:length(r)
        norm_r(ii) = norm(r{ii});
    end
    norm_r = norm(norm_r);
else
    norm_r = norm(r);
end