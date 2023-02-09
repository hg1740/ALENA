function [] = plotremesh(g,p,clr)

y = struct2mat(g,'coord',2);
halfSpan = max(y);
y = y / halfSpan;

idx = 1;

for i = 1:numel(p)
    plot([y(i);y(i + 1)],[p(i).j(idx);p(i).j(idx)],clr);
    hold on
end

