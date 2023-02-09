function val = bool2offon(b, vals)
if nargin < 2
    vals = {'off', 'on'};
end
val = vals{b + 1};
% if b
%     val = 'on';
% else
%     val = 'off';
% end
% end
