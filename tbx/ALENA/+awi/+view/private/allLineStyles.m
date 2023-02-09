function [linStyle, mkrStyle] = allLineStyles(n)
%allLineStyles Returns a cell array of all the line styles (combinations of
%line style and marker style) available in the default MATLAB installation.
%
% see also LineSpec
%
% Author    : Christopher Szczyglowski
% Username  : cs1414
% Email     : cs1414@bristol.ac.uk
% Timestamp : 27-Jul-2018 10:03:20
%
% Copyright (c) 2018 Christopher Szczyglowski
% All Rights Reserved
%
%
% Revision: 1.0 27-Jul-2018 10:03:20
%   - Initial function :
%

if nargin < 1
    n = [];
else
    if isnan(n) || n == 0
        n = [];
    else        
        validateattributes(n, {'numeric'}, {'positive', 'scalar', 'integer'}, ...
            'allLineStyles', 'n');
    end
end

%Default line and marker styles
lin = {'-', '--', '-.', '.'};
mkr = {'none', 'o', 's', 'd', '^', 'p', 'h', 'v', '>', '<', 'x', '.', '*'};

linStyle = repmat(lin , [numel(mkr), 1]);
mkrStyle = repmat(mkr', [1, numel(lin)]);

linStyle = linStyle(:);
mkrStyle = mkrStyle(:);

%Index?
if ~isempty(n)
    if n > numel(linStyle)
        
        %Pad with existing points
        linStyle = [linStyle ; linStyle(1 : n - numel(linStyle))];
        mkrStyle = [mkrStyle ; mkrStyle(1 : n - numel(mkrStyle))];
        
    else
        
        linStyle = linStyle(1 : n);
        mkrStyle = mkrStyle(1 : n);
        
    end
end

end

