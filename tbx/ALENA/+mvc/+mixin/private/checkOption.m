function [val, argout] = checkOption(name, def, varargin)

%Start here
val = def;
argout = varargin;

%Look for name in varargin
idx = find(cellfun(@(x)ischar(x) && strcmp(x, name), varargin));

%Anything ?
if isempty(idx)
    
    %No
    return;
    
end
    
%If we are looking for a logical
if islogical(val)
    
    %Then if the option is found, simply toggle the default
    val = ~val;
    
    %And remove option from return list
    argout(idx) = [];
    
else
    
    %If the option is found, the value must follow
    val = varargin{idx + 1};
    
    %And remove BOTH from return list
    argout(idx + [0,1]) = [];
    
end
