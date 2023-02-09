function val = customCodeFolders(val)

%Rather tedious but need to avoid losing original path settings,
% even after "clear all"
mlock;
persistent VAL;
if isempty(VAL)
    VAL = {};
end

%Do what ?
if nargin == 0
    
    %Retrieve
    val = VAL;
    
else
    
    %Assign
    VAL = val;
    
end

end
