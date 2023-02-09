function filename = getfile(filespec, dlg)
%getfile Thin wrapper around the in-built function 'uigetfile'. 
%
% Allows specification of the file specification and the modal box title.

if nargin < 1
    filespec = {'*', 'All Files'};
end
if nargin < 2
    dlg = 'Select a file';
end

[filename, filepath] = uigetfile(filespec, dlg) ;
if isnumeric(filename) && isnumeric(filepath)
    filename = [];
    return
else
    filename = fullfile(filepath, filename);
end

end

