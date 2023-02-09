function extractFuselageFromGA(filename)
%extractFuselageFromGA Extracts the fuselage data from a General
%Arrangement (GA) drawing. 

close all

validImageExt = {'.jpeg', '.png'};

%Parse
if nargin < 1
    temp       = strcat('*', validImageExt);
    imageStr   = strjoin(temp, ';');
    imageNam   = strjoin(temp, ', ');
    clear temp
    [file, path] = uigetfile({imageStr,  ...
        ['Image Files (', imageNam, ')']}, 'Pick an image file');
    if isnumeric(path) || isnumeric(file)
        return
    end
    filename = fullfile(path, file);
end
[~, ~, ext] = fileparts(filename);
validatestring(ext, validImageExt, 'extractFuselageFromGA', 'the file extension');

%Open the image and stash graphics handles
hIm = imshow(filename);
hAx = hIm.Parent;
hF  = hAx.Parent;

%Open up a seperate options window

%Want to be able to save a mouseclick....
%   - ButttonDownFcn
end

