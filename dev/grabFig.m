function varargout = grabFig(varargin)
%grabFig Saves the full figure window of a MATLAB Figure. The saved image
%will include the toolbar and the menubar.
%
% Parameter Inputs:
%   - 'handle' : Handle to the graphics object of type 'matlab.ui.Figure'.
%   - 'target' : Name of the file where the image will be saved to. Must
%                include a valid image file extension. See imwrite for more
%                details.
%   - 'crop'   : Defines the number of pixels on the 
%               [top, bottom, left, right] that will be cropped from the
%               figure.
%
%   N.B. This function uses the FEX function 'screencapture' which is
%   authored by Yair Altman.

p = inputParser;
addParameter(p, 'handle', gcf, @(x)isa(x, 'matlab.ui.Figure'));
addParameter(p, 'target', 'test.png'  , @ischar);
addParameter(p, 'crop'  , [0, 0, 0, 0], @(x)validateattributes(x, {'numeric'}, {'integer', 'numel', 4}));
parse(p, varargin{:});

% %Grab figure handle
% if isempty(p.Results.handle)
%     hFig = gcf;
% else
%     hFig = p.Results.handle;
% end

%Grab crop data
c = p.Results.crop; % [top, bottom, left, right]

%Use Yair Altman's function
im = screencapture('handle', p.Results.handle);

%Remove border from left, bottom & top of image
image = arrayfun(@(i) im(c(1) + 1 : end - c(2), c(3) + 1 : end - c(4), i), 1 : 3, 'Unif', false);
image = cat(3, image{:});

%Write the file
imwrite(image, p.Results.target);

if nargout >  1
    varargout{1} = image;
end

end

