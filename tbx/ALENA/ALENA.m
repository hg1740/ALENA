function [hF, Model] = ALENA(varargin)
%ALENA Executes the Aircraft Loads Environment for Nonlinear Aeroelastics
%(ALENA) framework.
%
% Inputs:
%   * Required: 
%       - 
%   * Optional:
%       - 'Model' -> An instance of an aircraft model.
%   * Parameters:
%       - 'Mode' -> 'Active' The framework will respond to changes in any
%                   'SetObservable' properties by establishing a range of
%                   listeners for properties, events, etc.
%                -> 'Passive' The framework will not respond to any changes
%                   in observable properties and no listeners will be
%                   generated.
%
% (c) University of Bristol, 2018
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
% TODO - Check dependencies

%Configure the workspace 
%   - Check version of MATLAB
checkRelease('2016b');
%   - Check GUI Layout Toolbox is installed

%Create the splash page
hAl = ALENA_splash_page;

%Add all folders to path?
%addpath(genpath(fullfile(pwd, 'UoB')));
%addpath(genpath(fullfile(pwd, 'FEX')));
%gltDir = fullfile(pwd, 'GLT'); %Folder containing the GUI Layout Toolbox
%if isfolder(gltDir)
%    addpath(genpath(gltDir));
%end

% %Parse inputs
% p = parseInputs(varargin);

%Pass on inputs to the 'awi.view.Framework' class.
Fmwrk = awi.view.Framework( ...
    awi.model.Framework('EntryPointFcn', mfilename, ... % so that 'deploy' knows where to start
    varargin{:}));

%Return the function handle
hF = Fmwrk.Parent;

%Return the model belonging to the ALENA framework
Model = Fmwrk.Model;

%Close the ALENA splash page
close(hAl);

end

function checkRelease(requiredVer)
%checkRelease Checks that the version of MATLAB meets the minimum
%requirments.

%Which year? Is it 'a' or 'b'?
yr   = str2double(requiredVer(isstrprop(requiredVer, 'digit')));
ver_ = requiredVer(isstrprop(requiredVer, 'alpha'));

%Check we are running minumum version or higher
matlabVersion = version('-release');
releaseNum    = sscanf(matlabVersion, '%f');
releaseVer    = matlabVersion(isstrprop(matlabVersion, 'alpha'));
if releaseNum < yr
    error('ALENA requires MATLAB release %s or higher.', requiredVer);
end
if releaseNum == yr && releaseVer < ver_
    error('ALENA requires MATLAB release %s or higher.', requiredVer);
end

end

function hAl = ALENA_splash_page
%ALENA_splash_page Generates a 'splash-page' for the ALENA Framework so 
%that the user knows the code is running.

%Path to the ALENA logo
%   - Image stored in the 'ico' folder in the 'awi.model' package 
dn     = fullfile(fileparts(which('awi.model.Framework')), 'ico');
impath = fullfile(dn, 'ALENA_logo.png');

%Set the figure in the centre of the primary monitor
mp  = get(0, 'MonitorPositions');
w   = 500; %[px], width of the window
h   = 200; %[px], height of the window
pos = [ ...
    mp(1, 1) + mp(1, 3) / 2 - w / 2, ...
    mp(1, 2) + mp(1, 4) / 2 - h / 2, ...
    w, h];

%Generate the figure and add an axes - Hide for the moment
hAl = figure( ...
    'Name'       , 'ALENA' , ...
    'NumberTitle', 'off'   , ...
    'MenuBar'    , 'none'  , ...
    'ToolBar'    , 'none'  , ...
    'Units'      , 'pixels', ...
    'Position'   , pos     , ...
    'Visible'    , 'off'    , ...
    'WindowStyle', 'normal' , ...
    'Resize'     , 'off');
hAx = axes( ...
    'Parent'  , hAl         , ...
    'Visible' , 'off'       , ...
    'Units'   , 'normalized', ...
    'Position', [0 0 1 1]   , ...
    'XTick'   , []          , ...
    'YTick'   , []);

%Import the image and add it to the axes
im  = imread(impath);
hIm = image(im, 'Parent', hAx);

%Resize the figure in case there is a mismatch between the image & figure
figpos       = hAl.Position;
figpos(3:4)  = [hIm.XData(2) hIm.YData(2)];
hAl.Position = figpos;

%Force the axes to have no ticks
set(hAx, 'XTick', [], 'YTick', []);
set(hAx, 'Visible', 'off');
hAl.Visible  = 'on';

movegui(hAl, 'center')

end

% function p = parseInputs(varargin)
% %parseInputs Parses the inputs for the function 'AWI' and returns an
% %inputParser object 'p'.
% 
% p = inputParser;
% addOptional(p, 'Model', []);
% parse(p, 'Model', varargin);
% 
% 
% end

