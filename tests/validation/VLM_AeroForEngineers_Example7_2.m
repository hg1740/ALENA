function VLM_AeroForEngineers_Example7_2
%VLM_AeroForEngineers_Example7_2 Runs a VLM analysis using the geometry
%from Example 7.2, pp301 of Aerodynamics for Engineers, Third Edition, 
%Bertin & Smith.
%
% Description
%   - The wing has an aspect ratio of 5 with a leading edge sweep of 45 deg
%   and no taper (i.e. constant chord).
%   - This unit test will run the following VLM solvers:
%       * 
%
% Outputs

close all

%% Parameters

%Independent
b     = 10;     %[m], span
swpLE = 45;     %[deg], leading edge sweep

%Ananlysis Data
RunData.U_FreeStream = [30, 0, 0];  %[m/s], Freestream velocity components in global system
RunData.AoA          = 5;           %[deg], Angle of attack in degrees

%Dependent
s           = b/2;       %[m], semi-span
c           = 0.2 * b;   %[m], chord
panelWidth  = s / 4;     %[m], aero panel width (in span direction)
panelLength = c;         %[m], aero panel length (in flow direction)
panelAR     = panelWidth / panelLength;    

%% Make the geometry model

%Use the 'Sweep Set' (sSet) parameterisation scheme
Wing = awi.model.LiftingSurface;
Wing.ActiveSet   = 'sSet';
Wing.SpanVector  = 'Y';
Wing.Span        = s;
Wing.RootChord   = c;
Wing.LESweep     = [swpLE, swpLE];
Wing.LESweep_eta = [0, 1];
Wing.TESweep     = [swpLE, swpLE];
Wing.TESweep_eta = [0, 1];
build(Wing);

%Plot it to check
hF = figure('Name', 'Geometry Model');
hAx = axes('Parent', hF, 'NextPlot', 'add');
drawElement(Wing, hAx);

%% Convert to FE
%Set up the FE properties
Wing.AeroPanelLength = panelLength;
Wing.AeroPanelAR     = panelAR;

%Convert to FE
FEM = convertToFE(Wing);

%Draw to check
hF(2)  = figure('Name', 'Analysis Model');
hAx(2) = axes('Parent', hF(2), 'NextPlot', 'add');
drawElement(FEM.AeroPanels, hAx(2));   %To draw the aero panels only
%draw(FEM); %To draw the whole model

%Format axes
axis(hAx, 'equal');
set([hAx.XLabel], 'String'  , 'X [m]');
set([hAx.YLabel], 'String'  , 'Y [m]');
set([hAx.ZLabel], 'String'  , 'Z [m]');
set([hAx.Title] , {'String'}, {'Geometry Model' ; 'Analysis Model'});

%% Run VLM
VLMResults = awi.methods.vlm.hs_vlm(FEM.AeroPanels, RunData);

end

