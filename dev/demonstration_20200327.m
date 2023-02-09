%% demonstration_20200327
% Demonstration for Airbus personnel.
%   - Basic scripting 
%   - FEM conversion
%   - Parameteric modelling
%   - Exporting FEM

close all

%What have we got?
hF  = figure('Name', Alena.Name);
hAx = axes('Parent', hF, 'NextPlot', 'add');
draw(Alena.Aircraft, hAx);
axis(hAx, 'equal');

%% Basic scripting
%   - Finding parts of the model
%   - Plotting the geometry

%Can find objects based on their attributes (supports any property)
LS        = findall(Alena, 'Type', 'LiftingSurface');
Fuse      = findall(Alena, 'Name', 'Fuselage'      , '-and', 'Type', 'BluffBody');
StbdStrut = findall(Alena, 'Type', 'LiftingSurface', '-and', 'Name', 'Starboard Strut');
Spars     = findall(Alena, 'Type', 'Spar');

%Can also use conditional statements
%   - Same syntax as built-in MATLAB 'findobj' function
SW  = findall(Alena, 'Type', 'LiftingSurface', '-and', 'Name', 'Starboard Wing');
hF  = figure('Name', 'Starboard Wing and Truss Structure');
hAx = axes('Parent', hF, 'NextPlot', 'add');
draw([SW, StbdStrut], hAx);

hF  = figure('Name', 'Starboard Truss');
hAx = axes('Parent', hF, 'NextPlot', 'add');
draw(StbdStrut, hAx);

%% Modifying the object properties

%Set the strut to join the wing at the tip
StbdStrut.TipNormPos = 1;
StbdStrut.TipXOffset = 0;
build(StbdStrut);
draw([SW, StbdStrut]);

%% Converting to FEM representation

%Can make individual FEMs based on specific objects
FEM_Truss = convertToFE(StbdStrut);
FEM_Wing  = convertToFE(SW);
FEM       = convertToFE(Alena.Aircraft);

draw(FEM_Truss);
draw(FEM_Wing);
draw(FEM);

%% Exporting the FEM

export(FEM);


