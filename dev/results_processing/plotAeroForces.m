function plotAeroForces(FEM, AeroForce, AeroPressure)
%plotAeroForces Draws markers for the aerodynamic force on each panel in
%FEM'.

close all

%Collect data
FEMs = flatlist(FEM);

%Set up figure environment
hF     = figure('Name', 'Aero Results');
hAx(1) = subplot(2, 1, 1, 'Parent', hF);
hAx(2) = subplot(2, 1, 2, 'Parent', hF);
axis(hAx, 'equal');
set(hAx, 'NextPlot', 'add');

title(hAx(1), 'Aerodynamic Normal Force');
title(hAx(2), 'Aerodynamic Pressures');

%% Force vector plot

%Draw the'awi.fe.AeroPanel' objects only
AP = [FEMs.AeroPanels];
[hPanel, PanelData] = drawElement(AP, hAx(1));

%Calculate normal vector of each panel
n = awi.fe.AeroPanel.calculateFaceNormals(hPanel);

%Grab coordinates of the centre of the panel
centre = squeeze(PanelData.centre);

%Determine the scale factor for the aerodynamic forces
Fz   = AeroForce.T3;
Fz_  = Fz ./ abs(max(Fz)); %normalise

%Maximum length of the quiver
maxVec = 5;

%Scale the normal vector
n = maxVec .* n .* repmat(Fz_', [1, 3]);

%Plot the quiver plot
hQ = quiver3(hAx(1), ...
    centre(:, 1), centre(:, 2), centre(:, 3), n(:, 1), n(:, 2), n(:, 3));

%% Aerodynamic pressures


hg = patch(hAx(2), ...
    PanelData.coord(:, :, 1), ...
    PanelData.coord(:, :, 2), ...
    PanelData.coord(:, :, 3), ...
    AeroPressure.PressCoeff);
colorbar

end

