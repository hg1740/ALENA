function Results = hs_vlm(AeroPanels, runData)
%hs_vlm Calculates the circulation distribution using the HORSESHOE VLM
%method for a set of panels ('AeroPanels') using the aerodynamic state data
%in 'runData'.
%
% Inputs
%   - 'AeroPanels' : Array of 'awi.fe.AeroPanel' objects describing the
%   aerodynamic geometry of the panels (not the vortex elements).
%   - 'RunData' : Structure containing the aerodynamic state data. 
%       TODO - Update to an object.
%
% Outputs
%
%

%Draw?
%hF = figure;
%hAx = axes('Parent', hF, 'NextPlot', 'add');
%drawElement(AeroPanels, hAx);

% %Get panel geometry
% AeroPanelGeometry = definePanels(AeroPanels);

%Calculate geometry for VLM mesh
AeroMesh = defineAeroMeshGeometry(AeroPanels);

end

