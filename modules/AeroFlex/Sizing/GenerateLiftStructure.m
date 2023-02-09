function data = GenerateLiftStructure(data,Bar_prop,Aircraft_param,Part,Parameters)
%   - C.Szczyglowski (08/10/2019) Creastes the 'Nastran' bulk data objects
%   for the analysis - Could be replaced by my FEModel object.

% Recover the indices
NodeIdx    = Aircraft_param.Planform.NodeIdx;
LENodeIdx  = Aircraft_param.Planform.LENodeIdx;
TENodeIdx  = Aircraft_param.Planform.TENodeIdx;
RBE0Idx    = Aircraft_param.Planform.RBE0Idx;
CBEAMIdx   = Aircraft_param.Planform.CBEAMIdx;
PBEAMIdx   = Aircraft_param.Planform.PBEAMIdx;
CONM2Idx   = Aircraft_param.Planform.CONM2Idx;
SET1Idx    = Aircraft_param.Planform.SET1Idx;
SPLINE1Idx = Aircraft_param.Planform.SPLINE1Idx;
CAEROIdx   = Aircraft_param.Planform.CAEROIdx;

Bar_prop.Mat = data.PartIDs.(Part).Mat;

%% Generate Grid node for carrythrough section
if isfield(Aircraft_param,'CT')
    [data.Node,data.Info,~]  = ...
        NodeCreator(data.Node,data.Optim.(Part).nodex(1), ...
        data.Optim.(Part).nodey(1), data.Optim.(Part).nodez(1) - 2, data.Info, 100000);
end

%% Generate Grid nodes of the beam ends
[data.Node,data.Info,data.PartIDs.(Part).BeamNodes]  = ...
    NodeCreator(data.Node,Bar_prop.nodex, ...
    Bar_prop.nodey, Bar_prop.nodez, data.Info, NodeIdx);

%% Generate the required Aeronodes
[data.RBE0,data.Info,data.PartIDs.(Part).RBE0]  = ...
    RBE0Creator(data.RBE0,data.Info,RBE0Idx,NodeIdx,LENodeIdx,TENodeIdx,Bar_prop.Nsec+1);

%% Generate the leading edge grid nodes
[data.Node,data.Info,data.PartIDs.(Part).LENodes]  = ...
    NodeCreator(data.Node,data.Aero.(Part).LE_ends, ...
    Bar_prop.nodey, Bar_prop.nodez, data.Info, LENodeIdx);

%% Generate the trailing edge grid nodes
[data.Node,data.Info,data.PartIDs.(Part).TENodes]  = ...
    NodeCreator(data.Node,data.Aero.(Part).TE_ends, ...
    Bar_prop.nodey, Bar_prop.nodez, data.Info, TENodeIdx);

if Parameters.BeamType == 1
    %% Generate Bar Connections and orientation vectors
    [data.Beam,data.Info,data.PartIDs.(Part).CBeam]   = ...
        CBEAMCreator(data.Beam,Bar_prop,data.Info,CBEAMIdx,data.Node);
    
    %% Generate Bar properties
    [data.PBeam,data.Info,data.PartIDs.(Part).PBeam]  = ...
        PBEAMCreator(data.PBeam,Bar_prop,data.Info,PBEAMIdx);
    
else
    %% Generate Bar Connections and orientation vectors
    [data.Bar,data.Info,data.PartIDs.(Part).CBar]   = ...
        CBARCreator(data.Bar,Bar_prop,data.Info,CBEAMIdx,data.Node);
    %% Generate Bar properties
    [data.PBar,data.Info,data.PartIDs.(Part).PBar]  = ...
        PBARCreator(data.PBar,Bar_prop,data.Info,PBEAMIdx);
end


%%
data.PartIDs.(Part).Nodes = [data.PartIDs.(Part).BeamNodes,...
    data.PartIDs.(Part).LENodes,data.PartIDs.(Part).TENodes];

if strcmp(Parameters.SplineType,'Surface')
    
    [data.Aero,data.Info,data.PartIDs.(Part).SetIDs] = ...
        SET1Creator(data.PartIDs.(Part).Nodes,data.Aero,data.Info,data.RBE0,SET1Idx);
    
    if Aircraft_param.Planform.SymPlane == 0
        [data.Aero,data.Info,data.PartIDs.(Part).SplineIDs] = ...
            SPLINECreator(data.Aero,data.Aero.(Part),Bar_prop,...
            data.Info,CAEROIdx,SPLINE1Idx);
    end
    
elseif strcmp(Parameters.SplineType,'Beam')
    
    [data.Aero,data.Info,data.PartIDs.(Part).SetIDs] = ...
        SET1CreatorLinear(data.PartIDs.(Part).Nodes,data.Aero,data.Info,data.RBE0,SET1Idx);
    if Aircraft_param.Planform.SIMXZ == 0
        [data.Aero,data.Info,data.PartIDs.(Part).SplineIDs] = ...
            SPLINECreatorLinear(data.Aero,data.Aero.(Part),Bar_prop,...
            data.Info,CAEROIdx,SPLINE1Idx);
    end
end


[data.Mat,data.Info]   = ...
    MAT1Creator(data.Mat,Aircraft_param.Box,data.Info,data.PartIDs.(Part).Mat);

%% GENERATE CONM2s or CONM1s for the Structure mass instead of using the bar/beam density
% Requires having to calculate the inertias about the grid nodes or cg
% positions rotated around to the global frame
if Parameters.ElementMass == 0
    [data.Conm2,data.Info,data.PartIDs.(Part)] = Beam2Conm2(data.Conm2,data.Beam,data.PBeam,data.Mat,data.Node,data.Info,CONM2Idx,1,data.PartIDs.(Part));
    data.Mat.Rho(end) = 0;
end
end