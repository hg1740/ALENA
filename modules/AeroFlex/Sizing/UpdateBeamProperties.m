function data = UpdateBeamProperties(data,Bar_prop,Aircraft_param,Part,param)

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

Bar_prop.Mat = data.PartIDs.(Part).Mat;

PBeamID = find(data.PBeam.ID == PBEAMIdx+1);

if param.BeamType == 0
    data.Info.npbar = PBeamID - 1;
    
    %% Generate Bar properties
    [data.PBar,data.Info]  = ...
        PBARCreator(data.PBar,Bar_prop.Wing,data.Info,110000);
    
    if ~isempty(Aircraft_param.HTP)
        [data.PBar,data.Info]  = ...
            PBARCreator(data.PBar,Bar_prop.HTP,data.Info,120000);
    end
    
    if ~isempty(Aircraft_param.VTP)
        [data.PBar,data.Info]  = ...
            PBARCreator(data.PBar,Bar_prop.VTP,data.Info,130000);
    end
    
    [data.PBar,data.Info]     = ...
        sortPBarArray(data.PBar,data.Mat,data.Param,data.Info);
    
    data.Bar = set_cbar_database(data.Info.nbar, data.Bar, data.PBar, data.Mat, data.Node, []);
    
else
    
    % Fudge to correctly invoke the beam property
    
    oldpbeam = data.Info.npbeam;
    oldMat   = data.PBeam.Mat;
    
    data.Info.npbeam = PBeamID - 1;
    
    %% Generate Bar properties
    [data.PBeam,data.Info,data.PartIDs.(Part).PBeam]  = ...
        PBEAMCreator(data.PBeam,Bar_prop,data.Info,PBEAMIdx);
    
    if param.ElementMass == 0
        [data.Conm2,data.Info,data.PartIDs.(Part)] = Beam2Conm2(data.Conm2,data.Beam,data.PBeam,data.Mat,data.Node,data.Info,CONM2Idx,1,data.PartIDs.(Part));
        
        data.Mat.Rho(end) = 0;
        
        [data.ConM,data.Info] = ...
            sortCONMArray(data.Conm1,data.Conm2,data.Node,data.Coord,data.Param,data.Info);
    end
    
    
    data.Info.npbeam = oldpbeam;
end

data.PBeam.Mat = oldMat;
data.WB = WBnCG(data.Node,data.Param,data.ConM,data.WB,data.Bar,data.Beam,data.Info);

end