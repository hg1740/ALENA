function beam_model = AWI2AeroFlex(beam_model,LC)

% Definition needed
beam_model.SolParam.Aeroelastic = 1;
beam_model.SolParam.Trim = 1;

%% Pre-process some of the workspace

% Setup the vortex lattice mesh
[beam_model.Aero.lattice_vlm, beam_model.Aero.ref] = ...
    vlm_setup(1, beam_model.Aero.geo, beam_model.Aero.state, beam_model.Aero.ref);

% Sort some of the arrays 
beam_model = sortArrays(beam_model);

[beam_model.Aero,beam_model.Info] = ...
    sortINTERPArray(beam_model.Aero,beam_model.Param,beam_model.Info);

%% Trim analysis module

% Get a weight breakdown
beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,[],beam_model.Bar,beam_model.Beam,beam_model.Info);

% Let's deal withe some of the load cases here (This should be done outside
% to preserve the code)
idx = 1;

beam_model.SolParam.trimI = idx;
beam_model.SolParam.NGusts = 0;
beam_model.Param.GRAV = [0,0,9.81]';
beam_model.Param.SUBCASE{1}.ID = 1;
beam_model.SolParam.WingAngle = 0;

%% LOADCASE CONVERSION SCRIPT

LoadCases.M         = LC.Mach;
LoadCases.Alt       = LC.Altitude;
LoadCases.URDD3     = LC.LoadFactor*9.81;

% Assign Control Surface name, dof, and deflection

% Search for the control surface names and deflections
lc_cs_name   = {LC.ControlSurfaces.Name};

% Ignore flaps and slats for the time being
CsIdx = and(cellfun(@isempty,regexpi(lc_cs_name,'fla')),cellfun(@isempty,regexpi(lc_cs_name,'sla')));

% Update the control surface list to only include what we need
Cs      = LC.ControlSurfaces(CsIdx);
CsDef   = LC.CsDeflections.Value(CsIdx);
CsType  = LC.CsDeflections.Type(CsIdx);

% Which one of those control surfaces are fixed?
cs_fixed_idx = strcmpi(CsType,'fixed');

% Isolate the fixed control surfaces
FixedCS = Cs(cs_fixed_idx);

% Find out if any have negative span (this fudge is needed as the control
% surface deflection in the framework is hinge independent
Yend = [];
for i = 1:numel(FixedCS)
    Yend = [Yend,FixedCS(i).Coords.LE.Y(2)];
end

% Which have negative span
NegSpan = find(Yend<0);

% Update the Loadcase structure needed by AeroFlex ... this specifies fixed
% control surfaces
LoadCases.CS.Label = {Cs(cs_fixed_idx).Name};
LoadCases.CS.Value = CsDef(cs_fixed_idx)';

% Account for the negative span
LoadCases.CS.Value(NegSpan) = -LoadCases.CS.Value(NegSpan);

LoadCases.PayloadCG = LC.CgMac;

beam_model.Aero.lattice_vlm.Control.Name = beam_model.Aero.Control.Name;

beam_model.Aero.geo.nc = numel(beam_model.Aero.lattice_vlm.Control.Name);
beam_model.twist_corr = 1;
beam_model.camber_corr = 1;
beam_model.Param.G = LoadCases.URDD3;

%%
% Setup the load cases
beam_model = SetupLoadCases(beam_model,idx,LoadCases(idx));
beam_model.Aero.state.SIMXZ = 0;

% Find the Wing CAERO cards

IDType = {beam_model.PartId.Type};
IDPart = {beam_model.PartId.Part};

WingPartIdx = and(ismember(IDType,'CAERO'),ismember(IDPart,'StbdWing'));

wingidx = [];
for i = 1:length(beam_model.PartId(WingPartIdx).data)
    wingidx = [wingidx,find(beam_model.Aero.ID == beam_model.PartId(WingPartIdx).data(i))];
end

% Calculate the aerodynamic centres and store somewhere.
[mac, mac_LE, mac_AC] = CalculateMAC(beam_model.Aero,wingidx);

beam_model.Aero.ref.MAC       = mac;
beam_model.Aero.ref.MAC_LE_x  = mac_LE;
beam_model.Aero.ref.MAC_ac    = mac_AC;

% TODO - GET RID OF THIS HARDCODED LINE!!!
Payload_Idx = find(beam_model.Conm2.ID == 1002);

% Get the current mass of the beam model
CurrentMass = beam_model.WB.MCG(1);

% Recall the target Aircraft Mass
TargetMass = LC.AcMass;

% Calculate the delta needed to get the correct mass
DeltaMass = TargetMass - CurrentMass;

% Adjust the payload mass
beam_model.Conm2.M(1:3,1:3,Payload_Idx) = beam_model.Conm2.M(1,1,Payload_Idx) + DeltaMass;

beam_model = adjustPayloadMass(beam_model,LoadCases(idx),Payload_Idx);

end