function beam_model = FameObj2BeamModel(Fame)

%TODO:  This function only takes converts a subset of the object and
%populates a subset of the beam model. More work has to be done on this to
%get this running in Aeroflex.

beam_model = InitRead;

% Deal the Grid nodes
beam_model.Node.ID      = [Fame.Grid.id]';
beam_model.Node.CS      = [Fame.Grid.cs]';
beam_model.Node.Coord   = vertcat(Fame.Grid.coord);
beam_model.Node.CD      = [Fame.Grid.cd]';
beam_model.Info.ngrid   = length(beam_model.Node.ID);

% Deal the Beam Connectivity
beam_model.Beam.ID      = [Fame.Cbar.id]';
beam_model.Beam.PID     = [Fame.Cbar.pid]'; 
beam_model.Beam.Conn(:,[1,2])   = vertcat(Fame.Cbar.conn);
beam_model.Beam.Conn(:,3)       = 0;
beam_model.Beam.Orient  = vertcat(Fame.Cbar.orient);
beam_model.Beam.Offset  = zeros(numel(beam_model.Beam.ID),9);
beam_model.Beam.beamg0  = zeros(numel(beam_model.Beam.ID),1);
beam_model.Beam.OffsetT = ones(numel(beam_model.Beam.ID),1); % GGG
beam_model.Info.nbeam   = length(beam_model.Beam.ID);

% Deal the Beam Properties
beam_model.PBeam.ID     = [Fame.Pbar.id]';
beam_model.PBeam.Mat    = [Fame.Pbar.mat]';      
for i = 1:numel(Fame.Pbar)
    beam_model.PBeam.A(i).data   = [Fame.Pbar(i).a,Fame.Pbar(i).a];
    beam_model.PBeam.I(i).data   = [Fame.Pbar(i).i',Fame.Pbar(i).i'];
    beam_model.PBeam.J(i).data   = [Fame.Pbar(i).j,Fame.Pbar(i).j];
    beam_model.PBeam.RhoNS(i).data  = [0,0];
    beam_model.PBeam.Kshear(i,:)    = [0,0];
    beam_model.PBeam.X_L(i).data    = [0,1];
    beam_model.PBeam.NSI(i,:)       = zeros(1,2);
    beam_model.PBeam.NSCG(i,:)       = zeros(1,4);
    beam_model.PBeam.NA(i,:)       = zeros(1,4);
    beam_model.PBeam.DA(:, :, i) = zeros(1,2);
    beam_model.PBeam.DI(:, :, i) = zeros(3,2);
    beam_model.PBeam.DJ(:, :, i) =  zeros(1,2);
    beam_model.PBeam.DRhoNS(:, :, i) =  zeros(1,2);
    
    beam_model.PBeam.DA(1, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.A(i).data,	  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    beam_model.PBeam.DI(1, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.I(i).data(1,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    beam_model.PBeam.DI(2, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.I(i).data(2,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    beam_model.PBeam.DI(3, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.I(i).data(3,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    beam_model.PBeam.DJ(1, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.J(i).data,	  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    beam_model.PBeam.DRhoNS(1, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
        beam_model.PBeam.RhoNS(i).data,  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
%     beam_model.PBeam.DNAOff(1, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
%         beam_model.PBeam.NAOff(i).data(1,:),  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
%     beam_model.PBeam.DNAOff(2, 1:2, i) = interp1(beam_model.PBeam.X_L(i).data,...
%         beam_model.PBeam.NAOff(i).data(2,:),  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
end
beam_model.Info.npbeam   = length(beam_model.PBeam.ID);

% Deal the Aerodynamic Properties
beam_model.Aero.ID          = [Fame.Caero.id]';
beam_model.Aero.CP          = [Fame.Caero.cid]';
beam_model.Aero.geo.dihed   = pi/180*[Fame.Caero.dih]';
beam_model.Aero.geo.ny      = [Fame.Caero.ny]';
beam_model.Aero.geo.nx      = [Fame.Caero.nx]';
beam_model.Aero.geo.flapped = [Fame.Caero.flapped]';
beam_model.Aero.geo.fnx     = [Fame.Caero.fnx]';
fc = vertcat(Fame.Caero.fc);
beam_model.Aero.geo.fc(:,1,1)      = fc(:,1);
beam_model.Aero.geo.fc(:,1,2)      = fc(:,2);
beam_model.Aero.geo.startx  = [Fame.Caero.startX]';
beam_model.Aero.geo.starty  = [Fame.Caero.startY]';
beam_model.Aero.geo.startz  = [Fame.Caero.startZ]';
beam_model.Aero.geo.c       = [Fame.Caero.c]';
beam_model.Aero.geo.b       = [Fame.Caero.b]';
beam_model.Aero.geo.T       = [Fame.Caero.t]';
beam_model.Aero.geo.SW      = pi/180*[Fame.Caero.sw]';
beam_model.Aero.geo.TW(:,1,1)      = pi/180*[Fame.Caero.twist1]';
beam_model.Aero.geo.TW(:,1,2)      = pi/180*[Fame.Caero.twist2]';
beam_model.Aero.geo.nelem   = ones(size(beam_model.Aero.geo.startx));
beam_model.Aero.geo.fsym    = zeros(size(beam_model.Aero.geo.startx));
beam_model.Aero.geo.foil(:,1,1) = {Fame.Caero.rootAirfoil}';
beam_model.Aero.geo.foil(:,1,2) = {Fame.Caero.tipAirfoil}';
beam_model.Aero.geo.symetric    = [Fame.Caero.symetric]';
beam_model.Aero.geo.flap_vector = zeros(size(beam_model.Aero.geo.startx));
beam_model.Aero.geo.meshtype    = [Fame.Caero.meshType]';

% List all of the individual control surfaces
CsNames = {Fame.Caero.csId};
% List the parent control surface
CsType  = {Fame.Caero.csType};
% List control surface links
CsLink  = [Fame.Caero.csLink];
% Remove any empty entries
CsIdx   = ~cellfun(@isempty,CsNames);
% Overwrite the list of control surfaces
CsType  = CsType(CsIdx);
CsNames = CsNames(CsIdx);
CsLink  = CsLink(CsIdx);
beam_model.Aero.Control.Name = CsNames;

% Let's separate out the port and starboad cs
for i = 1:numel(CsType)
    CsType{i} = [CsType{i} num2str(CsLink(i))];
end

UniqueCS = unique(CsType,'Stable');

% Let's check whether or not any of the control surfaces need to be slaved
beam_model.Aero.Trim.Link.ID   = [];
beam_model.Aero.Trim.Link.Master = {};
beam_model.Aero.Trim.Link.Slave  = {};
beam_model.Aero.Trim.Link.Coeff  = [];
LinkCount = 0;
for i = 1:numel(UniqueCS)
    
    RepeatedCS = ismember(CsType,UniqueCS{i});
    RepeatedCSIdx = find(RepeatedCS);
    
    for j = 2:numel(RepeatedCSIdx)
        LinkCount = LinkCount + 1;
        beam_model.Aero.Trim.Link.ID     = [beam_model.Aero.Trim.Link.ID,LinkCount];
        beam_model.Aero.Trim.Link.Master = [beam_model.Aero.Trim.Link.Master,CsNames(RepeatedCSIdx(1))];
        beam_model.Aero.Trim.Link.Slave  = [beam_model.Aero.Trim.Link.Slave,CsNames(RepeatedCSIdx(j))];
        beam_model.Aero.Trim.Link.Coeff  = [beam_model.Aero.Trim.Link.Coeff,1];%CsLink(RepeatedCSIdx(j))
    end
end

% Update some of the counts
beam_model.Info.nlink = numel(beam_model.Aero.Trim.Link.ID );
beam_model.Info.ncaero= length(beam_model.Aero.ID);

% Dummy values
beam_model.Aero.state.AS    = 100;
beam_model.Aero.state.alpha = 0;
beam_model.Aero.state.betha = 0;
beam_model.Aero.state.P     = 0;
beam_model.Aero.state.Q     = 0;
beam_model.Aero.state.R     = 0;
beam_model.Aero.state.SIMXZ = 0;

beam_model.Aero.ref.b_ref = 2*max(beam_model.Node.Coord(:,2));
beam_model.Aero.ref.AR    = [];
beam_model.Aero.ref.S_ref = [];
beam_model.Aero.ref.C_mgc = Fame.Caero(1).c;

outl = vertcat(Fame.Caero.outline);
np   = numel(Fame.Caero);

beam_model.Aero.lattice.XYZ = zeros(np,4,3);

beam_model.Aero.lattice.XYZ(:,:,1) = reshape(outl(:,1),4,np,1)';
beam_model.Aero.lattice.XYZ(:,:,2) = reshape(outl(:,2),4,np,1)';
beam_model.Aero.lattice.XYZ(:,:,3) = reshape(outl(:,3),4,np,1)';

beam_model.Aero.body.ID = [];

%Deal the lumped masses (Conm2)

%Ignore the fuel masses for the time being
Massconfig = findobj(Fame.Conm2, 'type','Fuel','-xor');
%Fuel   = findobj(Fame.Conm2,'type','RightWingfuel_6596');
% Massconfig = cat(1, noFuel, Fuel);

beam_model.Conm2.ID     = [Massconfig.id]';
beam_model.Conm2.Node   = [Massconfig.grid]';
beam_model.Conm2.CID    = [Massconfig.cid]';
for i = 1:numel(beam_model.Conm2.ID)
    beam_model.Conm2.M(:,:,i) = [Massconfig(i).m*eye(3),zeros(3);zeros(3),Massconfig(i).i];
end
beam_model.Conm2.Offset = vertcat(Massconfig.offset);
beam_model.Info.ncom2   = length(beam_model.Conm2.ID);


% Deal the material properties
beam_model.Mat.ID       = [Fame.Mat.id]';
beam_model.Mat.E        = [Fame.Mat.e]';
beam_model.Mat.G        = [Fame.Mat.g]';
beam_model.Mat.nu       = [Fame.Mat.nu]';
beam_model.Mat.Rho      = [Fame.Mat.rho]';
beam_model.Info.nmat    = length(beam_model.Mat.ID);

% Deal the part ids
for i = 1:numel(Fame.PartId)
    beam_model.PartId(i).Id     = Fame.PartId(i).id;
    beam_model.PartId(i).Part   = Fame.PartId(i).part;
    beam_model.PartId(i).Type   = Fame.PartId(i).type;
    beam_model.PartId(i).data   = Fame.PartId(i).data;
end

PartNames = {beam_model.PartId.Part};

WingIdx         = ~cellfun(@isempty,regexpi(PartNames,'wing'));
WingCaeroPartIdx    = and(WingIdx,ismember({beam_model.PartId.Type},'CAERO'));

part            = beam_model.PartId(WingCaeroPartIdx);
WingCaeroIdx = [];
for i = 1:numel(part)
   for j = 1:length(part(i).data)
    WingCaeroIdx = [WingCaeroIdx,find(beam_model.Aero.ID == part(i).data(j))];
   end
end

beam_model.Aero.ref.S_ref  = sum(tarea(beam_model.Aero.lattice.XYZ(WingCaeroIdx,:,:)));

% Deal the Rbe2s
beam_model.RBE2.ID   = vertcat([Fame.Rbe2.id]');
beam_model.RBE2.IDM  = vertcat([Fame.Rbe2.gn]');
for i = 1:numel(beam_model.RBE2.IDM)
    beam_model.RBE2.GDL(i).data  = num2str(Fame.Rbe2(i).dof);
    beam_model.RBE2.IDS(i).data  = Fame.Rbe2(i).gm;
end

beam_model.Info.nrbe2 = length(beam_model.RBE2.ID);

% Deal SPCs
constr = num2str(Fame.Spc.dof);
beam_model.SPC.ID    = Fame.Spc.id;
beam_model.SPC.DOF.list = zeros(1, length(constr));
for i = 1:length(constr)
    beam_model.SPC.DOF.list(i) = int32(str2double(constr(i)));
end
beam_model.SPC.Nodes.list = Fame.Spc.grids;
beam_model.Param.SPC = beam_model.SPC.ID;

beam_model.Info.nspc = length(beam_model.SPC.ID);

% Deal RBE0s
beam_model.RBE0.ID      = [Fame.Rbe0.id]';
beam_model.RBE0.Master  = [Fame.Rbe0.grid]';
for i = 1:numel(Fame.Rbe0)
    beam_model.RBE0.Node(i).data = Fame.Rbe0(i).data';
end
beam_model.Info.nrbe0 = length(beam_model.RBE0.Master);

% Deal Sets
beam_model.Aero.Set.ID = [Fame.Sets.id]';
for i = 1:numel(Fame.Sets)
    beam_model.Aero.Set.Node(i).data = Fame.Sets(i).data';
end
beam_model.Info.nset = length(beam_model.Aero.Set.ID);

% Deal Splines
beam_model.Aero.Interp.ID    = [Fame.Spline.id]';
beam_model.Aero.Interp.Index = [[Fame.Spline.p1]',[Fame.Spline.p2]'];
beam_model.Aero.Interp.Set   = [Fame.Spline.set]';
beam_model.Aero.Interp.Param = [ones(size(beam_model.Aero.Interp.Set)),[Fame.Spline.w]',zeros(size(beam_model.Aero.Interp.Set)),[Fame.Spline.rmx]',[Fame.Spline.cond]'];
beam_model.Aero.Interp.Type  = 2*ones(size(beam_model.Aero.Interp.Set));
beam_model.Aero.Interp.Patch = [Fame.Spline.id]';
beam_model.Info.ninterp      = length(beam_model.Aero.Interp.Set);
beam_model.Info.spline_type  = 2;
beam_model.MeshType = 3;
beam_model.Info.amesh_av_vlm = 1;
beam_model.Param.GRDPNT = 0;

beam_model.WingNodes  = find(and(ismember({Fame.Grid.part},'StbdWing'),ismember({Fame.Grid.type},'Beam')));

end