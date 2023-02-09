function [beam_model,PartIDs] = reflectwing(beam_model,PartIDs,WingCAERO)

ngrid = beam_model.Info.ngrid;
INFO   = beam_model.Info;
NODE   = beam_model.Node;
BAR    = beam_model.Bar;
CONM2  = beam_model.Conm2;
RBE0   = beam_model.RBE0;
SPC    = beam_model.SPC;
RBE2   = beam_model.RBE2;
CAERO  = beam_model.Aero;
BEAM   = beam_model.Beam;

BeamNodeIDs  = PartIDs.BeamNodes;
LENodeIDs    = PartIDs.LENodes;
TENodeIDs    = PartIDs.TENodes;
NodeIDs      = PartIDs.Nodes;
if isfield(PartIDs,'BeamConm2')
    BeamConmIDs  = PartIDs.BeamConm2;
else
    BeamConmIDs = [];
end

SecConmIDs   = PartIDs.SecConm2IDs;
if isfield(PartIDs,'CBar')
    CBarIDs  = PartIDs.CBar;
else
    CBarIDs  = [];
end
if isfield(PartIDs,'CBeam')
    CBeamIDs = PartIDs.CBeam;
else
    CBeamIDs = [];
end
RBE0IDs      = PartIDs.RBE0;
SetIDs       = PartIDs.SetIDs;
SplineIDs    = PartIDs.SplineIDs;
CAEROIDs     = PartIDs.CAEROIDs;
if isfield(PartIDs,'RBE2IDs')
    RBE2IDs      = PartIDs.RBE2IDs;
else
    RBE2IDs = [];
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Find indices corresponding to the IDs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
beam_node_index = [];
for j = 1:length(BeamNodeIDs)
    index = find(beam_model.Node.ID == BeamNodeIDs(j)) ;
    beam_node_index = [beam_node_index,index];
end

LE_node_index = [];
for j = 1:length(LENodeIDs)
    index = find(beam_model.Node.ID == LENodeIDs(j)) ;
    LE_node_index = [LE_node_index,index];
end

TE_node_index = [];
for j = 1:length(TENodeIDs)
    index = find(beam_model.Node.ID == TENodeIDs(j)) ;
    TE_node_index = [TE_node_index,index];
end

node_index = [];
for j = 1:length(NodeIDs)
    index = find(beam_model.Node.ID == NodeIDs(j)) ;
    node_index = [node_index,index];
end

cbar_index = [];
for j = 1:length(CBarIDs)
    index = find(beam_model.Bar.ID == CBarIDs(j)) ;
    cbar_index = [cbar_index,index];
end

cbeam_index = [];
for j = 1:length(CBeamIDs)
    index = find(beam_model.Beam.ID == CBeamIDs(j)) ;
    cbeam_index = [cbeam_index,index];
end

rbe0_index = [];
for j = 1:length(RBE0IDs)
    index = find(beam_model.RBE0.ID == RBE0IDs(j)) ;
    rbe0_index = [rbe0_index,index];
end

set_index = [];
for j = 1:length(SetIDs)
    index = find(beam_model.Aero.Set.ID == SetIDs(j)) ;
    set_index = [set_index,index];
end

spline_index = [];
for j = 1:length(SplineIDs)
    index = find(beam_model.Aero.Interp.ID == SplineIDs(j)) ;
    spline_index = [spline_index,index];
end

% Structural Masses
conm2st_index = [];
for j = 1:length(BeamConmIDs)
    index = find(beam_model.Conm2.ID == BeamConmIDs(j)) ;
    conm2st_index = [conm2st_index,index];
end

conm2st_index= sort(conm2st_index);

% Secondary Masses
conm2se_index = [];
for j = 1:length(SecConmIDs)
    index = find(beam_model.Conm2.ID == SecConmIDs(j)) ;
    conm2se_index = [conm2se_index,index];
end

conm2se_index= sort(conm2se_index);

aero_index = [];
for j = 1:length(CAEROIDs)
    index = find(beam_model.Aero.ID == CAEROIDs(j)) ;
    aero_index = [aero_index,index];
end

rbe2_index = [];
for j = 1:length(RBE2IDs)
    index = find(beam_model.RBE2.ID == RBE2IDs(j)) ;
    rbe2_index = [rbe2_index,index];
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Define ID offsets
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
NodeOff = 5000;
BarOff  = 5000;
RBE0Off = 5000;
AeroOff = 5000;
SetOff  = 5000;
IntOff  = 5000;
ConmOff = 5000;
RBE2Off = 5000;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect Nodes
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
ncount = 0;
ngrid = INFO.ngrid;
for i = node_index
    ncount = ncount + 1;
    if NODE.ID(i) == NODE.ID(beam_node_index(1)) || NODE.ID(i) == NODE.ID(TE_node_index(1)) || NODE.ID(i) == NODE.ID(LE_node_index(1))
        % Do not reflect root nodes
    else
        count = count+1;
        NODE.ID(ngrid+count)      = NODE.ID(i) + NodeOff;
        NODE.CS(ngrid+count)      = NODE.CS(i);
        NODE.Coord(ngrid+count,:) = NODE.Coord(i,:)*diag([1,-1,1]);
        NODE.CD(ngrid+count,:)    = NODE.CD(i);
        PartIDs.Nodes = [PartIDs.Nodes,NODE.ID(ngrid+count)];
        if find(PartIDs.BeamNodes == NODE.ID(i))
            PartIDs.BeamNodes = [PartIDs.BeamNodes,NODE.ID(ngrid+count)];
        end
        if find(PartIDs.LENodes == NODE.ID(i))
            PartIDs.LENodes = [PartIDs.LENodes,NODE.ID(ngrid+count)];
        end
        if find(PartIDs.TENodes == NODE.ID(i))
            PartIDs.TENodes = [PartIDs.TENodes,NODE.ID(ngrid+count)];
        end
    end
end
INFO.ngrid = ngrid + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect CBARs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count=0;
nbar = INFO.nbar;
for i = cbar_index
    count = count+1;
    BAR.ID(nbar+count)       = BAR.ID(i) + BarOff;
    BAR.PID(nbar+count)      = BAR.PID(i);
    if count == 1
        BAR.Conn(nbar+count,:)   = [find(NODE.ID == NODE.ID(BAR.Conn(i,1))), 0, find(NODE.ID == NODE.ID(BAR.Conn(i,3))+NodeOff)];
    else
        BAR.Conn(nbar+count,:)   = [find(NODE.ID == NODE.ID(BAR.Conn(i,1))+NodeOff), 0, find(NODE.ID == NODE.ID(BAR.Conn(i,3))+NodeOff)];
    end
    BAR.barg0(nbar+count)    = BAR.barg0(i);
    BAR.Orient(nbar+count,:) = [BAR.Orient(i,1),-BAR.Orient(i,2),BAR.Orient(i,3)];
    BAR.OffsetT(nbar+count)  = BAR.OffsetT(i);
    BAR.Offset(nbar+count,:) = BAR.Offset(i,:);
    PartIDs.CBar = [PartIDs.CBar,BAR.ID(nbar+count)];
end
INFO.nbar = INFO.nbar + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect CBEAMs
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count=0;
nbeam = INFO.nbeam;
for i = cbeam_index
    count = count+1;
    BEAM.ID(nbeam+count)      = BEAM.ID(i) + BarOff;
    BEAM.PID(nbeam+count)     = BEAM.PID(i);
    if count == 1
        BEAM.Conn(nbeam+count,:)   = [find(NODE.ID == NODE.ID(BEAM.Conn(i,1))), 0, find(NODE.ID == NODE.ID(BEAM.Conn(i,3))+NodeOff)];
    else
        BEAM.Conn(nbeam+count,:)   = [find(NODE.ID == NODE.ID(BEAM.Conn(i,1))+NodeOff), 0, find(NODE.ID == NODE.ID(BEAM.Conn(i,3))+NodeOff)];
    end
    BEAM.beamg0(nbeam+count)   = BEAM.beamg0(i);
    BEAM.Orient(nbeam+count,:) = [BEAM.Orient(i,1),-BEAM.Orient(i,2),BEAM.Orient(i,3)];
    BEAM.OffsetT(nbeam+count)  = BEAM.OffsetT(i);
    BEAM.Offset(nbeam+count,:) = BEAM.Offset(i,:);
    PartIDs.CBeam = [PartIDs.CBeam,BEAM.ID(nbeam+count)];
end
INFO.nbeam = INFO.nbeam + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect RBE0s
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
ncount = 0;
nrbe0 = INFO.nrbe0;
for i = rbe0_index
    ncount = ncount +1;
    if ncount == 1
        % Do not duplicate the root RBE0
    else
        count = count+1;
        RBE0.ID(nrbe0 + count)        = RBE0.ID(i) + RBE0Off;
        RBE0.Master(nrbe0 + count)    = find(NODE.ID == NODE.ID(RBE0.Master(i)) + NodeOff);
        for j =1:length(NODE.ID(RBE0.Node(i).data))
            RBE0.Node(nrbe0 + count).data(j,1) = find(NODE.ID == NODE.ID(RBE0.Node(i).data(j,1)) + NodeOff);
        end
        PartIDs.RBE0 = [PartIDs.RBE0,RBE0.ID(nrbe0+count)];
    end
    
end
INFO.nrbe0 = INFO.nrbe0 + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect CAERO Cards & AELINKS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
ncaero = INFO.ncaero;
nlink  = INFO.nlink;
nwcaero = length(WingCAERO.geo.ny);
if isfield(beam_model.Aero.lattice_vlm.Control,'Patch')
    npatch = length(beam_model.Aero.lattice_vlm.Control.Patch);
end
pcount = 0;
%pcount = length(CAERO.lattice_vlm.Control.Name);
for i = 1:nwcaero
    count = count + 1;
    CAERO.ID(ncaero+i)              = WingCAERO.ID(i) + AeroOff;
    CAERO.CP(ncaero+i)              = WingCAERO.CP(i);
    CAERO.geo.startx(ncaero+i)      = WingCAERO.geo.startx(i);
    CAERO.geo.starty(ncaero+i)      = -WingCAERO.geo.starty(i);
    CAERO.geo.startz(ncaero+i)      = WingCAERO.geo.startz(i);
    CAERO.geo.c(ncaero+i)           = WingCAERO.geo.c(i);
    CAERO.geo.b(ncaero+i)           = -WingCAERO.geo.b(i);
    CAERO.geo.T(ncaero+i)           = WingCAERO.geo.T(i);
    CAERO.geo.SW(ncaero+i)          = -WingCAERO.geo.SW(i);
    CAERO.geo.foil(ncaero+i,:,1)    = WingCAERO.geo.foil(i,:,1);
    CAERO.geo.foil(ncaero+i,:,2)    = WingCAERO.geo.foil(i,:,2);
    CAERO.geo.dihed(ncaero+i)       = -WingCAERO.geo.dihed(i);
    CAERO.geo.TW(ncaero+i,:,1)      = WingCAERO.geo.TW(i,:,1);
    CAERO.geo.TW(ncaero+i,:,2)      = WingCAERO.geo.TW(i,:,2);
    CAERO.geo.nx(ncaero+i)          = WingCAERO.geo.nx(i);
    CAERO.geo.ny(ncaero+i)          = WingCAERO.geo.ny(i);
    CAERO.geo.flapped(ncaero+i)     = WingCAERO.geo.flapped(i);
    CAERO.geo.nc                    = 2*WingCAERO.geo.nc;
    CAERO.geo.fnx(ncaero+i)         = WingCAERO.geo.fnx(i);
    CAERO.geo.meshtype(ncaero+i)    = WingCAERO.geo.meshtype(i);
    CAERO.geo.nelem(ncaero+i)       = WingCAERO.geo.nelem(i);
    CAERO.geo.fsym(ncaero+i)        = WingCAERO.geo.fsym(i);
    CAERO.geo.flap_vector(ncaero+i) = WingCAERO.geo.flap_vector(i);
    CAERO.geo.fc(ncaero+i,1,1)      = WingCAERO.geo.fc(i,1,1);
    CAERO.geo.fc(ncaero+i,1,2)      = WingCAERO.geo.fc(i,1,2);
    CAERO.geo.symetric(ncaero+i)    = WingCAERO.geo.symetric(i);
    CAERO.geo.TWIST(ncaero+i)       = WingCAERO.geo.TWIST(i);
    if WingCAERO.geo.fc(i)
        pcount = pcount + 1;
        npatch = npatch + 1;
        CAERO.lattice_vlm.Control.Patch(npatch) = ncaero+i;
        CAERO.lattice_vlm.Control.Name{npatch} = ...
            [WingCAERO.Control.name{pcount}(1:end-1) 'l'];
    end
    PartIDs.CAEROIDs = [PartIDs.CAEROIDs,CAERO.ID(ncaero+i)];
end
INFO.ncaero = ncaero + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect AELINKS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

nlink = INFO.nlink;
if pcount>0
    % A lot of fudging is needed here to correctly reflect the links
    [~,Mstr_Idx] = unique(WingCAERO.Trim.Link.Master,'stable');
    if ~isempty(Mstr_Idx)
        NumMaster = length(Mstr_Idx);
        AddSlave = (0:NumMaster-1)';
        NewMstrIdx = Mstr_Idx + AddSlave;
        DiffMstrIdx = diff(NewMstrIdx);
        
        MasterIndex = zeros(pcount,1);
        count = 1;
        for i = 1:length(DiffMstrIdx)
            MasterIndex(count:count + DiffMstrIdx(i)-1,1) = Mstr_Idx(i);
            count = count + DiffMstrIdx(i);
        end
        MasterIndex(count:pcount,1) = Mstr_Idx(end);
        count = 0;
        indxcount = 0;
        indxcmpr  = 1;
        for i = 1:pcount
            if indxcmpr ~= MasterIndex(i)
                indxcount = 0;
                indxcmpr = MasterIndex(i);
            end
            indxcount = indxcount + 1;
            count = count + 1;
            CAERO.Trim.Link.ID(nlink+i)   =  nlink + i;
            CAERO.Trim.Link.Master{nlink+i} =  WingCAERO.Trim.Link.Master{MasterIndex(i)};
            CAERO.Trim.Link.Slave{nlink+i}  =  [WingCAERO.Trim.Link.Slave{MasterIndex(i)}(1:end-2) num2str(indxcount) 'l'];
            CAERO.Trim.Link.Coeff(nlink+i)  = -WingCAERO.Trim.Link.Coeff(1);
        end
        INFO.nlink = nlink + count;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect Control surface labels
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect SET cards
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
nset = INFO.nset;
for i = set_index
    count = count + 1;
    CAERO.Set.ID(count + nset) =  CAERO.Set.ID(i) + SetOff;
    for j = 1: length(CAERO.Set.Node(i).data)
        if NODE.ID(CAERO.Set.Node(i).data(j)) == NODE.ID(beam_node_index(1)) || NODE.ID(CAERO.Set.Node(i).data(j)) == NODE.ID(TE_node_index(1)) || NODE.ID(CAERO.Set.Node(i).data(j)) == NODE.ID(LE_node_index(1))
            CAERO.Set.Node(count + nset).data(j,1) = ...
                find(NODE.ID == NODE.ID(CAERO.Set.Node(i).data(j)));
        else
            CAERO.Set.Node(count + nset).data(j,1) = ...
                find(NODE.ID == NODE.ID(CAERO.Set.Node(i).data(j)) + NodeOff);
        end
    end
    PartIDs.SetIDs = [PartIDs.SetIDs,CAERO.Set.ID(count + nset)];
end
INFO.nset = nset + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect Spline Cards
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
ninterp = INFO.ninterp;
InterpSet = max(CAERO.Interp.Set) + 1;
for i = spline_index
    count = count + 1;
    CAERO.Interp.ID(1,count+ninterp) = CAERO.Interp.ID(1,i) + IntOff;
    CAERO.Interp.Patch(1,count + ninterp) = CAERO.Interp.Patch(1,count-1 + ninterp) + 1;
    CAERO.Interp.Param(count + ninterp,:) = CAERO.Interp.Param(i,:);
    CAERO.Interp.Index(count + ninterp,1) = CAERO.Interp.Index(i,1);
    CAERO.Interp.Index(count + ninterp,2) = CAERO.Interp.Index(i,2);
    CAERO.Interp.Set(count + ninterp) = InterpSet;
    PartIDs.SplineIDs = [PartIDs.SplineIDs,CAERO.Interp.ID(1,count+ninterp)];
end
INFO.ninterp = ninterp + count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect Structural Masses
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
ncom2 = INFO.ncom2;
count=0;
for i = conm2st_index
    count=count+1;
    CONM2.ID(count+ncom2)     = CONM2.ID(i) + ConmOff;
    if CONM2.Node(i) == NODE.ID(node_index(1))
        CONM2.Node(count+ncom2)   = CONM2.Node(i);
    else
        CONM2.Node(count+ncom2)   = CONM2.Node(i) + NodeOff;
    end
    CONM2.CID(count+ncom2)    = CONM2.CID(i);
    CONM2.M(:,:,count+ncom2)  = CONM2.M(:,:,i);
    CONM2.Offset(count+ncom2,:) = CONM2.Offset(i,:)*diag([1,-1,1]);
    PartIDs.BeamConm2 = [PartIDs.BeamConm2,CONM2.ID(count+ncom2)];
end
INFO.ncom2= ncom2+count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect Secondary Masses
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
ncom2 = INFO.ncom2;
count=0;
for i = conm2se_index
    count=count+1;
    CONM2.ID(count+ncom2)     = CONM2.ID(i) + ConmOff;
    if CONM2.Node(i) == NODE.ID(node_index(1))
        CONM2.Node(count+ncom2)   = CONM2.Node(i);
    else
        CONM2.Node(count+ncom2)   = CONM2.Node(i) + NodeOff;
    end
    CONM2.CID(count+ncom2)    = CONM2.CID(i);
    CONM2.M(:,:,count+ncom2)  = CONM2.M(:,:,i);
    CONM2.Offset(count+ncom2,:) = CONM2.Offset(i,:)*diag([1,-1,1]);
    PartIDs.SecConm2IDs = [PartIDs.SecConm2IDs,CONM2.ID(count+ncom2)];
end
INFO.ncom2= ncom2+count;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Reflect RBE2 Cards
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
count = 0;
nrbe2 = INFO.nrbe2;
for i = rbe2_index
    % Only reflect the RBE2 as long as the slave node is not the centreline
    % node
    if RBE2.IDS(i).data ~= NODE.ID(node_index(1))
        count = count + 1;
        %         if RBE2.ID(i) == NODE.ID(node_index(1))
        %             RBE2.ID(count + nrbe2) = RBE2.ID(i);
        %         else
        RBE2.ID(count + nrbe2) = RBE2.ID(i) + RBE2Off;
        %         end
        RBE2.IDM(count + nrbe2)  = RBE2.IDM(i);
        RBE2.GDL(count + nrbe2).data  = RBE2.GDL(i).data;
        RBE2.IDS(count + nrbe2).data  = RBE2.IDS(i).data + NodeOff;
    end
end
INFO.nrbe2 = nrbe2 + count;

beam_model.Info  = INFO;
beam_model.Node  = NODE;
beam_model.Bar   = BAR;
beam_model.RBE0  = RBE0;
beam_model.SPC   = SPC;
beam_model.Aero  = CAERO;
beam_model.Conm2 = CONM2;
beam_model.RBE2  = RBE2;
beam_model.Beam  = BEAM;
end