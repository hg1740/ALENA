function data = Combine_Aero(data,Parts)

% if idx == 0
%     Parts = {'Wing','HTP'};
% else
%     Parts = {'Wing','HTP','VTP'};
% end

index = [];
for i = 1:length(Parts)
    index_n = find(strcmp(fieldnames(data.Aero), Parts{i})==1);
    if ~isempty(index_n)
        index_n = i;
    end
    index = [index,index_n];
end

data.Aero.Control.Name = {};
data.Aero.Trim.Link.ID = [];
data.Aero.Trim.Link.Master = {}; 
data.Aero.Trim.Link.Slave  = {};
data.Aero.Trim.Link.Coeff = [];
data.Aero.ID = [];
data.Aero.CP = [];
data.Aero.geo.ny = [];
data.Aero.geo.nx = [];
data.Aero.geo.c = [];
data.Aero.geo.b = [];
data.Aero.geo.startx = [];
data.Aero.geo.starty = [];
data.Aero.geo.startz = [];
data.Aero.geo.T = [];
data.Aero.geo.SW = [];
data.Aero.geo.dihed = [];
data.Aero.geo.TW = [];
data.Aero.geo.TWIST = [];
data.Aero.geo.foil = {};
data.Aero.geo.fc = [];
data.Aero.geo.nwing = 0;
data.Aero.geo.nelem = [];
data.Aero.geo.flapped = [];
data.Aero.geo.nc = 0;
data.Aero.geo.fnx = [];
data.Aero.geo.fsym = [];
data.Aero.geo.symetric = [];
data.Aero.geo.flap_vector = [];
data.Aero.geo.meshtype = [];
data.Aero.Control.Name = {};
data.Aero.Trim.Link.ID = [];
data.Aero.Trim.Link.Master = {};
data.Aero.Trim.Link.Slave = {};
data.Aero.Trim.Link.Coeff = [];
data.Info.nlink = 0;

for j = 1:length(index)
    
    ID              = [data.Aero.ID,data.Aero.(Parts{index(j)}).ID];
    CP              = [data.Aero.CP;data.Aero.(Parts{index(j)}).CP];
    ny              = [data.Aero.geo.ny;data.Aero.(Parts{index(j)}).geo.ny];
    nx              = [data.Aero.geo.nx;data.Aero.(Parts{index(j)}).geo.nx];
    c               = [data.Aero.geo.c;data.Aero.(Parts{index(j)}).geo.c];
    b               = [data.Aero.geo.b;data.Aero.(Parts{index(j)}).geo.b];
    startx          = [data.Aero.geo.startx;data.Aero.(Parts{index(j)}).geo.startx];
    starty          = [data.Aero.geo.starty;data.Aero.(Parts{index(j)}).geo.starty];
    startz          = [data.Aero.geo.startz;data.Aero.(Parts{index(j)}).geo.startz];
    T               = [data.Aero.geo.T;data.Aero.(Parts{index(j)}).geo.T];
    SW              = [data.Aero.geo.SW;data.Aero.(Parts{index(j)}).geo.SW];
    dihed           = [data.Aero.geo.dihed;data.Aero.(Parts{index(j)}).geo.dihed];
    TW1             = [data.Aero.geo.TW(:,:,1);data.Aero.(Parts{index(j)}).geo.TW(:,:,1)];
    TWIST           = [data.Aero.geo.TWIST;data.Aero.(Parts{index(j)}).geo.TWIST];
    foil1           = [data.Aero.geo.foil(:,:,1);data.Aero.(Parts{index(j)}).geo.foil(:,:,1)];
    fc1             = [data.Aero.geo.fc(:,:,1);data.Aero.(Parts{index(j)}).geo.fc(:,:,1)];
    
    if size(data.Aero.geo.TW,3) == 1
        TW2   = data.Aero.(Parts{index(j)}).geo.TW(:,:,2);
        foil2 = data.Aero.(Parts{index(j)}).geo.foil(:,:,2);
        fc2   = data.Aero.(Parts{index(j)}).geo.fc(:,:,2);
    else
        TW2   = [data.Aero.geo.TW(:,:,2);data.Aero.(Parts{index(j)}).geo.TW(:,:,2)];
        foil2 = [data.Aero.geo.foil(:,:,2);data.Aero.(Parts{index(j)}).geo.foil(:,:,2)];
        fc2   = [data.Aero.geo.fc(:,:,2);data.Aero.(Parts{index(j)}).geo.fc(:,:,2)];
    end

    nwing       = data.Aero.geo.nwing + data.Aero.(Parts{index(j)}).geo.nwing;
    nelem       = [data.Aero.geo.nelem;data.Aero.(Parts{index(j)}).geo.nelem];
    flapped     = [data.Aero.geo.flapped;data.Aero.(Parts{index(j)}).geo.flapped];
    nc          = data.Aero.geo.nc + data.Aero.(Parts{index(j)}).geo.nc;
    
    fnx         = [data.Aero.geo.fnx;data.Aero.(Parts{index(j)}).geo.fnx];
    fsym        = [data.Aero.geo.fsym;data.Aero.(Parts{index(j)}).geo.fsym];
    symetric    = [data.Aero.geo.symetric,data.Aero.(Parts{index(j)}).geo.symetric];
    flap_vector = [data.Aero.geo.flap_vector;data.Aero.(Parts{index(j)}).geo.flap_vector];
    meshtype    = [data.Aero.geo.meshtype;data.Aero.(Parts{index(j)}).geo.meshtype];
    
    ControlName = [data.Aero.Control.Name,data.Aero.(Parts{index(j)}).Control.name];
    
    LinkID      = [data.Aero.Trim.Link.ID,data.Aero.(Parts{index(j)}).Trim.Link.ID + data.Info.nlink];
    Master      = [data.Aero.Trim.Link.Master,data.Aero.(Parts{index(j)}).Trim.Link.Master];
    Slave       = [data.Aero.Trim.Link.Slave,data.Aero.(Parts{index(j)}).Trim.Link.Slave];
    Coeff       = [data.Aero.Trim.Link.Coeff,data.Aero.(Parts{index(j)}).Trim.Link.Coeff];
    nlink       = length(data.Aero.(Parts{index(j)}).Trim.Link.Coeff);
    
    data.Aero.ID              = ID;
    data.Aero.CP              = CP;
    data.Aero.geo.ny          = ny;
    data.Aero.geo.nx          = nx;
    data.Aero.geo.c           = c;
    data.Aero.geo.startx      = startx;
    data.Aero.geo.starty      = starty;
    data.Aero.geo.startz      = startz;
    data.Aero.geo.b           = b;
    data.Aero.geo.T           = T;
    data.Aero.geo.SW          = SW;
    data.Aero.geo.dihed       = dihed;
    data.Aero.geo.TW          = [];
    data.Aero.geo.TW(:,:,1)   = TW1;
    data.Aero.geo.TW(:,:,2)   = TW2;
    data.Aero.geo.TWIST       = TWIST;
    data.Aero.geo.foil        = {};
    data.Aero.geo.foil(:,:,1) = foil1;
    data.Aero.geo.foil(:,:,2) = foil2;
    data.Aero.geo.nwing       = nwing;
    data.Aero.geo.nelem       = nelem;
    data.Aero.geo.flapped     = flapped;
    data.Aero.geo.nc          = nc;
    data.Aero.geo.fc          = [];
    data.Aero.geo.fc(:,:,1)   = fc1;
    data.Aero.geo.fc(:,:,2)   = fc2;
    data.Aero.geo.fnx         = fnx;
    data.Aero.geo.fsym        = fsym;
    data.Aero.geo.symetric    = symetric;
    data.Aero.geo.flap_vector = flap_vector;
    data.Aero.geo.meshtype    = meshtype;
    data.Aero.Control.Name    = ControlName;
    data.Aero.Trim.Link.ID    = LinkID;
    data.Aero.Trim.Link.Master   = Master;
    data.Aero.Trim.Link.Slave    = Slave;
    data.Aero.Trim.Link.Coeff    = Coeff;
    data.Info.nlink = data.Info.nlink + nlink;
end

data.Aero.Trim.MasterSurf = unique(data.Aero.Trim.Link.Master);

WingIdx = ~cellfun(@isempty,regexp(Parts,'Wing'));
WingPart = Parts(WingIdx);

data.Aero.ref.b_ref = data.Aero.(WingPart{1}).ref.b_ref;
data.Aero.ref.AR    = data.Aero.(WingPart{1}).ref.AR;
data.Aero.ref.S_ref = data.Aero.(WingPart{1}).ref.S_ref;
data.Aero.ref.C_mgc = data.Aero.(WingPart{1}).ref.C_mgc;

end