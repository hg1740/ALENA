function beam_model = InitRead(filename,fid)

if nargin == 0
    filename = '';
    fid = 1;
end

% Initialise card count
beam_model.Info = [];
beam_model.Info.ngrid    = 0;
beam_model.Info.ncord1   = 0;
beam_model.Info.ncord2   = 0;
beam_model.Info.nbar     = 0;
beam_model.Info.nbeam    = 0;
beam_model.Info.npbar    = 0;
beam_model.Info.npbeam   = 0;
beam_model.Info.ncom1    = 0;
beam_model.Info.ncom2    = 0;
beam_model.Info.ncaero   = 0;
beam_model.Info.nset     = 0;
beam_model.Info.nrbe0    = 0;
beam_model.Info.ninterp  = 0;
beam_model.Info.nlink    = 0;
beam_model.Info.nmat     = 0;
beam_model.Info.ninclude = 1;
beam_model.Info.naeros   = 0;
beam_model.Info.cc_nspc  = 0;
beam_model.Info.ntrim    = 0;
beam_model.Info.cc_trim  = 0;
beam_model.Info.ngust    = 0;
beam_model.Info.nrbe2    = 0;
beam_model.Info.ncelas   = 0;
beam_model.Info.nspc     = 0;
beam_model.Info.nspcd    = 0;
beam_model.Info.nbaero   = 0;
beam_model.Info.naesurf  = 0;
beam_model.Info.ndof     = 0;
beam_model.Info.ndof2    = 0;
beam_model.Info.nload    = 0;
beam_model.Info.ndload   = 0;
beam_model.Info.nf       = 0;
beam_model.Info.nm       = 0;
beam_model.Info.nflw     = 0;
beam_model.Info.ncbush   = 0;
beam_model.Info.npbush   = 0;
beam_model.Info.nrjoint  = 0;
beam_model.Info.spline_type = 0;
beam_model.Info.amesh_av_vlm = 0;
beam_model.Info.nsurfdef = 0;
beam_model.Info.nsubcases = 0;
beam_model.Info.ndamp    = 0;
beam_model.Info.ngrav    = 0;
beam_model.Info.nrbar    = 0;
beam_model.Info.nparts   = 0;

%% Initialise card structures

% beam_model.SPC struct
beam_model.SPC = [];
beam_model.SPC.ID = [];
beam_model.SPC.DOF = [];
beam_model.SPC.Nodes = [];

% beam_model.SPCD Struct
beam_model.SPCD = [];
beam_model.SPCD.ID = [];
beam_model.SPCD.Node = [];
beam_model.SPCD.DOF = [];
beam_model.SPCD.VAD = [];

% beam_model.Param struct
beam_model.Param = [];

% DAMP
beam_model.Param.SDAMP = 0;
beam_model.Param.KDAMP = 1;

% OPTIMIZATION
beam_model.Param.OBJ = '';
beam_model.Param.DESVAR = '';
beam_model.Param.DLINK = '';
beam_model.Param.ST_CONSTR = '';
beam_model.Param.DY_CONSTR = '';
beam_model.Param.DIRDER = 1.0e-5;
beam_model.Param.OPTDVAR = 1.0;

% DLM
beam_model.Param.DLM_ORDER = 2;
beam_model.Param.DLM_NP    = 12;
beam_model.Param.DLM_KMAX  = 0;
beam_model.Param.DLM_AR    = 3;

% flutter plot
beam_model.Param.NVSTEP = 50;
beam_model.Param.VMAX = 0;
beam_model.Param.RHO_VG = 0;
beam_model.Param.FMODES = [];
beam_model.Param.UMODES = [];

% divergence in SOL 144
beam_model.Param.DIVERG = 0;

% external derivatives in SOL 144
beam_model.Param.DER_FILE = '';
beam_model.Param.DER_TYPE = 0;
beam_model.Param.AEROINT = [];

% eigenvalues
beam_model.Param.EIG_FILE = '';
beam_model.Param.MSELECT = [];
beam_model.Param.SUPORT = [];
beam_model.Param.MINF = 0;
beam_model.Param.MAXF = 0;
beam_model.Param.NROOTS = 10;
beam_model.Param.MSCALE = 'MASS';
beam_model.Param.MG = 0;
beam_model.Param.MC = 0;
beam_model.Param.SUP_MAMPL = 1.0;

beam_model.Param.FUSE_DP = 0;

beam_model.Param.FILE = filename;
beam_model.Param.INCLUDE = {};
beam_model.Param.INCLUDE{1} = filename; % dummy
beam_model.Param.AUTOPLOT = false;
beam_model.Param.GRDPNT = int32(0);
beam_model.Param.LOAD = int32(0);
beam_model.Param.Gust = [];
beam_model.Param.SURFDEF = [];
beam_model.Param.SPC = int32(0);
beam_model.Param.SOL = int32(0);
beam_model.Param.MSOL = [];

% landing gear
beam_model.Param.LANDG = []; % node IDs

% non linear beam settings
beam_model.Param.EPS = 1.0e-3;
beam_model.Param.NSTEP = 10;
beam_model.Param.NITER = 5;
beam_model.Param.RES_TOL = 1.0e-6;

% aerodynamic relax factor
beam_model.Param.REL_FAC = 0.5;
beam_model.Param.WTMASS = 1.0;
beam_model.Param.Grav = zeros(3,1);
beam_model.Param.G = 9.81;
beam_model.Grav = []; % dummy variable
beam_model.Param.FID = fid;

beam_model.Param.TIMESTEP = [];
beam_model.Param.NUMSTEPS = [];

beam_model.Param.SUBCASE{1}.ID =[];
beam_model.Param.SUBCASE{1}.LOAD =[];
beam_model.Param.SUBCASE{1}.SPC =[];

% responces settings
beam_model.Param.RHOREF = 0.0;
beam_model.Param.MACH   = 0.0;
beam_model.Param.VREF   = 0.0;
beam_model.Param.MODACC = -1;
beam_model.Param.BCOU = 1;
beam_model.Param.ACCELERATION = [];
beam_model.Param.VELOCITY = [];
beam_model.Param.DISP = [];
beam_model.Param.IFORCE = [];      % internal loads on Bar
beam_model.Param.IFORCEBE = [];    % internal loads on Beam
beam_model.Param.AEROFORCE = [];   % Aerodinamic Loads
beam_model.Param.HINGEFORCE = [];  % Hinge Moments

% FORCE struct
beam_model.F = [];
beam_model.F.ID =   [];
beam_model.F.Type = [];
beam_model.F.Node = [];
beam_model.F.Mag = [];
beam_model.F.Orient = [];
beam_model.F.CID = [];
beam_model.F.Offset = [];

% FORCE struct
beam_model.M = [];
beam_model.M.ID =   [];
beam_model.M.Type = [];
beam_model.M.Node = [];
beam_model.M.Mag = [];
beam_model.M.Orient = [];
beam_model.M.CID = [];

% FOLLOWER struct
beam_model.F_FLW = [];
beam_model.F_FLW.ID =   [];
beam_model.F_FLW.Node = [];
beam_model.F_FLW.Mag = [];
beam_model.F_FLW.Orient = [];
beam_model.F_FLW.CID = [];
beam_model.F_FLW.Offset = [];

% dummy beam_model.COORD1 struct
beam_model.COORD1 = [];
beam_model.COORD1.ID = [];
beam_model.COORD1.Nodes = [];
beam_model.COORD1.R = [];
beam_model.COORD1.Origin = [];

% dummy beam_model.COORD2 struct
beam_model.COORD2 = [];
beam_model.COORD2.ID = [];
beam_model.COORD2.RID = [];
beam_model.COORD2.Nodes = [];
beam_model.COORD2.R = [];
beam_model.COORD2.Origin = [];

% beam_model.Coord struct
beam_model.Coord = [];
beam_model.Coord.ID = [];
beam_model.Coord.Origin = [];
beam_model.Coord.R = [];

% beam_model.Node struct
beam_model.Node  = [];
beam_model.Node.ID = [];
beam_model.Node.CS = [];
beam_model.Node.Coord = [];
beam_model.Node.CD = [];
beam_model.Node.Index = [];
beam_model.Node.R = [];
beam_model.Node.DOF = [];
beam_model.Node.Aero = [];

% Beam struct
beam_model.Beam  = [];
beam_model.Beam.ID = [];       % element ID
beam_model.Beam.PID = [];      % Property ID
beam_model.Beam.Conn = [];     % Connectivity
beam_model.Beam.Orient = [];   % Orientation vector
beam_model.Beam.OffsetT = [];  % Offset type
beam_model.Beam.Offset = [];   % offset values
beam_model.Beam.OffsetN = [];  % offset beam_model.Node
beam_model.Beam.Colloc = [];   % collocation point global coordinates
beam_model.Beam.R = [];        % 5 Reference frame matrices
beam_model.Beam.D = [];        % 2 section stiffness matrices
beam_model.Beam.M = [];        % 3 nodes mass matrix

% PBar struct
beam_model.PBeam = [];
beam_model.PBeam.ID = [];
beam_model.PBeam.Mat = [];
beam_model.PBeam.A = [];
beam_model.PBeam.I = [];
beam_model.PBeam.J = [];
beam_model.PBeam.RhoNS = [];
beam_model.PBeam.Kshear = [];
beam_model.PBeam.X_L = [];
beam_model.PBeam.NSI = [];
beam_model.PBeam.NSCG = [];
beam_model.PBeam.NA = [];
beam_model.PBeam.DA = [];
beam_model.PBeam.DI = [];
beam_model.PBeam.DJ = [];
beam_model.PBeam.DRhoNS = [];

% Bar struct
beam_model.Bar  = [];
beam_model.Bar.ID = [];       % element ID
beam_model.Bar.PID = [];      % Property ID
beam_model.Bar.Conn = [];     % Connectivity
beam_model.Bar.Orient = [];   % Orientation vector
beam_model.Bar.OffsetT = [];  % Offset type
beam_model.Bar.Offset = [];   % offset values
beam_model.Bar.Colloc = [];   % collocation point global coordinates
beam_model.Bar.R = [];        % 5 Reference frame matrices
beam_model.Bar.D = [];        % 2 section stiffness matrices
beam_model.Bar.M = [];        % 3 nodes mass matrix
beam_model.Bar.barg0 = [];    % dummy variable

% PBar struct
beam_model.PBar = [];
beam_model.PBar.ID = [];
beam_model.PBar.Mat = [];
beam_model.PBar.A = [];
beam_model.PBar.I = [];
beam_model.PBar.J = [];
beam_model.PBar.RhoNS = [];
beam_model.PBar.Kshear = [];
beam_model.PBar.Str_point = [];
beam_model.PBar.Type = [];
beam_model.PBar.Section = [];
beam_model.PBar.SI = [];

% RBar1 Struct
beam_model.RBar = [];
beam_model.RBar.ID = [];
beam_model.RBar.GA = [];
beam_model.RBar.GB = [];
beam_model.RBar.CB = [];

% MAT1 struct
beam_model.Mat = [];
beam_model.Mat.ID = [];
beam_model.Mat.E = [];
beam_model.Mat.G = [];
beam_model.Mat.nu = [];
beam_model.Mat.Rho = [];
beam_model.Mat.ST = [];
beam_model.Mat.SC = [];
beam_model.Mat.SS = [];

% Conm1 struct
beam_model.Conm1 = [];
beam_model.Conm1.ID = [];
beam_model.Conm1.Node = [];
beam_model.Conm1.CID = [];
beam_model.Conm1.M = [];

% Conm2 struct
beam_model.Conm2 = [];
beam_model.Conm2.ID = [];
beam_model.Conm2.Node = [];
beam_model.Conm2.CID = [];
beam_model.Conm2.M = [];
beam_model.Conm2.Offset = [];

% Aero struct
beam_model.Aero.ID = [];
beam_model.Aero.CP = [];
beam_model.Aero.IS = [];
beam_model.Aero.INT = [];
beam_model.Aero.AESURF = [];
beam_model.Aero.AESURF.ID = [];
beam_model.Aero.AESURF.Name = {};
beam_model.Aero.AESURF.CID = [];
beam_model.Aero.AESURF.AERID = [];

% BAero struct
beam_model.BAero.ID = [];
beam_model.BAero.CP = [];
beam_model.BAero.SET = [];
beam_model.BAero.Interp = [];

%
beam_model.Aero.geo.ny = [];
beam_model.Aero.geo.nx = [];
beam_model.Aero.geo.startx = [];
beam_model.Aero.geo.starty = [];
beam_model.Aero.geo.startz = [];
beam_model.Aero.geo.c = [];
beam_model.Aero.geo.foil = cell(0);
beam_model.Aero.geo.b = [];
beam_model.Aero.geo.T = [];
beam_model.Aero.geo.SW = [];
beam_model.Aero.geo.TW = [];
beam_model.Aero.geo.dihed = [];
beam_model.Aero.geo.meshtype = [];
beam_model.Aero.geo.flapped = [];
beam_model.Aero.geo.nwing = [];
beam_model.Aero.geo.nelem = [];
beam_model.Aero.geo.CG = [];
beam_model.Aero.geo.ref_point = [];
beam_model.Aero.geo.symetric = [];
beam_model.Aero.geo.fc = [];
beam_model.Aero.geo.fnx = [];
beam_model.Aero.geo.fsym = [];
beam_model.Aero.geo.flap_vector = [];
beam_model.Aero.geo.nc = 0;
beam_model.Aero.geo.TWIST = [];

beam_model.Aero.state.AS = 1;
beam_model.Aero.state.alpha = 0;
beam_model.Aero.state.betha =0;
beam_model.Aero.state.P = 0;
beam_model.Aero.state.Q = 0;
beam_model.Aero.state.R = 0;
beam_model.Aero.state.ALT = 0;
beam_model.Aero.state.rho = 0;
beam_model.Aero.state.pgcorr = 0; % set Prandtl-Glauert correction

beam_model.Aero.body.geo.ref_point = [];
beam_model.Aero.body.geo.fs    = [];
beam_model.Aero.body.geo.Rs    = [];
beam_model.Aero.body.geo.L     = [];
beam_model.Aero.body.geo.R     = [];
beam_model.Aero.body.geo.Rx    = [];
beam_model.Aero.body.geo. x    = [];
beam_model.Aero.body.geo.Nelem = [];
beam_model.Aero.body.geo.CAERO_INT = [];

beam_model.Aeros = [];

% addition to Tornado
beam_model.Aero.state.SIMXZ = 0;
beam_model.Aero.state.SIMXY = 0;
beam_model.Aero.state.Mach     = [];
beam_model.Aero.state.Mach_qhh = [];
beam_model.Aero.state.Kfreq = [];
beam_model.Aero.ref.S_ref = 0;
beam_model.Aero.ref.C_mac = 0;
beam_model.Aero.ref.C_mgc = 0;
beam_model.Aero.ref.b_ref = 0;
beam_model.Aero.lattice = [];
beam_model.Aero.lattice_defo = [];
beam_model.Aero.lattice_vlm = [];
beam_model.Aero.lattice_dlm = [];

% interpolation beam_model.Params
beam_model.Aero.Interp = [];
beam_model.Aero.Interp.ID = [];
beam_model.Aero.Interp.Patch = [];
beam_model.Aero.Interp.Index = [];
beam_model.Aero.Interp.Param = [];
beam_model.Aero.Interp.Type = [];

% strctural interpolation set
beam_model.Aero.Set = [];
beam_model.Aero.Set.ID = [];
beam_model.Aero.Set.Node = [];

% TRIM section
beam_model.Aero.Trim = [];
beam_model.Aero.Trim.ID = [];
beam_model.Aero.Trim.Type = [];
beam_model.Aero.Trim.Select = -1;
beam_model.Aero.Trim.Man_index = [];
beam_model.Aero.Trim.CID = [];
beam_model.Aero.Trim.Mach = [];
beam_model.Aero.Trim.ALT = [];
beam_model.Aero.Trim.FM = [];
beam_model.Aero.Trim.CS = [];
beam_model.Aero.Trim.FM.Fixed = [];
beam_model.Aero.Trim.FM.Value = [];
beam_model.Aero.Trim.CS.Fixed = [];
beam_model.Aero.Trim.CS.Value = [];
beam_model.Aero.Trim.Link = [];
beam_model.Aero.Trim.NC = [];
beam_model.Aero.Trim.Symm = [];
beam_model.Aero.Trim.Param = [];
beam_model.Aero.Trim.Value = [];
beam_model.Aero.Trim.Ext = [];
beam_model.Aero.Trim.MINDEX = [];

% Aero node connections
beam_model.RBE0 = [];
beam_model.RBE0.ID = [];
beam_model.RBE0.Master = [];
beam_model.RBE0.Node = [];

% dummy control surface cell variables
beam_model.Control_Name = cell(0);

% beam_model.RJoint
beam_model.RJoint.ID = [];
beam_model.RJoint.Node = [];
beam_model.RJoint.DOF = [];
beam_model.RJoint.jointg0 = [];
beam_model.RJoint.Orient = [];
beam_model.RJoint.R = [];

%beam_model.Celas
beam_model.Celas.ID = [];
beam_model.Celas.Node = [];

% beam_model.Celas.DOF
beam_model.Celas.K = [];

% beam_model.CBush element
beam_model.CBush.ID = [];
beam_model.CBush.PID = [];
beam_model.CBush.Node = [];
beam_model.CBush.bushg0 = [];
beam_model.CBush.Orient = [];
beam_model.CBush.S = [];
beam_model.CBush.OCID = [];
beam_model.CBush.Offset = [];
beam_model.CBush.R = [];

% beam_model.PBush entries
beam_model.PBush.ID = [];
beam_model.PBush.K = [];
beam_model.PBush.B = [];
beam_model.PBush.GE = [];
beam_model.PBush.RCV = [];

% beam_model.RBE2
beam_model.RBE2.ID = [];
beam_model.RBE2.IDM =[];

% beam_model.Gust
beam_model.Gust.ID = [];
beam_model.Gust.DIR = [];
beam_model.Gust.Amp = [];
beam_model.Gust.Tmax = [];
beam_model.Gust.X0 = [];
beam_model.Gust.funs = {}; % gust space profile
beam_model.Gust.fun = {};  % gust time profile

% beam_model.SURFDEF
beam_model.SURFDEF.ID = [];
beam_model.SURFDEF.Label = {};
beam_model.SURFDEF.Amp = [];
beam_model.SURFDEF.Tmax = [];
beam_model.SURFDEF.X0 = [];
beam_model.SURFDEF.fun = {};
beam_model.SURFDEF.Ti = [];
beam_model.SURFDEF.Npiece = [];

% beam_model.Dextload
beam_model.Dextload.ID = [];
beam_model.Dextload.Node = [];
beam_model.Dextload.NDOF = [];
beam_model.Dextload.Amp = [];
beam_model.Dextload.Tmax = [];
beam_model.Dextload.Tmin = [];
beam_model.Dextload.FLW = [];
beam_model.Dextload.MOM = [];
beam_model.Dextload.fun = {};
beam_model.Dextload.Orient = [];
beam_model.Dextload.Offset = [];

% beam_model.Damping
beam_model.Damping.ID = [];
beam_model.Damping.G = [];
beam_model.Damping.ALPHA1 = [];
beam_model.Damping.ALPHA2 = [];
beam_model.Damping.HYBRID = [];
beam_model.Damping.GEFACT = [];
beam_model.Damping.W3 = [];
beam_model.Damping.W4 = [];
beam_model.Damping.WH = [];

% beam_model.GravITY
beam_model.Grav.ID = [];
beam_model.Grav.COORD = [];
beam_model.Grav.Scale = [];
beam_model.Grav.Orient = [];

% beam_model.PartIdS
beam_model.PartId.Id = [];
beam_model.PartId.Part = '';
beam_model.PartId.Type = '';
beam_model.PartId.data = [];

beam_model.Control_Name = cell(0);

% Other Initialisations
beam_model.Param.INCLUDE = {};

beam_model.WB     = [];
end