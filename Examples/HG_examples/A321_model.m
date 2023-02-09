
%% sized parameters

% spar thickness = x(1:25)
% skin thickness = x(26:50)
% stringer thickness = x(51:75)

x=[0.0191135910000000,0.0193953730000000,0.0201137600000000,0.0213823650000000,0.0229999860000000,0.0255251240000000,0.0285624380000000,0.0372113470000000,0.0373207110000000,0.0372799900000000,0.0368848460000000,0.0378039760000000,0.0379687300000000,0.0379303110000000,0.0378464170000000,0.0368959070000000,0.0365928770000000,0.0358741000000000,0.0351157990000000,0.0325204110000000,0.0321366360000000,0.0317223020000000,0.0280889430000000,0.0226448800000000,0.0118821390000000,0.00707914700000000,0.00718361400000000,0.00730636500000000,0.00745189900000000,0.00763438400000000,0.00783787200000000,0.00810313600000000,0.00823100800000000,0.00803824900000000,0.00781619300000000,0.00760075000000000,0.00736034800000000,0.00709415700000000,0.00680991600000000,0.00650766900000000,0.00618441800000000,0.00597655500000000,0.00574729800000000,0.00545697300000000,0.00516124700000000,0.00473319700000000,0.00418339100000000,0.00343190400000000,0.00256919700000000,0.000896999000000000,6.20000000000000e-05,6.38000000000000e-05,6.58000000000000e-05,6.82000000000000e-05,7.14000000000000e-05,7.49000000000000e-05,7.99000000000000e-05,8.21000000000000e-05,7.78000000000000e-05,7.33000000000000e-05,6.92000000000000e-05,6.48000000000000e-05,6.03000000000000e-05,5.55000000000000e-05,5.06000000000000e-05,4.55000000000000e-05,4.24000000000000e-05,3.92000000000000e-05,3.52000000000000e-05,3.13000000000000e-05,2.63000000000000e-05,2.05000000000000e-05,1.40000000000000e-05,8.02000000000000e-06,9.15000000000000e-07];

%% Run folder


if ~isfolder(fullfile(pwd,'bin'))
   mkdir(fullfile(pwd,'bin'))
else
   delete(fullfile('bin','*'))
end
run_folder = fullfile(pwd,'bin');
%% Generate wing parameters

% Aspect_ratio is the only input you need to define the wing geometry
% with increasing Aspect_ratio, the wing geometry will be streched in a
% way that the total wing area, leading edge sweep angle, chord ratio
% to be constant and the wing span, trailing edge sweeps are changing
% based on the selected AR. Engine at fixed position.


Aspect_ratio=10.172; % Aspect ratio = 10.172 for A321 model

Total_area=126;         % include two wing surface areas + floor size on the fuselage
Fuselage_width=4;       % dimeter of the fuselage

Wing_span = sqrt(Aspect_ratio*Total_area);
BeamLoc = 0.4;          % choose a spar location: 0 --> 1
Semi_span=(Wing_span-Fuselage_width)/2; % length of one wing: 16m for A321 model

Root_chord =  Total_area/(1.064*Semi_span + 4);
LE_sweep=27;            % deg

Wing_area = (Total_area - Fuselage_width*Root_chord)/2;

Mid_chord=0.63685*Root_chord;
Tip_chord=0.2248*Root_chord;

X0=Root_chord; 
X1=0.27*Semi_span*tan(27*pi/180) + Mid_chord;
X2=Semi_span*tan(27*pi/180) + Tip_chord;

tan_TE_sweep1=(X1-X0)/(0.27*Semi_span);
tan_TE_sweep2=(X2-X1)/(0.73*Semi_span);

TE_sweep1=atan(tan_TE_sweep1)*180/pi; % deg
TE_sweep2=atan(tan_TE_sweep2)*180/pi; % deg


Taper_ratio=Tip_chord/Root_chord;

Mean_cord_coefficient=(2/3)*(1+Taper_ratio+Taper_ratio^2)/(1+Taper_ratio);

% wing twist parameters
twist_eta = [0,0.75,1];
twist = deg2rad([0,0,0]);


%% Mass configurations

Payload_max=25000; % kg
Payload_percentage=0.5;

Fuel_fraction=0.5; % percentage of fuel in the tank
Fuel_capacity=32940; % L
Fuel_density=840; % g/L
Fuel_mass=Fuel_capacity*Fuel_fraction*Fuel_density/1000; %kg

MTOW=93500; % maximum take off mass
OWE=48500;  % Operating empty mass
MWE=44057;  % Manufacture's empty mass
MZF=73000;  % Maxumum zero fuel mass


Engine_mass=7362/2; % kg
Pylon=1239/2; % kg
Horizontal_tail=682; % kg
Vertical_tail=522; % kg

Secondary_mass=900;% kg

Fuselage_empty_mass=33000; % kg
Fuselage_total_mass=Fuselage_empty_mass+Payload_max*Payload_percentage;


%% A321 Fuselage length, wing positions, engine position

Fuselage_length=45;
Wing_position=20;
Horizontal_tail_position=42;
Vertical_tail_position=41;
Engine_position=4.29;

%% Model generation based on parameters selected above

Connector_right = awi.model.LiftingSurface;
Connector_right.Name = 'Connector_Right';
Connector_right.Origin=[Wing_position,0,0]; %15

%Use the Leading/Trailing edge sweep to define the planform
Connector_right.ActiveSet = 'sSet';

%Tail wing dimensions
Connector_right.SpanVector  = 'Y';
Connector_right.Span        = 2;
Connector_right.LESweep     = [0, 0];
Connector_right.LESweep_eta = [0, 1];
Connector_right.TESweep     = [0,  0];
Connector_right.TESweep_eta = [0,  1];
Connector_right.RootChord   = Root_chord;

%Make sure the beam is at the midchord
all_eta           = Connector_right.Eta_;
Connector_right.BeamLoc     = repmat(BeamLoc, size(all_eta));
Connector_right.BeamLoc_eta = all_eta;


% Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_root_right = awi.model.Spar;
FrontSpar_root_right.XLoc = [0.15, 0.15];
FrontSpar_root_right.Eta  = [0   , 1];
RearSpar_root_right = awi.model.Spar;
RearSpar_root_right.XLoc = [0.65, 0.65];
RearSpar_root_right.Eta  = [0   , 1];
Connector_right.add([FrontSpar_root_right, RearSpar_root_right]);

%Define internal layout
Connector_right.RibPitch      = 0.65;
Connector_right.StringerPitch = 0.15;

%Make the connector material
E0  = 15.4e10; %[N/m^2], typical YM of IM7 composite
nu0 = 0.333;
rho0=1550; 
Mat_conn = awi.model.Material;
Mat_conn.E  = E0;
Mat_conn.Nu = nu0;
Mat_conn.G  = E0 / (2 * (1 + nu0));
Mat_conn.Rho=rho0;

% material properties
Connector_right.Material_eta = [0, 1];
Connector_right.Material     = [Mat_conn, Mat_conn];

% Define box beam corss section
Connector_box_right=awi.model.BoxBeam;
Connector_box_right.BoxType='SymmetricBox';
Connector_box_right.Height=1;
Connector_box_right.Width=3;
Connector_box_right.CoverThickness=0.08;
Connector_box_right.SparThickness=0.08;
getGeometricProps(Connector_box_right)

Connector_right.BoxBeam = Connector_box_right;
Connector_right.A   = Connector_box_right.Abb;
Connector_right.I11 = Connector_box_right.Ixx;
Connector_right.I22 = Connector_box_right.Izz;
Connector_right.J   = Connector_box_right.Jbb;


for i=1:1:3
    handle_connectorR=strcat('PM_tail_R','i');
    handle_connectorR=awi.model.PointMass;
    handle_connectorR.SOffset=-0.1+i*0.2;
    handle_connectorR.Mass=1;
    handle_connectorR.MassGroup='Group3';
    Connector_right.add(handle_connectorR);

end

% Aeropanel definition
Connector_right.AeroPanelLength=0.5;

build(Connector_right);


%% Left root connector

Connector_left = awi.model.LiftingSurface;
Connector_left.Name = 'Connector_Left';
Connector_left.Origin=[Wing_position,0,0]; %15

%Use the Leading/Trailing edge sweep to define the planform
Connector_left.ActiveSet = 'sSet';

%Tail wing dimensions
Connector_left.SpanVector  = 'Y';
Connector_left.Span        = -2;
Connector_left.LESweep     = [0, 0];
Connector_left.LESweep_eta = [0, 1];
Connector_left.TESweep     = [0,  0];
Connector_left.TESweep_eta = [0,  1];
Connector_left.RootChord   = Root_chord;

%Make sure the beam is at the midchord
all_eta           = Connector_left.Eta_;
Connector_left.BeamLoc     = repmat(BeamLoc, size(all_eta));
Connector_left.BeamLoc_eta = all_eta;
%     Connector_right.XOffset=35;
%     Tailwing_right.YOffset=1;

% Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_root_left = awi.model.Spar;
FrontSpar_root_left.XLoc = [0.15, 0.15];
FrontSpar_root_left.Eta  = [0   , 1];
RearSpar_root_left = awi.model.Spar;
RearSpar_root_left.XLoc = [0.65, 0.65];
RearSpar_root_left.Eta  = [0   , 1];
Connector_left.add([FrontSpar_root_left, RearSpar_root_left]);

%Define internal layout
Connector_left.RibPitch      = 0.65;
Connector_left.StringerPitch = 0.15;

% material properties
Connector_left.Material_eta = [0, 1];
Connector_left.Material     = [Mat_conn, Mat_conn];

% Define box beam corss section
Connector_box_left=awi.model.BoxBeam;
Connector_box_left.BoxType='SymmetricBox';
Connector_box_left.Height=1;
Connector_box_left.Width=3;
Connector_box_left.CoverThickness=0.08;
Connector_box_left.SparThickness=0.08;
getGeometricProps(Connector_box_left)

Connector_left.BoxBeam = Connector_box_left;
Connector_left.A   = Connector_box_left.Abb;
Connector_left.I11 = Connector_box_left.Ixx;
Connector_left.I22 = Connector_box_left.Izz;
Connector_left.J   = Connector_box_left.Jbb;


for i=1:1:3
    handle_connectorL=strcat('PM_tail_R','i');
    handle_connectorL=awi.model.PointMass;
    handle_connectorL.SOffset=-0.1+i*0.2;
    handle_connectorL.Mass=1;
    handle_connectorL.MassGroup='Group3';
    Connector_left.add(handle_connectorL);

end

% Aeropanel definition
Connector_left.AeroPanelLength=0.5;

build(Connector_left);


%% Wingbox 1 - right and control surf.

Wingbox_right = awi.model.LiftingSurface;
Wingbox_right.Name = 'A320Wing_right';
Wingbox_right.Origin=[Wing_position,2,0]; %15
%Use the Leading/Trailing edge sweep to define the planform
Wingbox_right.ActiveSet = 'sSet';

% Num of element 
Wingbox_right.NumBeamElem = 23;

%Wing dimensions
Wingbox_right.SpanVector  = 'Y';
Wingbox_right.Span        = Semi_span;   %34.1/2;
Wingbox_right.LESweep     = [LE_sweep, LE_sweep];
Wingbox_right.LESweep_eta = [0, 1];
Wingbox_right.TESweep     = [TE_sweep1, TE_sweep2, TE_sweep2];
Wingbox_right.TESweep_eta = [0, 0.27, 1];
Wingbox_right.RootChord   = Root_chord;


%Dihedral 
Wingbox_right.Dihedral=[5,5];
Wingbox_right.Dihedral_eta=[0,1];


%Make sure the beam is at the midchord
all_eta           = Wingbox_right.Eta_;
Wingbox_right.BeamLoc     = repmat(BeamLoc, size(all_eta));

Wingbox_right.BeamLoc_eta = all_eta;

%Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_right = awi.model.Spar;
FrontSpar_right.XLoc = [0.15, 0.15];
FrontSpar_right.Eta  = [0   , 1];
RearSpar_right = awi.model.Spar;
RearSpar_right.XLoc = [0.65, 0.65];
RearSpar_right.Eta  = [0   , 1];

Wingbox_right.add([FrontSpar_right, RearSpar_right]);

%Define internal layout
Wingbox_right.RibPitch      = 0.65;
Wingbox_right.StringerPitch = 0.15;

%Make the material
E_wing  = 70e9; %[N/m^2], typical YM of aluminium
nu_wing = 0.333;
rho_wing=2810;
Mat_wing = awi.model.Material;
Mat_wing.E  = E_wing;
Mat_wing.Nu = nu_wing;
Mat_wing.G  = E_wing / (2 * (1 + nu_wing));
Mat_wing.Rho=rho_wing;

Wingbox_right.Material_eta = [0, 1];
Wingbox_right.Material     = [Mat_wing, Mat_wing];

build(Wingbox_right)
Wingbox_right.AeroPanelLength = [];


%% Create discretised boxbeam with varied cross section prperties along the span 

NumSec=Wingbox_right.NumBeamElem+2;

%%sizing variables ---------------------------------------------

thickness1=x(1:NumSec);
thickness2=x(NumSec+1:NumSec*2);
Astrg=x(NumSec*2+1:NumSec*3);

d_strg=sqrt(Astrg/0.36);
t_strg=0.12*d_strg;
% -------------------------------------------------------------

% etaS=linspace(0,Wingbox_right.Span,NumSec);

% set width and height array 
YData=Wingbox_right.YData;
SparWidth=Wingbox_right.Chord*0.5;

RootH=Wingbox_right.Chord(1)*0.15; % root thickness/chord = 0.15
MidH=Wingbox_right.Chord(2)*0.12;  % middle thickness/chord = 0.12
TipH=Wingbox_right.Chord(end)*0.11;% tip thickness/chord = 0.11


% set up eta values
elnum=Wingbox_right.NumBeamElem + 1; % total number of beam elements along the wing
Num_seg1=ceil(elnum*0.27); % number of elements in the inboard section
Num_seg2=elnum - Num_seg1; % number of elements in the outboard section

Num_sec1=Num_seg1+1;    % number of sections in the inboard section
Num_sec2=Num_seg2+1;    % number of sections in the outboard section

eta1_=linspace(0,0.27, Num_sec1);
eta2_=linspace(0.27,1,Num_sec2);
etaS=[eta1_(1:end-1),eta2_(1:end)];

RData=Wingbox_right.RData;
eta_R=RData/RData(end);
eta_Y=YData/YData(end);
etaRS=interp1(eta_Y,eta_R,etaS);

Bwidth=interp1(RData/RData(end),SparWidth,etaRS);
Bheight=interp1(RData/RData(end),0.79*[RootH,MidH,TipH],etaRS);

% stringer pitch 
strg_n=0.24;

%intialise data array
A_val=zeros(1,NumSec);
Ixx_val=zeros(1,NumSec);
Izz_val=zeros(1,NumSec);
J_val=zeros(1,NumSec);

%NSM - no use !
NSM_val=zeros(1,NumSec);
NSI_val=zeros(1,NumSec);

%offset from shear center - no use! 
SCy_val=zeros(1,NumSec);
SCz_val=zeros(1,NumSec);
NAy_val=zeros(1,NumSec);
NAz_val=zeros(1,NumSec);
CMy_val=zeros(1,NumSec);
CMz_val=zeros(1,NumSec);

%offset from CoG - no use !
xOff=interp1(YData/YData(end),Wingbox_right.XData,etaS);
xOff_1=xOff(1:end-1);
xOff_2=xOff(2:end);
xOff_val=[0,xOff_2-xOff_1];


for ii=1:NumSec

    boxname=strcat('Box',string(ii));
    boxname=awi.model.BoxBeam;
    boxname.BoxType='SymmetricBox';
    boxname.Height=Bheight(ii);
    boxname.Width=Bwidth(ii);
    boxname.CoverThickness=thickness2(ii);
    boxname.SparThickness=thickness1(ii);

    NumStrg=floor(Bwidth(ii)/strg_n);

    ts=t_strg(ii);
    ds=d_strg(ii);
    hs=Bheight(ii);
    ws=Bwidth(ii);

    Istrg_xx=((ts*ds^3/12)+(ds*ts^3/12 + ds*ts*(ds/2)^2)*2 + 3*ds*ts*(hs/2-ds/2)^2)*NumStrg*2;
    Istrg_zz_=(ds*ts^3/12 + (ts*ds^3/12 + ts*ds*(ds/2)^2)*2);

    if mod(NumStrg,2)==0
        offset=0.12:strg_n:ws/2;
        Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

    elseif mod(NumStrg,2)==1
        offset=0:strg_n:ws/2;
        Istrg_zz=(NumStrg*Istrg_zz_+3*ds*ts*(sum(offset.^2)*2))*2;

    end

    getGeometricProps(boxname)
    A_val(ii)=boxname.Abb+0;
    Ixx_val(ii)=boxname.Ixx+Istrg_xx;
    Izz_val(ii)=boxname.Izz+Istrg_zz;
    J_val(ii)=boxname.Jbb;

    % NSM - didn't use for now!

    NSM_val(ii)=boxname.NSM;
    NSI_val(ii)=boxname.NSI;

    % offset - didn't use for now!
    SCy_val(ii)=boxname.xSC;
    SCz_val(ii)=boxname.zSC;
    NAy_val(ii)=boxname.xNA;
    NAz_val(ii)=boxname.zNA;
    CMy_val(ii)=boxname.xCM;
    CMz_val(ii)=boxname.zCM;

end


eta_=etaRS;
Wingbox_right.A   =  A_val;
Wingbox_right.A_eta=eta_;

Wingbox_right.I11 = Izz_val;
Wingbox_right.I11_eta=eta_;

Wingbox_right.I22 = Ixx_val;
Wingbox_right.I22_eta = eta_;

Wingbox_right.J   = J_val;
Wingbox_right.J_eta= eta_;


% Aeropanel definition 
% TODO: need better mesh espatially when control surfaces are included
% in the model

% AeroPanelLength
%     NumAeroPanel
%     Wingbox_right.NumAeroPanel=20;
Wingbox_right.AeroPanelLength=0.4;

Wingbox_right.Twist = deg2rad([0,5,-10]);
Wingbox_right.Twist_eta = [0,0.75,1];

build(Wingbox_right)


%% Mass definition


wingmass_eta=etaRS;
Mwidth=interp1(RData/RData(end),SparWidth,wingmass_eta);

mass_set=(Secondary_mass+Fuel_mass)*(Mwidth)/sum(Mwidth);

%     m=total_mass/19;

for i=1:1:25
    handle=strcat('PM_right','i');
    handle=awi.model.PointMass;
    handle.SOffset=wingmass_eta(i);
%         handle.SOffset=0+i*0.2;
    handle.Mass=mass_set(i);
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
    handle.MassGroup='Group1';
    Wingbox_right.add(handle);

end



%% attachments - engine

Engine=awi.model.BluffBody;
Engine.Name='Engine';

% cylinder body
Engine.Radius=[1.4, 1.4, 1];
Engine.Eta =  [0, 0.6, 1];
Engine.Length = 3.5;    


% Engine location - user defined
Y_Engine=Engine_position;

X_Engine=interp1(Wingbox_right.YData,Wingbox_right.XData,Y_Engine);
Z_Engine=interp1(Wingbox_right.YData,Wingbox_right.ZData,Y_Engine);

Engine.Origin = [X_Engine-Engine.Length+Wing_position, Y_Engine + 2, Z_Engine];


%Make engine material
E1  = 76e9; %[N/m^2],set as a rigid body
nu = 0.333;
Engine_Mat = awi.model.Material;
Engine_Mat.E  = E1;
Engine_Mat.Nu = nu;
Engine_Mat.G  = E1 / (2 * (1 + nu));
Engine_Mat.Rho=1; % using lumped mass instead


% use the strong material
Engine.Material_eta = [0, 1];
Engine.Material     = [Engine_Mat, Engine_Mat];

% Engine stiffness
Engine_radius=1;
Engine_thickness=0.015;
Engine_A=2*pi*Engine_radius*Engine_thickness;

Engine_I11=pi*Engine_radius^3*Engine_thickness;
Engine_I22=pi*Engine_radius^3*Engine_thickness;
Engine_J=2*pi*Engine_radius^3*Engine_thickness;
%     Engine_Inertia = Mat1.Rho*Engine_A*Engine_radius^2;

Engine.A   = Engine_A;
Engine.I11 = Engine_I11;
Engine.I22 = Engine_I22;
Engine.J   = Engine_J;


%Aeropanel althoufh it is useless now
Engine.AeroPanelLength=0.5;

% add engine mass
engine_mass=awi.model.PointMass;   
engine_mass.SOffset=0.1;
engine_mass.Mass=Engine_mass;
Engine.add(engine_mass);

% add pylon
pylon_mass=awi.model.PointMass;   
pylon_mass.SOffset=0.9;
pylon_mass.Mass=Pylon;
Engine.add(pylon_mass);

build(Engine)

Wingbox_right.add(Engine);


%Control surfaces - flaps
flap_R=awi.model.ControlSurface;
flap_R.Eta=[0, 0.24];
flap_R.xLE=[0.8,0.8];
flap_R.xTE=[1,1];
flap_R.Max_def=0.1;
flap_R.Max_rate=0.1;
flap_R.HingeLine='LE';
flap_R.Label='FlapR';
flap_R.FaceColor='m';

%     flap_R.NumAeroPanel=10;
flap_R.AeroPanelLength=0.4;

build(flap_R)
Wingbox_right.add(flap_R);

Wingbox_right.ModelControlSurf = 1;


build(Wingbox_right);



%% Wingbox 2 - left and control surf.

Wingbox_left = awi.model.LiftingSurface;
Wingbox_left.Name = 'A320Wing_left';
Wingbox_left.Origin=[Wing_position,-2,0];
%Use the Leading/Trailing edge sweep to define the planform
Wingbox_left.ActiveSet = 'sSet';

%Wing dimensions
Wingbox_left.SpanVector  = 'Y';
Wingbox_left.Span        = -Semi_span;  
Wingbox_left.LESweep     = [-LE_sweep, -LE_sweep];
Wingbox_left.LESweep_eta = [0, 1];
Wingbox_left.TESweep     = [-TE_sweep1, -TE_sweep2, -TE_sweep2];
Wingbox_left.TESweep_eta = [0, 0.27, 1];
Wingbox_left.RootChord   = Root_chord;   

%Dihedral 
Wingbox_left.Dihedral=[5,5];
Wingbox_left.Dihedral_eta=[0,1];

%Make sure the beam is at the midchord
all_eta           = Wingbox_left.Eta_;
Wingbox_left.BeamLoc     = repmat(BeamLoc, size(all_eta));
%     Wingbox_right.BeamLoc     = [0.34,0.4,0.4];
Wingbox_left.BeamLoc_eta = all_eta;

%Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_left = awi.model.Spar;
FrontSpar_left.XLoc = [0.15, 0.15];
FrontSpar_left.Eta  = [0   , 1];
RearSpar_left = awi.model.Spar;
RearSpar_left.XLoc = [0.65, 0.65];
RearSpar_left.Eta  = [0   , 1];

Wingbox_left.add([FrontSpar_left, RearSpar_left]);

%Define internal layout
Wingbox_left.RibPitch      = 0.65;
Wingbox_left.StringerPitch = 0.15;

Wingbox_left.Material_eta = [0, 1];
Wingbox_left.Material     = [Mat_wing, Mat_wing];

Wingbox_left.Twist = -1*deg2rad([0,5,-10]);
Wingbox_left.Twist_eta = [0,0.75,1];
build(Wingbox_left)

%% Create discretised boxbeam with varied cross section prperties along the span 

eta_=etaRS;
Wingbox_left.A   =  A_val;
Wingbox_left.A_eta=eta_;

Wingbox_left.I11 = Izz_val;
Wingbox_left.I11_eta=eta_;

Wingbox_left.I22 = Ixx_val;
Wingbox_left.I22_eta = eta_;

Wingbox_left.J   = J_val;
Wingbox_left.J_eta= eta_;


% Aeropanel definition

% AeroPanelLength
%     NumAeroPanel
%     Wingbox_left.NumAeroPanel=20;
Wingbox_left.AeroPanelLength=0.4;

build(Wingbox_left)


%% Mass definition


for i=1:1:25
    handle=strcat('PM_left','i');
    handle=awi.model.PointMass;
    handle.SOffset=wingmass_eta(i);
%         handle.SOffset=0+i*0.2;
    handle.Mass=mass_set(i);
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
    handle.MassGroup='Group1';
    Wingbox_left.add(handle);

end




%% attachments 2  - engine_left

Engine2=awi.model.BluffBody;
Engine2.Name='Engine_left';

% cylinder body
Engine2.Radius=[1.4, 1.4, 1];
Engine2.Eta =  [0, 0.6, 1];
Engine2.Length = 3.5;


Engine2.Origin = [X_Engine-Engine.Length+Wing_position, -(Y_Engine + 2), Z_Engine];


% use the strong material
Engine2.Material_eta = [0, 1];
Engine2.Material     = [Engine_Mat, Engine_Mat];
Engine2.A   = Engine_A;
Engine2.I11 = Engine_I11;
Engine2.I22 = Engine_I22;
Engine2.J   = Engine_J;


%Aeropanel althoufh it is useless now
Engine2.AeroPanelLength=0.5;

% add engine mass
engine2_mass=awi.model.PointMass;   
engine2_mass.SOffset=0.1;
engine2_mass.Mass=Engine_mass;
Engine2.add(engine2_mass);

% add pylon
pylon2_mass=awi.model.PointMass;   
pylon2_mass.SOffset=0.9;
pylon2_mass.Mass=Pylon;
Engine2.add(pylon2_mass);

build(Engine2)

Wingbox_left.add(Engine2);

%Control surfaces - flaps
flap_L=awi.model.ControlSurface;
flap_L.Eta=[0, 0.24];
flap_L.xLE=[0.8,0.8];
flap_L.xTE=[1,1];
flap_L.Max_def=0.1;
flap_L.Max_rate=0.1;
flap_L.HingeLine='LE';
flap_L.Label='FlapL';
flap_L.FaceColor='m';

flap_L.AeroPanelLength=0.4;

build(flap_R)
Wingbox_left.add(flap_L);

Wingbox_left.ModelControlSurf = 1;

build(Wingbox_left);


%% Create a BluffBody

Body=awi.model.BluffBody;
Body.Name='Fuselage';
% cylinder body
% Body.Radius=[2,2];
% Body.Eta=[0,1];

% real body
Body.Eta = [0;0.005268;0.010536;0.015805;0.021073;...
    0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
    0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
    0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
    0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';

Body.Radius = [0.01;0.3844030;0.565081;0.707928;0.830682;0.940375;...
    1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
    1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
    2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
    1.75;1.6;1.4;1.2;1.0;0.01]';
Body.Origin = [0, 0, 0];
Body.Length=Fuselage_length;


%Make the material
E1  = 76e9; %[N/m^2],set as a rigid body
nu = 0.333;
Body_Mat = awi.model.Material;
Body_Mat.E  = E1;
Body_Mat.Nu = nu;
Body_Mat.G  = E1 / (2 * (1 + nu));
Body_Mat.Nu = nu;
Body_Mat.Rho = 2800; 

% use the strong material

Body.Material_eta = [0, 1];
Body.Material     = [Body_Mat, Body_Mat];

%define  panel size
% Body.NumAeroPanel=5;
Body.AeroPanelLength=0.5;


Body_radius=2; 
Body_thickness=0.004;
CS_A=2*pi*Body_radius*Body_thickness;

CS_I11=pi*Body_radius^3*Body_thickness;
CS_I22=pi*Body_radius^3*Body_thickness;
CS_J=2*pi*Body_radius^3*Body_thickness;


Body.A   = CS_A;
Body.I11 = CS_I11;
Body.I22 = CS_I22;
Body.J   = CS_J;

Mass_val=Fuselage_total_mass/11;


for i=1:1:11
    handle=strcat('PM_body','i');
    handle=awi.model.PointMass;
    handle.SOffset=-0.1+i*0.1;
    handle.Mass=Mass_val;
%         handle.Inertia11 =  0;
%         handle.Inertia22 =  0;
%         handle.Inertia33 =  0;
%         handle.Inertia12 =  0;
%         handle.Inertia23 =  0;
%         handle.Inertia13 =  0;
    handle.MassGroup='Group2';
    Body.add(handle);

end

build(Body)


%% Generate tailwing Right and control surf.

Tailwing_right = awi.model.LiftingSurface;
Tailwing_right.Name = 'Tail_Wing_Right';

%Use the Leading/Trailing edge sweep to define the planform
Tailwing_right.ActiveSet = 'sSet';

%Tail wing dimensions
Tailwing_right.SpanVector  = 'Y';
Tailwing_right.Span        = 12.45/2;
Tailwing_right.LESweep     = [32, 32];
Tailwing_right.LESweep_eta = [0, 1];
Tailwing_right.TESweep     = [15,  15];
Tailwing_right.TESweep_eta = [0,  1];
Tailwing_right.RootChord   = 3.31;

%Dihedral
Tailwing_right.Dihedral=[5,5];
Tailwing_right.Dihedral_eta=[0,1];

%Make sure the beam is at the midchord
all_eta           = Tailwing_right.Eta_;
Tailwing_right.BeamLoc     = repmat(0.5, size(all_eta));
Tailwing_right.BeamLoc_eta = all_eta;
Tailwing_right.XOffset = Horizontal_tail_position;
%     Tailwing_right.YOffset=1;

% Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_tail_right = awi.model.Spar;
FrontSpar_tail_right.XLoc = [0.15, 0.15];
FrontSpar_tail_right.Eta  = [0   , 1];
RearSpar_tail_right = awi.model.Spar;
RearSpar_tail_right.XLoc = [0.65, 0.65];
RearSpar_tail_right.Eta  = [0   , 1];
Tailwing_right.add([FrontSpar_tail_right, RearSpar_tail_right]);

%Define internal layout
Tailwing_right.RibPitch      = 0.65;
Tailwing_right.StringerPitch = 0.15;

%Make the material
Et  = 76e9; %[N/m^2],set as a rigid body
nut = 0.333;
Tail_Mat = awi.model.Material;
Tail_Mat.E  = Et;
Tail_Mat.Nu = nut;
Tail_Mat.G  = E1 / (2 * (1 + nut));
Tail_Mat.Nu = nu;
Tail_Mat.Rho = 2800;

% material properties
Tailwing_right.Material_eta = [0, 1];
Tailwing_right.Material     = [Tail_Mat, Tail_Mat];

% Define box beam corss section
tailbox_right=awi.model.BoxBeam;
tailbox_right.BoxType='SymmetricBox';
tailbox_right.Height=0.5;
tailbox_right.Width=1;
tailbox_right.CoverThickness=0.006;
tailbox_right.SparThickness=0.006;
getGeometricProps(tailbox_right)
Tailwing_right.BoxBeam = tailbox_right;
Tailwing_right.A   = tailbox_right.Abb;
Tailwing_right.I11 = tailbox_right.Ixx;
Tailwing_right.I22 = tailbox_right.Izz;
Tailwing_right.J   = tailbox_right.Jbb;


% Aeropanel definition
Tailwing_right.AeroPanelLength=0.5;

%Control surfaces - elevators
myelevator_right=awi.model.ControlSurface;
myelevator_right.Eta=[0, 1];
myelevator_right.xLE=[0.6,0.6];
myelevator_right.xTE=[1,1];
myelevator_right.Max_def=0.1;
myelevator_right.Max_rate=0.1;
myelevator_right.HingeLine='LE';
myelevator_right.Label='elevatR';
myelevator_right.FaceColor='m';
myelevator_right.AeroPanelLength=0.5;
build(myelevator_right)
Tailwing_right.add(myelevator_right);

Tailwing_right.ModelControlSurf = 1;


build(Tailwing_right);


  %% Generate tailwing Left and control surf.

Tailwing_left = awi.model.LiftingSurface;
Tailwing_left.Name = 'Tail_Wing_Left';

%Use the Leading/Trailing edge sweep to define the planform
Tailwing_left.ActiveSet = 'sSet';

%Tail wing dimensions
Tailwing_left.SpanVector  = 'Y';
Tailwing_left.Span        = -12.45/2;
Tailwing_left.LESweep     = [-32, -32];
Tailwing_left.LESweep_eta = [0, 1];
Tailwing_left.TESweep     = [-15,  -15];
Tailwing_left.TESweep_eta = [0,  1];
Tailwing_left.RootChord   = 3.31;

%Dihedral
Tailwing_left.Dihedral=[5,5];
Tailwing_left.Dihedral_eta=[0,1];

%Make sure the beam is at the midchord
all_eta           = Tailwing_left.Eta_;
Tailwing_left.BeamLoc     = repmat(0.5, size(all_eta));
Tailwing_left.BeamLoc_eta = all_eta;
Tailwing_left.XOffset=Horizontal_tail_position;
%     Tailwing_left.YOffset=-1;

% Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
FrontSpar_tail_left = awi.model.Spar;
FrontSpar_tail_left.XLoc = [0.15, 0.15];
FrontSpar_tail_left.Eta  = [0   , 1];
RearSpar_tail_left = awi.model.Spar;
RearSpar_tail_left.XLoc = [0.65, 0.65];
RearSpar_tail_left.Eta  = [0   , 1];
Tailwing_left.add([FrontSpar_tail_left, RearSpar_tail_left]);

%Define internal layout
Tailwing_left.RibPitch      = 0.65;
Tailwing_left.StringerPitch = 0.15;

% material properties
Tailwing_left.Material_eta = [0, 1];
Tailwing_left.Material     = [Tail_Mat, Tail_Mat];

% Define box beam corss section
tailbox_left=awi.model.BoxBeam;
tailbox_left.BoxType='SymmetricBox';
tailbox_left.Height=0.5;
tailbox_left.Width=1;
tailbox_left.CoverThickness=0.006;
tailbox_left.SparThickness=0.006;
getGeometricProps(tailbox_left)
Tailwing_left.BoxBeam = tailbox_left;
Tailwing_left.A   = tailbox_left.Abb;
Tailwing_left.I11 = tailbox_left.Ixx;
Tailwing_left.I22 = tailbox_left.Izz;
Tailwing_left.J   = tailbox_left.Jbb;


% Aeropanel definition
Tailwing_left.AeroPanelLength=0.5;

%Control surfaces - elevators
myelevator_left=awi.model.ControlSurface;
myelevator_left.Eta=[0, 1];
myelevator_left.xLE=[0.6,0.6];
myelevator_left.xTE=[1,1];
myelevator_left.Max_def=0.1;
myelevator_left.Max_rate=0.1;
myelevator_left.HingeLine='LE';
myelevator_left.Label='elevatL';
myelevator_left.FaceColor='m';

myelevator_left.AeroPanelLength=0.5;

build(myelevator_left)
Tailwing_left.add(myelevator_left);

Tailwing_left.ModelControlSurf = 1;


build(Tailwing_left);


%% Generate vertical wing and rudder

Verticalwing=awi.model.LiftingSurface;
Verticalwing.Name = 'Vertical_wing';
Verticalwing.ActiveSet = 'pSet';
Verticalwing.Chord     = [3.31, 1.5];
Verticalwing.Chord_eta = [0, 1];
Verticalwing.Span      = 12.45/2;

Verticalwing.SpanVector = 'Z';
Verticalwing.Sweep = [30, 30];
Verticalwing.Dihedral = [0,0];


all_eta           = Verticalwing.Eta_;
Verticalwing.BeamLoc     = repmat(0.5, size(all_eta));
Verticalwing.BeamLoc_eta = all_eta;
Verticalwing.XOffset=Vertical_tail_position;


Verticalwing.Material_eta = [0, 1];
Verticalwing.Material     = [Tail_Mat, Tail_Mat];

% % Aeropanel definition
Verticalwing.NumAeroPanel=8;

% Define box beam corss section
Verticalbox=awi.model.BoxBeam;
Verticalbox.BoxType='SymmetricBox';
Verticalbox.Height=0.5;
Verticalbox.Width=1;
Verticalbox.CoverThickness=0.005;
Verticalbox.SparThickness=0.005;
getGeometricProps(Verticalbox)
Verticalwing.BoxBeam = Verticalbox;
Verticalwing.A   = Verticalbox.Abb;
Verticalwing.I11 = Verticalbox.Ixx;
Verticalwing.I22 = Verticalbox.Izz;
Verticalwing.J   = Verticalbox.Jbb;


build(Verticalwing);


%% Build aircraft model
Aircraft = awi.model.Aircraft;

Aircraft.add(Body);


Body.add(Connector_right)
Body.add(Connector_left)

Connector_right.add(Wingbox_right)
Connector_left.add(Wingbox_left)

Body.add(Tailwing_right)
Body.add(Tailwing_left)
Body.add(Verticalwing)


%The analysis methods require an 'awi.model.Aircraft' object
% This is because some information is only known at the aircraft level,
% e.g. all-up mass, reference span, area, etc.
% Aircraft = awi.model.Aircraft;
% Aircraft.add(LS);

Aircraft.RefArea  = sum([Wingbox_right.SurfaceArea, Wingbox_left.SurfaceArea,...
    Connector_right.SurfaceArea,  Connector_left.SurfaceArea]);


Aircraft.RefSpan  = Wingbox_right.Span*2+Connector_right.Span*2;
Aircraft.RefChord = Wingbox_right.RootChord*Mean_cord_coefficient; %mean aerodynamic chord = 0.697 for A321 wing;
%     Aircraft.RefChord = Aircraft.RefArea/Aircraft.RefSpan; 

build(Aircraft)


%% Generate the FEM 

% Convert to a finite element model and draw it
FEM_full = convertToFE(Aircraft);
FEM_full.updateDmiEntry();
% %Export it to a file
export(FEM_full, run_folder);

%% NASTRAN method - RUN SOL 144

% Generate the loadcase object

TrimLoadcase1 = awi.model.LoadCase;

acMass = 94000;
altitude          = 36000;
mach_number       = 0.78;
aircraft_velocity = mach_number*340;
flap_angle=0;

TrimLoadcase1.Name = 'A321_cruise_g';
TrimLoadcase1.Altitude   = altitude;
TrimLoadcase1.Mach       = mach_number;
TrimLoadcase1.AcVelocity = aircraft_velocity;
TrimLoadcase1.AcMass = acMass;

TrimLoadcase1.LoadFactor = 1;

% CS deflection - flap deflection angle in radian
TrimLoadcase1.CsDeflection=flap_angle*pi/180;

build(TrimLoadcase1)

% Write input and run analysis

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_full;
MassCases=[];

% AELINK Cards
% get control surfaces
controlsurfs = awi.methods.Nastran.getControlSurfs(FEM_full);
elev_surfs = controlsurfs(cellfun(@(x)contains(x,'elev'),controlsurfs));
if length(elev_surfs) == 2
    label_d = elev_surfs{1};
    aelink = mni.printing.cards.AELINK(elev_surfs{1},{{elev_surfs{2},1}});
else
    error('Expecting only two elevator control surfaces')
end

%Write trim file
trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase1,...
    MassCases,run_folder,'DatFilename','A321_cruise_g',...
    'ExtraCards',aelink);

NastranMethods1.runNastran(trimFile1);

%% plot results

model = mni.import_matran(fullfile(run_folder,'A321_cruise_g.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'A321_cruise_g.f06'));
res_disp =  f06.read_disp;
res_aeroP = f06.read_aeroP;
res_aeroF = f06.read_aeroF;

% apply deformation result
[~,i] = ismember(model.GRID.GID,res_disp.GP);
model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];

% apply aero pressure
model.CAERO1.PanelPressure = res_aeroP.Cp;

%apply aero forces
f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
    res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
model.CAERO1.PanelForce = f';

% update the plot to apply deformations and aero pressures + forces
model.update('Scale',1)

%% NASTRAN method - RUN SOL 146 - Gust analysis Case: crusing altitude + 1MC

GustLoadcase1 = awi.model.LoadCase;
GustLoadcase1.Altitude   = 36000;
GustLoadcase1.AcVelocity = 0.78*340;
GustLoadcase1.AcMass = 500;
GustLoadcase1.Mach = 0.78;
GustLoadcase1.GustLength = linspace(18,214,2);

% Gust direction: positive or negative hit
GustLoadcase1.GustDirection=1;

FlightPoint1=awi.model.FlightPoint;
FlightPoint1.Mach=0.78;
FlightPoint1.AcVelocity=FlightPoint1.Mach*340;
FlightPoint1.Altitude = 36000;
getFlightPointData(FlightPoint1,'ISA');

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM_full;
MassCases=[];

gustfile1=NastranMethods1.writeGustFile(Aircraft, GustLoadcase1, MassCases, FlightPoint1, run_folder,'DatFilename','gust_analysis_g_36000ft_pos_1MC');

NastranMethods1.runNastran(gustfile1);



