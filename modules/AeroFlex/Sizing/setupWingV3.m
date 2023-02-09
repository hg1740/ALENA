%% SETUPWING The following function takes the input properties and sets up
% the stick beam discretisation, the optimsation variables, the aerodynamic
% mesh geometry etc.

% Assumptions. Eta values are in the global system!!
% 
% 
%   Author: Dario Calderon 

function  [Aero,BoxGeo,InputGeo] = setupWingV3(BoxType,InputGeo,BoxGeo)

% Recall some of the global planform properties
SurfaceArea = InputGeo.S_ref;
AR          = InputGeo.AR;
Span        = InputGeo.b_ref;
Rootchord   = InputGeo.c_ref;

% Recall the discretisation of the aero/struc mesh
NumCpanels  = InputGeo.NumCPanels;
NWingBeams  = InputGeo.NBeams;

% Store some of the reference lengths
if isempty(AR) % Fixed Dimensions
    Aero.ref.b_ref = Span;
    Aero.ref.C_mgc = Rootchord;
    Aero.ref.S_ref = SurfaceArea;
else % Variable span and chord
    Aero.ref.b_ref = Span;
    Aero.ref.C_mgc = Rootchord;
    Aero.ref.AR    = AR;
    Aero.ref.S_ref = SurfaceArea; % This does not take into account sweep !!!!!!!!!
end

% Check the reference lengths
tapersec    = 1;
suma        = zeros(1,numel(InputGeo.taper.eta_beg));
Wingsecdiff = diff([InputGeo.taper.eta_beg;1])*Aero.ref.b_ref;

% Calculate the accumulated taper ratio at each section
for i = 1:numel(InputGeo.taper.eta_beg)
    tapersec = [tapersec;tapersec(end)*InputGeo.taper.value_beg(i)];
    suma(i) = 0.5*(tapersec(i) + tapersec(i+1))*Wingsecdiff(i);
end
Sum_Area = sum(suma);

if isempty(AR) % Fixed Dimensions
    Aero.ref.S_ref = Aero.ref.C_mgc*Sum_Area;
    Aero.ref.AR    = (2*Aero.ref.b_ref)^2/(2*Aero.ref.S_ref);
else % Variable span and chord
    Aero.ref.S_ref = SurfaceArea;
    Aero.ref.b_ref = sqrt(2*Aero.ref.S_ref*Aero.ref.AR)/2;
end

% Recalculate the root chord
Aero.ref.C_mgc = Aero.ref.S_ref/Sum_Area;

% Determine the chord values here
InputGeo.chord.eta   = [InputGeo.taper.eta_beg;1];
InputGeo.chord.value = tapersec*Aero.ref.C_mgc;

% Find the unique eta values for (t_c/frontspar/rearspar)
t_c_eta   = BoxGeo.height_fraction.eta;
fspar_eta = BoxGeo.fspar.eta;
aspar_eta = BoxGeo.aspar.eta;

eta_box = uniquetol([unique(t_c_eta);unique(fspar_eta);...
    unique(aspar_eta)],1e-5);

% Calculate the midpoint of eta_box
eta_mbox = 0.5*(eta_box(1:end-1)+eta_box(2:end));

% Interpolate some of the box properties onto a common eta position
t_c   = interp1(t_c_eta,BoxGeo.height_fraction.value,eta_box);
fspar_box = interp1(fspar_eta,BoxGeo.fspar.value,eta_box);
aspar_box = interp1(aspar_eta,BoxGeo.aspar.value,eta_box);

% Determine the thickness at the point of spar skin intersection
Ya = 5*t_c.*(0.2969*sqrt(aspar_box) + (-0.1260)*(aspar_box) + ...
    (-0.3516)*(aspar_box).^2 + 0.2843*(aspar_box).^3 + (-0.1015)*(aspar_box).^4);

Yf = 5*t_c.*(0.2969*sqrt(fspar_box) + (-0.1260)*(fspar_box) + ...
    (-0.3516)*(fspar_box).^2 + 0.2843*(fspar_box).^3 + (-0.1015)*(fspar_box).^4);

% Assume here that the box height is the mean of the two (note Ya & Yf are 
% half thickness)
box_height = (Ya+Yf);

% Recover some of the other eta inputs
dih_eta_beg     = InputGeo.dihedral.eta_beg;
dih_value_beg   = (pi/180)*InputGeo.dihedral.value_beg;
dih_eta_end     = InputGeo.dihedral.eta_end;
dih_value_end   = (pi/180)*InputGeo.dihedral.value_end;

qcsw_eta_beg    = InputGeo.QCsweep.eta_beg;
qcsw_value_beg  = (pi/180)*InputGeo.QCsweep.value_beg;
qcsw_eta_end    = InputGeo.QCsweep.eta_end;
qcsw_value_end  = (pi/180)*InputGeo.QCsweep.value_end;

tap_eta_beg     = [InputGeo.taper.eta_beg;1]; % Add one as taper applies to sec.
tap_value_beg   = InputGeo.taper.value_beg;

% Interpolate eta positions (dih/sweep/taper)
eta_plan = uniquetol([unique(dih_eta_beg);unique(dih_eta_end);...
    unique(qcsw_eta_beg); unique(qcsw_eta_end);...
    unique(tap_eta_beg)],1e-5);

% Calculate the midpoint of the et_plan positions
eta_mplan = 0.5*(eta_plan(1:end-1) + eta_plan(2:end));

% Determine how to discretise the beam elements for the wing sections
SecFrac = diff(eta_plan);

% Interpolate the dihedral/sweep to the new eta_plan positions

% Treat each dih section separately
dihedral = [];
dih_beg = [];
dih_end = [];
[~,idxu] = uniquetol(dih_eta_end,1e-5);
for i = 1:numel(idxu)
    % Identify which values to use
    idx_ib = find((eta_plan<= dih_eta_end(idxu(i))));
    idx_ob = find((eta_plan>= dih_eta_beg(idxu(i))));
    idx    = intersect(idx_ib,idx_ob);
    dihedral(i,:) = interp1([dih_eta_beg(idxu(i)),dih_eta_end(idxu(i))],[dih_value_beg(idxu(i)),dih_value_end(idxu(i))],eta_plan(idx)');
    dih_beg = [dih_beg;dihedral(i,1:end-1)'];
    dih_end = [dih_end;dihedral(i,2:end)'];
end

sweep = [];
qcsw_beg = [];
qcsw_end = [];
[~,idxu] = uniquetol(qcsw_eta_end,1e-5);

for i = 1:numel(idxu)
    % Identify which values to use
    idx_ib = find((eta_plan<= qcsw_eta_end(idxu(i))));
    idx_ob = find((eta_plan>= qcsw_eta_beg(idxu(i))));
    idx    = intersect(idx_ib,idx_ob);
    sweep(i,:) = interp1([qcsw_eta_beg(idxu(i)),qcsw_eta_end(idxu(i))],[qcsw_value_beg(idxu(i)),qcsw_value_end(idxu(i))],eta_plan(idx)');
    qcsw_beg = [qcsw_beg;sweep(i,1:end-1)'];
    qcsw_end = [qcsw_end;sweep(i,2:end)'];
end

% Calculate the mean value of the entries 
dih_mean = 0.5*(dih_beg + dih_end);
qsw_mean = 0.5*(qcsw_beg + qcsw_end);

% Temporarily create a new variable that overides a 90 degree panel! FIX
% THIS!
dih_mean_adj = dih_mean;
dih_mean_adj(dih_mean_adj==pi/2) = 0;

% Adjust the section fractions according the length of the sections not the
% y_eta values
SecFrac = SecFrac./cos(dih_mean_adj);
SecFrac = SecFrac/sum(SecFrac);
SecFrac = SecFrac./cos(qsw_mean);
SecFrac = SecFrac/sum(SecFrac);

% Divide the beams up by the fraction size of the sections
beam_div  = SecFrac* NWingBeams;

% Round the values up
beam_dist = round(beam_div);

% If any sections round to 0 add a single beam element to that section
nullentry = beam_dist == 0;
beam_dist(nullentry)  = 1;
beam_div(nullentry)   = 1;

while sum(beam_dist) < NWingBeams
    [~,mini] = min(beam_div);
    beam_div(mini) = beam_div(mini) + 1;
    beam_dist = round(beam_div);
end

while sum(beam_dist) > NWingBeams
    [~,mini] = max(beam_div);
    beam_div(mini) = beam_div(mini) - 1;
    beam_dist = round(beam_div);
end

 beam_dist(beam_dist == 0) = 1;

% Store the distribution of beam elements
BoxGeo.Beam_dist = beam_dist;
BoxGeo.Nsec      = sum(beam_dist);

% Calculate the eta positions of all of the nodes from the beam
% distribution
eta_beam = [];
for i = 1:numel(beam_dist)
    eta_beam = [eta_beam,linspace(eta_plan(i),eta_plan(i+1),beam_dist(i)+1)];
end

% Get rid of duplicate nodes that might be created.
[eta_beam,edgeindex,~] = unique(eta_beam);

% Calculate the mid point of the eta_beam
eta_mbeam = 0.5*(eta_beam(1:end-1)+eta_beam(2:end));

% Store the eta positions of the beam elements
BoxGeo.y_ends = Aero.ref.b_ref*eta_beam';

% Calculate the midpoint location of the beams
BoxGeo.y_mbox = Aero.ref.b_ref*eta_mbeam';

% Calculate the beam lengths
BoxGeo.y_lbox = diff(BoxGeo.y_ends);

% Store the Beam Indices
BoxGeo.SectionNodeIdx = [1;1+cumsum(beam_dist)];

%% Loop through different box types here!

% Interpolate the thickness onto the midpoints of the beam
t_skin = interp1(BoxGeo.skin.eta,BoxGeo.skin.thickness,eta_mbeam');
A_str  = interp1(BoxGeo.stringer.eta,BoxGeo.stringer.area,eta_mbeam');
t_spar = interp1(BoxGeo.spar.eta,BoxGeo.spar.thickness,eta_mbeam');

% Interpolate various properties onto the midpoint of the beam section
chord_beam = interp1(InputGeo.chord.eta,InputGeo.chord.value,eta_mbeam');
chord_ends = interp1(InputGeo.chord.eta,InputGeo.chord.value,eta_beam');

aspar_beam = interp1(aspar_eta,BoxGeo.aspar.value,eta_mbeam');
aspar_ends = interp1(aspar_eta,BoxGeo.aspar.value,eta_beam');

fspar_beam = interp1(fspar_eta,BoxGeo.fspar.value,eta_mbeam');
fspar_ends = interp1(fspar_eta,BoxGeo.fspar.value,eta_beam');

c_box_beam = (aspar_beam-fspar_beam).*chord_beam;
hspar_beam = interp1(eta_box,box_height,eta_mbeam').*chord_beam;

BoxGeo.X0     = [t_skin;t_spar;A_str];
BoxGeo.t_skin = t_skin;
BoxGeo.t_spar = t_spar;
BoxGeo.A_str  = A_str;
BoxGeo.c_box  = c_box_beam;
BoxGeo.h_spar = hspar_beam;
BoxGeo.a_spar = aspar_beam;
BoxGeo.f_spar = fspar_beam;

% Interpolate the sweep and dihedral onto the eta_plan
sweep_mid = [];
dihed_mid = [];
for i = 1:numel(beam_dist)
    sweep_mid   = [sweep_mid;linspace(qcsw_beg(i),qcsw_end(i),beam_dist(i))'];
    dihed_mid   = [dihed_mid;linspace(dih_beg(i),dih_end(i),beam_dist(i))'];
end

BoxGeo.sweep  = sweep_mid;
BoxGeo.dihed  = dihed_mid;

%%
BoxGeo.NS = floor(BoxGeo.c_box/BoxGeo.SP) + 1;
BoxGeo.NS(BoxGeo.NS < 2) = 2;

% Calculate where the leading and trailing edges are, based on the spar
% locations
if ~isempty(InputGeo.Beamsweep)
    beam_xdiff = [0;BoxGeo.y_lbox.*tan(BoxGeo.sweep)];
    beam_ydiff = [0;BoxGeo.y_lbox];
    beam_zdiff = [0;BoxGeo.y_lbox.*tan(BoxGeo.dihed)];
    Vert_index = find(BoxGeo.dihed == pi/2);
    beam_zdiff(Vert_index+1) = BoxGeo.y_lbox(Vert_index);
    beam_x = cumsum(beam_xdiff) + InputGeo.Offset(1);
    beam_y = cumsum(beam_ydiff) + InputGeo.Offset(2);
    beam_z = cumsum(beam_zdiff) + InputGeo.Offset(3);
    Aero.LE_ends = beam_x -(0.5*box_per + BoxGeo.fspar) .* ((1./box_per)*c_box_ends);
    Aero.QC_ends = beam_x -(0.5*(BoxGeo.aspar + BoxGeo.fspar) - 0.25) .* ((1./box_per)*c_box_ends);
    Aero.TE_ends = beam_x + (0.5*box_per + (1-BoxGeo.aspar)) .* ((1./box_per)*c_box_ends);
elseif ~isempty(InputGeo.QCsweep)
    
    SweepLoc = InputGeo.QCsweep.loc;
    
    % Calculate the true length of the beam elements
    xdiff    = [0;BoxGeo.y_lbox.*tan(BoxGeo.sweep)];

    % Calculate the x coordinate of the quarter chord line (remove the root
    % value of the beam axis and offset the whole thing
    Aero_x   = cumsum(xdiff) + InputGeo.Offset(1) - (0.5*(aspar_ends(1) + fspar_ends(1))-SweepLoc)*chord_ends(1);
    
    % Calculate the x-coordinates of the leading ege points
    Aero.LE_ends = Aero_x - SweepLoc*chord_ends;

    % Calculate the x-coordinates of the trailing ege points
    Aero.TE_ends = Aero_x + (1-SweepLoc)*chord_ends;

    % Store the x-coordinates of the quarter chord points
    Aero.QC_ends = Aero.LE_ends + 0.25*chord_ends;

    % Recover the x-coordinates of the beam points
    beam_x = Aero_x + (0.5*(aspar_ends + fspar_ends)-SweepLoc).* chord_ends;

    % Calculate the z difference of the beam points
    beam_zdiff = [0;BoxGeo.y_lbox.*tan(BoxGeo.dihed)];
    
    % Calculate the z difference of the beam points
    beam_ydiff = [0;BoxGeo.y_lbox];
    
    % Calculate the z difference of the beam points
    Vert_index = find(BoxGeo.dihed == pi/2);
    beam_ydiff(Vert_index+1) = 0*BoxGeo.y_lbox(Vert_index);
    beam_zdiff(Vert_index+1) = BoxGeo.y_lbox(Vert_index);
    beam_y = cumsum(beam_ydiff) + InputGeo.Offset(2);
    beam_z = cumsum(beam_zdiff) + InputGeo.Offset(3);
    beam_sweep = atan((beam_x(2:end)-beam_x(1:end-1))./(beam_y(2:end)-beam_y(1:end-1)));
    
elseif ~isempty(InputGeo.LEsweep)
    
    LE_xdiff =  [0;BoxGeo.y_lbox.*tan(BoxGeo.sweep)];
    LE_x = cumsum(LE_xdiff) + InputGeo.Offset(1) - (0.25*(BoxGeo.a_spar + BoxGeo.fspar) + BoxGeo.fspar).*((1./box_per)*c_box_ends(1));
    Aero.LE_ends = LE_x;
    Aero.TE_ends = LE_x + ((1./box_per)*c_box_ends);
    Aero.QC_ends = LE_x + 0.25*((1./box_per)*c_box_ends);
    beam_x = LE_x + (0.25*(BoxGeo.aspar + BoxGeo.fspar)+BoxGeo.fspar).* ((1./box_per)*c_box_ends);
    
    beam_zdiff = [0;BoxGeo.y_lbox.*tan(BoxGeo.dihed)];
    beam_ydiff = [0;BoxGeo.y_lbox];
    Vert_index = find(BoxGeo.dihed == pi/2);
    beam_ydiff(Vert_index+1) = 0*BoxGeo.y_lbox(Vert_index);
    beam_zdiff(Vert_index+1) = BoxGeo.y_lbox(Vert_index);
    beam_y = cumsum(beam_ydiff) + InputGeo.Offset(2);
    beam_z = cumsum(beam_zdiff) + InputGeo.Offset(3);
end

BoxGeo.c_aero = abs(Aero.LE_ends - Aero.TE_ends);

Outline = zeros(numel(beam_dist),4,3);
for iout = 1:numel(beam_dist)
    Outline(iout,1:4,1) = [Aero.LE_ends(BoxGeo.SectionNodeIdx(iout)),Aero.LE_ends(BoxGeo.SectionNodeIdx(iout+1)),Aero.TE_ends(BoxGeo.SectionNodeIdx(iout+1)),Aero.TE_ends(BoxGeo.SectionNodeIdx(iout))];
    Outline(iout,1:4,2) = [beam_y(BoxGeo.SectionNodeIdx(iout)),beam_y(BoxGeo.SectionNodeIdx(iout+1)),beam_y(BoxGeo.SectionNodeIdx(iout+1)),beam_y(BoxGeo.SectionNodeIdx(iout))];
    Outline(iout,1:4,3) = [beam_z(BoxGeo.SectionNodeIdx(iout)),beam_z(BoxGeo.SectionNodeIdx(iout+1)),beam_z(BoxGeo.SectionNodeIdx(iout+1)),beam_z(BoxGeo.SectionNodeIdx(iout))];
end
% Define some Vortex Lattice panel geometries
NumSPanels  = length(Aero_x)-1;

Aero.Outline = Outline;
Aero.geo.ny = ones(NumSPanels,1);
Aero.geo.nx = NumCpanels*ones(NumSPanels,1);
Aero.geo.c  = Aero.TE_ends-Aero.LE_ends; % Temporary
Aero.geo.startx = Aero.LE_ends(1:end-1);
Aero.geo.starty = beam_y(1:end-1);
Aero.geo.startz = beam_z(1:end-1);
Aero.geo.b  = BoxGeo.y_lbox;
Aero.geo.T  = Aero.geo.c(2:end)./Aero.geo.c(1:end-1);
Aero.geo.c  = Aero.TE_ends(1:end-1)-Aero.LE_ends(1:end-1);

% Calculate the 1/4 chord sweep
Aero.geo.SW = ...
    atan((Aero.LE_ends(2:end)+0.25*Aero.geo.T.*Aero.geo.c - ...
    (Aero.LE_ends(1:end-1)+0.25*Aero.geo.c))./(Aero.geo.b));

Aero.geo.dihed = BoxGeo.dihed;

Aero.geo.TW = zeros(NumSPanels,1,2);

% Interpolate the inboard twist onto the beam etas
jigtwist = (pi/180)*interp1(InputGeo.twist.eta,InputGeo.twist.value,eta_beam');

Aero.geo.TW(:,1,1) = jigtwist(1:end-1);
Aero.geo.TW(:,1,2) = jigtwist(2:end);

Aero.geo.TWIST = 0.5*(jigtwist(1:end-1) + jigtwist(2:end));

Aero.CP = zeros(NumSPanels,1);

% Foil geometry
t_c_ends = interp1(t_c_eta,BoxGeo.height_fraction.value,eta_beam');
if isfield(InputGeo,'Airfoils')
    for i = 1:length(InputGeo.Airfoils.y_eta)
        [EP, Zu(i,:), Zl(i,:), ~, ~] = airfoil_mean_line(InputGeo.Airfoils.label(i),0);
    end
    
    y_eta_new = BoxGeo.y_ends/BoxGeo.y_ends(end);
    for i = 1:size(Zu,2)
        Zu_interp(:,i) = interp1(InputGeo.Airfoils.y_eta',Zu(:,i),y_eta_new);
        Zl_interp(:,i) = interp1(InputGeo.Airfoils.y_eta',Zl(:,i),y_eta_new);
    end
    
    % Write the Airfoils into a temporary Folder
    if ~exist('TemporaryAirfoils','dir')
        mkdir('TemporaryAirfoils');
    end
    
    for i = 1:size(Zu_interp,1)
        dlmwrite(['TemporaryAirfoils/Airfoils_' num2str(i) '.dat'],[size(Zu_interp,2),size(Zu_interp,2);EP,Zu_interp(i,:)';EP,Zl_interp(i,:)'],'precision','%.6f','delimiter',' ');
        NewLabel{i} = ['Airfoils_' num2str(i) '.dat'];
    end
    
    % Interpolate the airfoils onto the required points
    for i = 1:length(t_c_ends)-1
        Aero.geo.foil(i,:,1) = NewLabel(i);
        Aero.geo.foil(i,:,2) = NewLabel(i+1);
    end
else
    for i = 1:length(t_c_ends)-1
        Aero.geo.foil(i,:,1) = {sprintf('%04d',floor(100*t_c_ends(i)))};
        Aero.geo.foil(i,:,2) = {sprintf('%04d',floor(100*t_c_ends(i+1)))};
    end
end


% Number of caero elements
Aero.geo.nwing = NumSPanels;
Aero.geo.nelem = ones(NumSPanels,1);

% Control surface definition
Aero.geo.flapped    = zeros(NumSPanels,1);
Aero.geo.fc         = zeros(NumSPanels,1,2);
Aero.geo.fnx        = zeros(NumSPanels,1);
Aero.geo.nc         = 0;
Aero.geo.fsym       = zeros(NumSPanels,1);
Aero.geo.symetric   = zeros(NumSPanels,1)';
Aero.geo.flap_vector = zeros(NumSPanels,1);
Aero.Control.name   = {};
Aero.Trim.Link.ID   = [];
Aero.Trim.Link.Master = {};
Aero.Trim.Link.Slave  = {};
Aero.Trim.Link.Coeff  = [];

if InputGeo.CS.exist == 1
    Control_count = 0;
    for icont = 1:length(InputGeo.CS.inboard)
        CS_in = InputGeo.CS.inboard(icont)*Aero.ref.b_ref;
        CS_out = InputGeo.CS.outboard(icont)*Aero.ref.b_ref;
        if dih_mean == pi/2
            [~,in_index] = min(abs(beam_z - (CS_in + beam_z(1))));
            [~,out_index] = min(abs(beam_z - (CS_out + beam_z(1))));
        else
            [~,in_index] = min(abs(beam_y - (CS_in + beam_y(1))));
            [~,out_index] = min(abs(beam_y - (CS_out + beam_y(1))));
        end
        if in_index == out_index
            error('Either mesh is too coarse or control surface is ill-defined')
        end
        
        in_pan_index = in_index;
        out_pan_index = out_index - 1;
        Aero.geo.flapped(in_pan_index:out_pan_index) = 1;
        Aero.geo.fc(in_pan_index:out_pan_index,1,:) = 1-InputGeo.CS.chord_le(icont);
        
        if NumCpanels == 1
            Aero.geo.nx(in_pan_index:out_pan_index) = 1;
            Aero.geo.fnx(in_pan_index:out_pan_index) = 1;
        else
            Aero.geo.nx(in_pan_index:out_pan_index) = round(InputGeo.CS.chord_le(icont)*NumCpanels);
            Aero.geo.fnx(in_pan_index:out_pan_index) = round((1-InputGeo.CS.chord_le(icont))*NumCpanels);
        end
        
        NumSpanPanels = out_pan_index - in_pan_index + 1;
        
        for i = 1:NumSpanPanels
            Aero.Control.name{Control_count + i} = [InputGeo.CS.name{icont} num2str(i) 'r']; % OK
        end
        
        Aero.Trim.Link.ID = [Aero.Trim.Link.ID, length(Aero.Trim.Link.ID) + (1:NumSpanPanels-1)];
        %AERO.Trim.Link.ID = [AERO.Trim.Link.ID, Control_count + (1:NumSpanPanels-1)]; % NOT OK
        
        for i = 1:NumSpanPanels-1
            Aero.Trim.Link.Slave  = [Aero.Trim.Link.Slave,Aero.Control.name{Control_count + i+1}];
        end
        
        Control_count = Control_count + out_pan_index - in_pan_index + 1;
        Aero.Trim.Link.Master  = [Aero.Trim.Link.Master,repmat({[InputGeo.CS.name{icont} num2str(1) 'r']},1,NumSpanPanels-1)];
        Aero.Trim.Link.Coeff   = [Aero.Trim.Link.Coeff,ones(1,NumSpanPanels-1)];
        Aero.geo.nc = Aero.geo.nc + NumSpanPanels;
    end
end

Aero.Trim.MasterSurf = unique(Aero.Trim.Link.Master);

panelcount = 1;
for i = 1:numel(Aero.geo.ny)
    subpanelcount = Aero.geo.ny(i)*(Aero.geo.nx(i)+Aero.geo.fnx(i));
    Aero.ID(i) = InputGeo.CAEROIdx + panelcount;
    panelcount = panelcount + subpanelcount;
end

Aero.geo.meshtype = ones(NumSPanels,1);

if BoxType == 2
    A = 2*(BoxGeo.t_skin.*BoxGeo.c_box + BoxGeo.t_spar.*(BoxGeo.h_spar-2*BoxGeo.t_spar));
elseif BoxType == 3
    A = 2*(BoxGeo.NS.*BoxGeo.A_str + BoxGeo.t_skin.*BoxGeo.c_box + BoxGeo.t_spar.*(BoxGeo.h_spar-2*BoxGeo.t_spar));
elseif BoxType == 4
    A = 2*(BoxGeo.t_str.*BoxGeo.c_box + BoxGeo.t_skin.*BoxGeo.c_box + BoxGeo.t_spar.*BoxGeo.h_spar);
end

BoxGeo.nodex = beam_x;
BoxGeo.nodey = beam_y;
BoxGeo.nodez = beam_z;

x_lbox = diff(beam_x);
y_lbox = diff(beam_y);
z_lbox = diff(beam_z);

BoxGeo.l_box = sqrt(x_lbox.^2 + y_lbox.^2 + z_lbox.^2);

BoxGeo.M = sum(BoxGeo.Rho*A.*BoxGeo.l_box);

BoxGeo.Type = BoxType;
BoxGeo.Colours = [0 0 0;1 0 0;0 1 0;0 0 1;1 0 1; 0 1 1;1 1 0; 0.5 0.5 0.5;0.25,1,0.25;1,0.75,0.25];
BoxGeo.PlotStress = 0;
BoxGeo.SAF_MARG = 1.5;

InputGeo.S_ref  = Aero.ref.S_ref;
InputGeo.AR     = Aero.ref.AR;
InputGeo.b_ref  = Aero.ref.b_ref;
InputGeo.c_ref  = Aero.ref.C_mgc;

end