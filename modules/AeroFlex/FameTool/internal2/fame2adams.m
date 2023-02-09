%% FAME2ADAMS Convert FAME model object into MSC.ADAMS format command file
%
% This function creates mass, beams and aerodynamic panels for NEOCASS and
% NASTRAN.
% 
% FAME2ADAMS(FILENAME,FAME2MATOBJECT)
%
%       FILENAME : Hard coded as |Fame2AdamsConversion.cmd| for now. 
%                  TODO : An ADAMS file name option needs to be created    
% FAME2MATOBJECT : Fame model object
%
%
%   See also |writeNeocassFiles| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function [] = fame2adams(fileName, obj )


fid = fopen(fileName,'w');

writeHeader(fid);
writeParts(fid,obj);
writeForces(fid,obj);
writeSimscript(fid);
writeDynamicgraphics(fid);
writeAcc(fid);
writeExpressdef(fid);

fclose(fid);

end



%% Write Header
function [] = writeHeader(fid)

fprintf(fid,'!\n');
fprintf(fid,'!-------------------------- Default Units for Model ---------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults units  &\n');
fprintf(fid,'   length = meter  &\n');
fprintf(fid,'   angle = deg  &\n');
fprintf(fid,'   force = newton  &\n');
fprintf(fid,'   mass = kg  &\n');
fprintf(fid,'   time = sec\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults units  &\n');
fprintf(fid,'   coordinate_system_type = cartesian  &\n');
fprintf(fid,'   orientation_type = body313\n');
fprintf(fid,'!\n');
fprintf(fid,'!------------------------ Default Attributes for Model ------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults attributes  &\n');
fprintf(fid,'   inheritance = bottom_up  &\n');
fprintf(fid,'   icon_visibility = on  &\n');
fprintf(fid,'   grid_visibility = off  &\n');
fprintf(fid,'   size_of_icons = 0.05  &\n');
fprintf(fid,'   spacing_for_grid = 1.0\n');
fprintf(fid,'!\n');
fprintf(fid,'!------------------------------ Adams/View Model ------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'model create  &\n');
fprintf(fid,'   model_name = WING_BEAM\n');
fprintf(fid,'!\n');
fprintf(fid,'view erase\n');
fprintf(fid,'!\n');
fprintf(fid,'!--------------------------------- Materials ----------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'material create  &\n');
fprintf(fid,'   material_name = .WING_BEAM.FAME_Material  &\n');
fprintf(fid,'   adams_id = 1  &\n');
fprintf(fid,'   density = 1.0  &\n');
fprintf(fid,'   youngs_modulus = 1.0  &\n');
fprintf(fid,'   poissons_ratio = 0.33\n');
fprintf(fid,'!\n');

end

%% Write Parts
% Each Conm2 object can be represented as a part with specific mass and
% inertia. The global coordinate system is used. 
%
% Note that the inertia terms are defined differently than NASTRAN, so the
% off-diagonal terms are multiplied by -1.
function [] = writeParts(fid,obj)

% Write Ground Part
fprintf(fid,'!-------------------------------- Rigid Parts ---------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'! Create parts and their dependent markers and graphics\n');
fprintf(fid,'!\n');
fprintf(fid,'!----------------------------------- ground -----------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'! ****** Ground Part ******\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults model  &\n');
fprintf(fid,'   part_name = ground\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults coordinate_system  &\n');
fprintf(fid,'   default_coordinate_system = .WING_BEAM.ground\n');
fprintf(fid,'!\n');
fprintf(fid,'! ****** Markers for current part ******\n');
fprintf(fid,'!\n');
fprintf(fid,'marker create  &\n');
fprintf(fid,'   marker_name = .WING_BEAM.ground.MARKER_REF  &\n');
fprintf(fid,'   adams_id = 100000  &\n');
fprintf(fid,'   location = 0.0, 0.0, 0.0  &\n');
fprintf(fid,'   orientation = 0.0d, 0.0d, 0.0d\n');
fprintf(fid,'!\n');
fprintf(fid,'part create rigid_body mass_properties  &\n');
fprintf(fid,'   part_name = .WING_BEAM.ground  &\n');
fprintf(fid,'   material_type = .WING_BEAM.FAME_Material\n');
fprintf(fid,'!\n');
fprintf(fid,'part attributes  &\n');
fprintf(fid,'   part_name = .WING_BEAM.ground  &\n');
fprintf(fid,'   name_visibility = off\n');
fprintf(fid,'!\n');

% Find structural nodes
idxNode = find([obj.Mdl.Grid.id]' < 130000);
Coord = struct2mat(obj.Mdl.Grid,'coord');
wingSemiSpan = max(Coord(idxNode,2));

% Loop through each node and build a part with markers at each node
for i = 1:numel(idxNode)
    
    ID      = obj.Mdl.Grid(idxNode(i)).id;
    Coord   = struct2mat(obj.Mdl.Grid,'coord');
    Coord   = Coord(idxNode(i),:);
    idxMass = find([obj.Mdl.Conm2.grid]' == ID);
    
    if ~isempty(idxMass)
        mass    = obj.Mdl.Conm2(idxMass).m;
        inertia = obj.Mdl.Conm2(idxMass).i;
    else
        mass    = 0;
        inertia = zeros(3);
    end
    
    % find orientation of beam and force vectors
    [rotBeamI,rotForceI] = beamOrientVec(obj,ID,'I');
    [rotBeamJ,rotForceJ] = beamOrientVec(obj,ID,'J');
    
    
fprintf(fid,'!-------------------------------- PART_%1.0f ---------------------------------!\n',ID); 
fprintf(fid,'!\n'); 
fprintf(fid,'!\n'); 
fprintf(fid,'defaults coordinate_system  &\n'); 
fprintf(fid,'   default_coordinate_system = .WING_BEAM.ground\n'); 
fprintf(fid,'!\n'); 
fprintf(fid,'part create rigid_body name_and_position  &\n'); 
fprintf(fid,'   part_name = .WING_BEAM.PART_%1.0f  &\n',ID); 
fprintf(fid,'   adams_id = %1.0f  &\n',ID); 
fprintf(fid,'   location = 0.0, 0.0, 0.0  &\n'); 
fprintf(fid,'   orientation = 0.0d, 0.0d, 0.0d\n'); 
fprintf(fid,'!\n'); 
fprintf(fid,'defaults coordinate_system  &\n'); 
fprintf(fid,'   default_coordinate_system = .WING_BEAM.PART_%1.0f\n',ID); 
fprintf(fid,'!\n'); 
fprintf(fid,'! ****** Markers for current part ******\n'); 
fprintf(fid,'!\n'); 
fprintf(fid,'marker create  &\n'); 
fprintf(fid,'   marker_name = .WING_BEAM.PART_%1.0f.MARKER_GRID  &\n',ID);  
fprintf(fid,'   adams_id = %1.0f  &\n',ID*10);  
fprintf(fid,'   location = %6.4f, %6.4f, %6.4f  &\n',Coord); 
fprintf(fid,'   orientation = 0.0d, 0.0d, 0.0d\n'); 
fprintf(fid,'!\n');

if ~isempty(rotBeamI)
fprintf(fid,'marker create  &\n'); 
fprintf(fid,'   marker_name = .WING_BEAM.PART_%1.0f.MARKER_BEAM_I  &\n',ID);  
fprintf(fid,'   adams_id = %1.0f  &\n',ID*10 + 1); 
fprintf(fid,'   location = %6.4f, %6.4f, %6.4f  &\n',Coord); 
fprintf(fid,'   orientation = %6.4fd, %6.4fd, %6.4fd\n',rotBeamI); 
fprintf(fid,'!\n'); 
end

if ~isempty(rotBeamJ)
fprintf(fid,'marker create  &\n'); 
fprintf(fid,'   marker_name = .WING_BEAM.PART_%1.0f.MARKER_BEAM_J  &\n',ID);  
fprintf(fid,'   adams_id = %1.0f  &\n',ID*10 + 2); 
fprintf(fid,'   location = %6.4f, %6.4f, %6.4f  &\n',Coord); 
fprintf(fid,'   orientation = %6.4fd, %6.4fd, %6.4fd\n',rotBeamJ); 
fprintf(fid,'!\n'); 
end

if ~isempty(rotForceI)
fprintf(fid,'marker create  &\n'); 
fprintf(fid,'   marker_name = .WING_BEAM.PART_%1.0f.MARKER_FORCE_I  &\n',ID);  
fprintf(fid,'   adams_id = %1.0f  &\n',ID*10 + 3); 
fprintf(fid,'   location = %6.4f, %6.4f, %6.4f  &\n',Coord); 
fprintf(fid,'   orientation = %6.4fd, %6.4fd, %6.4fd\n',rotForceI); 
fprintf(fid,'!\n'); 
end

if ~isempty(rotForceJ)
fprintf(fid,'marker create  &\n'); 
fprintf(fid,'   marker_name = .WING_BEAM.PART_%1.0f.MARKER_FORCE_J  &\n',ID);  
fprintf(fid,'   adams_id = %1.0f  &\n',ID*10 + 4); 
fprintf(fid,'   location = %6.4f, %6.4f, %6.4f  &\n',Coord); 
fprintf(fid,'   orientation = %6.4fd, %6.4fd, %6.4fd\n',rotForceJ); 
fprintf(fid,'!\n'); 
end

fprintf(fid,'part create rigid_body mass_properties  &\n'); 
fprintf(fid,'   part_name = .WING_BEAM.PART_%1.0f  &\n',ID);  

% Product of inertia terms need to be muliplied by -1 according to
% MSC.ADAMS definition.
Ixx =  abs(inertia(1,1));
Iyy =  abs(inertia(2,2));
Izz =  abs(inertia(3,3));
Ixy = -inertia(1,2);
Izx = -inertia(3,1);
Iyz = -inertia(2,3);

% Inertia checks around primary axes. The sum of the smaller two values
% needs to be equal or larger than the largest value. MSC.ADAMS does
% similar checks. Subtract small term just to make sure it works. There
% seems to be some rounding errors when reading from the command file.
if (Ixx < Iyy) && (Izz < Iyy) && (Ixx + Izz < Iyy)
    Ixx = floor(Ixx * 1000);
    Izz = floor(Izz * 1000);
    Iyy = (Ixx + Izz)/1000 - 0.001;
    Ixx = Ixx/1000;
    Izz = Izz/1000;
elseif (Ixx < Izz) && (Iyy < Izz) && (Ixx + Iyy < Izz)
    %Izz = Ixx + Iyy - 1e-6;
    Ixx = floor(Ixx * 1000);
    Iyy = floor(Iyy * 1000);
    Izz = (Ixx + Iyy)/1000 - 0.001;
    Ixx = Ixx/1000;
    Iyy = Izz/1000;    
elseif (Iyy < Ixx) && (Izz < Ixx) && (Iyy + Izz < Ixx)
    %Ixx = Iyy + Izz - 1e-6;
    Iyy = floor(Iyy * 1000);
    Izz = floor(Izz * 1000);
    Ixx = (Ixx + Izz)/1000 - 0.001;   
    Iyy = Iyy/1000;
    Izz = Izz/1000;    
end

% Check that when decomposed into principal axis that Inertia values are
% all positive.
[~,IM] = eig(inertia);

if any(find(IM < 0))   
   error('A principal axis negative detected for node ID %1.0f',ID);
end

if mass ~=  0
fprintf(fid,'   mass = %18.9f &\n',mass); 
fprintf(fid,'   center_of_mass_marker = .WING_BEAM.PART_%1.0f.MARKER_GRID  &\n',ID);  
fprintf(fid,'   inertia_marker = .WING_BEAM.PART_%1.0f.MARKER_GRID  &\n',ID);  
fprintf(fid,'   ixx = %18.9f   &\n',Ixx); 
fprintf(fid,'   iyy = %18.9f   &\n',Iyy); 
fprintf(fid,'   izz = %18.9f   &\n',Izz); 
fprintf(fid,'   ixy = %18.9f   &\n',Ixy); 
fprintf(fid,'   izx = %18.9f   &\n',Izx); 
fprintf(fid,'   iyz = %18.9f   \n',Iyz); 
fprintf(fid,'!\n'); 
else
fprintf(fid,'   mass = 0.0 &\n'); 
fprintf(fid,'   center_of_mass_marker = .WING_BEAM.PART_%1.0f.MARKER_GRID  &\n',ID);  
fprintf(fid,'   inertia_marker = .WING_BEAM.PART_%1.0f.MARKER_GRID  &\n',ID);  
fprintf(fid,'   ixx = 0.0 &\n'); 
fprintf(fid,'   iyy = 0.0 &\n');  
fprintf(fid,'   izz = 0.0 &\n');  
fprintf(fid,'   ixy = 0.0 &\n');  
fprintf(fid,'   izx = 0.0 &\n');  
fprintf(fid,'   iyz = 0.0 \n');  
fprintf(fid,'!\n');  
end

fprintf(fid,'!\n');
fprintf(fid,'! ****** Graphics for current part ******\n');
fprintf(fid,'!\n');

if   abs(Coord(2)) ~=  wingSemiSpan
idPt1 = ID + 101000;
idPt2 = ID + 101000 + 1;
idPt3 = ID + 102000 + 1;
idPt4 = ID + 102000;

Coord1 = obj.Mdl.Grid(find([obj.Mdl.Grid.id]  ==  idPt1)).coord;
Coord2 = obj.Mdl.Grid(find([obj.Mdl.Grid.id]  ==  idPt2)).coord;
Coord3 = obj.Mdl.Grid(find([obj.Mdl.Grid.id]  ==  idPt3)).coord;
Coord4 = obj.Mdl.Grid(find([obj.Mdl.Grid.id]  ==  idPt4)).coord;

fprintf(fid,'!\n');
fprintf(fid,'geometry create curve polyline  &\n');
fprintf(fid,'   polyline_name = .WING_BEAM.PART_%1.0f.AERO_OUTLINE  &\n',ID);
fprintf(fid,'   location = %d, %d, %d  &\n',Coord1);
fprintf(fid,'      , %d, %d, %d  &\n',Coord2);
fprintf(fid,'      , %d, %d, %d  &\n',Coord3);
fprintf(fid,'      , %d, %d, %d  &\n',Coord4);
fprintf(fid,'      , %d, %d, %d  &\n',Coord1);
fprintf(fid,'   close = yes\n');
fprintf(fid,'!\n');

fprintf(fid,'!\n');
fprintf(fid,'geometry attributes  &\n');
fprintf(fid,'   geometry_name = .WING_BEAM.PART_%1.0f.AERO_OUTLINE  &\n',ID);
fprintf(fid,'   color = CYAN\n');
fprintf(fid,'!\n');
end

fprintf(fid,'part attributes  &\n'); 
fprintf(fid,'   part_name = .WING_BEAM.PART_%1.0f  &\n',ID);  
fprintf(fid,'   color = MAGENTA  &\n'); 
fprintf(fid,'   name_visibility = off\n'); 
fprintf(fid,'!\n'); 

end

end

%% Create Beam elements
% Beam elements are defined as forces in ADAMS. The same definition as used
% in FAME. The E value is set to 1, and the EI and GJ terms are then
% written instead of the A and I terms. The G term is calculated.
function [] = writeForces(fid,obj)

% Obtain area and inertias of beam
coord = struct2mat(obj.Mdl.Grid,'coord');
wingSemiSpan = max(coord(:,2));

% Calculate Inertias along the span
etaBoxAreas = obj.Fame.BoxAreas.eta;
etaBoxStiffness = obj.Fame.BoxStiffness.eta;

% Interpolation function only works on monotonically increasing functions,
% hence add a slight offset to the eta value at points if they overlap
for i = 2:numel(etaBoxAreas)
    if etaBoxAreas(i)  ==  etaBoxAreas(i-1)
        etaBoxAreas(i) = etaBoxAreas(i) + 0.000001;
    end
end

for i = 2:numel(etaBoxStiffness)
    if etaBoxStiffness(i)  ==  etaBoxStiffness(i-1)
        etaBoxStiffness(i) = etaBoxStiffness(i) + 0.000001;
    end
end

fprintf(fid,'!\n');
fprintf(fid,'!----------------------------------- Forces -----------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');

for i = 1:numel(obj.Mdl.Cbar)
    
    % Find ID of beam
    ID = obj.Mdl.Cbar(i).id;
    
    % Obtain coordinates for end point of beam
    idEndPt = obj.Mdl.Cbar(i).conn(2);
    idx = find(idEndPt == [obj.Mdl.Grid.id]);
    coordEndPt = obj.Mdl.Grid(idx).coord;
    
    % Obtain coordinates for start point of beam
    idStartPt = obj.Mdl.Cbar(i).conn(1);
    idx = find(idStartPt == [obj.Mdl.Grid.id]);    
    coordStartPt = obj.Mdl.Grid(idx).coord;
    
    % Calculate length of beam
    beamLength = sqrt((coordEndPt(1) - coordStartPt(1))^2 + ...
                      (coordEndPt(2) - coordStartPt(2))^2 + ...
                      (coordEndPt(3) - coordStartPt(3))^2);
    
    spanPosEta = coordEndPt(2)/wingSemiSpan;
                  
    % Calculate area of material at eta position
    areaMat = interp1(etaBoxAreas,obj.Fame.BoxAreas.a, spanPosEta,'linear','extrap');
    
    EA = interp1(etaBoxStiffness,obj.Fame.BoxStiffness.EA, spanPosEta,'linear','extrap');
    E = 1;
    G = E/(2*(1 + obj.Mdl.Mat.nu));
    
    idx = find(obj.Mdl.Cbar(i).pid == [obj.Mdl.Pbar.id]);

    ixx = obj.Mdl.Pbar(idx).j/(2*(1 + obj.Mdl.Mat.nu));
    iyy = obj.Mdl.Pbar(idx).i(2);
    izz = obj.Mdl.Pbar(idx).i(1);
    
fprintf(fid,'force create element_like beam  &\n');
fprintf(fid,'   beam_name = .WING_BEAM.BEAM_%1.0f  &\n',ID);
fprintf(fid,'   adams_id = %1.0f  &\n',ID);
fprintf(fid,'   i_marker_name = .WING_BEAM.PART_%1.0f.MARKER_BEAM_I  &\n',idEndPt);
fprintf(fid,'   j_marker_name = .WING_BEAM.PART_%1.0f.MARKER_BEAM_J  &\n',idStartPt);
fprintf(fid,'   length = %d  &\n',beamLength);
fprintf(fid,'   area_of_cross_section = %d  &\n',EA);
fprintf(fid,'   y_shear_area_ratio = 1.0  &\n');
fprintf(fid,'   z_shear_area_ratio = 1.0  &\n');
fprintf(fid,'   youngs_modulus = %d  &\n',E);
fprintf(fid,'   shear_modulus = %d  &\n',G);
fprintf(fid,'   ixx = %d  &\n',ixx);
fprintf(fid,'   iyy = %d  &\n',iyy);
fprintf(fid,'   izz = %d  &\n',izz);
fprintf(fid,'   damping_ratio = 2.0E-002  &\n');
fprintf(fid,'   formulation = nonlinear\n');
fprintf(fid,'!\n');
end

end

%% Define Simulation script
function [] = writeSimscript(fid)

fprintf(fid,'!----------------------------- Simulation Scripts -----------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'simulation script create  &\n');
fprintf(fid,'   sim_script_name = .WING_BEAM.Last_Sim  &\n');
fprintf(fid,'   commands =   &\n');
fprintf(fid,'              "simulation single_run transient type = auto_select initial_static = no end_time = 10.0 number_of_steps = 1000 model_name = .WING_BEAM"\n');
fprintf(fid,'!\n');

end

%% Create some panels for visualisation
% A panel is created between two adjacent node points to help with
% visualisation.
function [] = writeDynamicgraphics(fid)

fprintf(fid,'!------------------------------ Dynamic Graphics ------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults coordinate_system  &\n');
fprintf(fid,'   default_coordinate_system = .WING_BEAM.ground\n');
fprintf(fid,'!\n');


end

%% Define Gravity
function [] = writeAcc(fid)

fprintf(fid,'!---------------------------------- Accgrav -----------------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'force create body gravitational  &\n');
fprintf(fid,'   gravity_field_name = gravity  &\n');
fprintf(fid,'   z_component_gravity = -9.80665\n');
fprintf(fid,'!\n');
fprintf(fid,'force attributes  &\n');
fprintf(fid,'   force_name = .WING_BEAM.gravity  &\n');
fprintf(fid,'   visibility = off\n');
fprintf(fid,'!\n');

end

function [] = writeExpressdef(fid)

fprintf(fid,'!--------------------------- Expression definitions ---------------------------!\n');
fprintf(fid,'!\n');
fprintf(fid,'!\n');
fprintf(fid,'defaults coordinate_system  &\n');
fprintf(fid,'   default_coordinate_system = ground\n');
fprintf(fid,'!\n');
fprintf(fid,'material modify  &\n');
fprintf(fid,'   material_name = .WING_BEAM.FAME_Material  &\n');
fprintf(fid,'   density = (1(kg/meter**3))  &\n');
fprintf(fid,'   youngs_modulus = (1(Newton/meter**2))\n');
fprintf(fid,'!\n');
fprintf(fid,'model display  &\n');
fprintf(fid,'   model_name = WING_BEAM\n');

end

%% Calculate beam orientation vectors
function [rotBeam,rotForce] = beamOrientVec(obj,ID,ijNode) 
    
rotBeam  = [0,0,0];
rotForce = [0,0,0];

% Find beam axis vector and work out orientations. Axis orientation is
% determined by using beam axis and vector defining xy-plane
if strcmp(ijNode,'I')
    colNo = 2;
else
    colNo = 1;
end

Conn = struct2mat(obj.Mdl.Cbar,'conn');
iConnId  = find(Conn(:,colNo) ==  ID);

if isempty(iConnId)
    return
end

if numel(iConnId) > 1
    warning('Non-unique node in structural beam definition discovered')
    iConnId = iConnId(1);
end

iEndPnt = find(Conn(iConnId,2) == [obj.Mdl.Grid.id]);
endPnt  = obj.Mdl.Grid(iEndPnt).coord;
iBegPnt = find(Conn(iConnId,1) == [obj.Mdl.Grid.id]);
begPnt  = obj.Mdl.Grid(iBegPnt).coord;

x = endPnt - begPnt;
p = obj.Mdl.Cbar(iConnId).orient;

x = x./(sqrt(dot(x,x)));
p = p./(sqrt(dot(p,p)));

z = cross(p,x);

if z(1) < 0
    z = cross(x,p);
end

y = cross(z,x);

z = z./(sqrt(dot(z,z)));
y = y./(sqrt(dot(y,y)));

dcm = [x;y;z];

[r(1),r(2),r(3)] = dcm2eul(dcm,'ZXZ');

rotBeam = r.*180/pi;

% orientation vector for force node
dcmf = [x;-z;y];

[r(1),r(2),r(3)] = dcm2eul(dcmf,'ZXZ');

rotForce = r.*180/pi;
    
        
end
