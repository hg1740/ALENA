function obj = defineextramasses(obj)

ReferenceId = obj.Opts.Struct.FuseStartID;

if obj.Opts.Geom.addTail
    
    if isfield(obj.Fame.Geometry,'HTP')
        
        % HTP
        HTPNodes = findobj(obj.Mdl.Grid,'part','StbdHTP','-and','type','Beam');
        % Assume that the HTP weighs 0.9% of the MTOW
        HTPMass = 0.45/100*obj.Fame.Geometry.Weights.mtow;
        npoints = numel(HTPNodes);
        Inertia = zeros(3,3);
        
        for i = 1:npoints
            Conm2Obj        = Conm2;
            Conm2Obj.id     = HTPNodes(i).id;
            Conm2Obj.grid   = HTPNodes(i).id;
            Conm2Obj.cid    = 0;
            Conm2Obj.m      = HTPMass/npoints;
            Conm2Obj.offset = [0,0,0];
            Conm2Obj.i      = Inertia;
            Conm2Obj.type   = 'Structural';
            Conm2Obj.part   = 'StbdHTP';
            obj.Mdl.Conm2   = cat(2,obj.Mdl.Conm2,Conm2Obj);
        end
    end
    if isfield(obj.Fame.Geometry,'VTP')
        % VTP
        VTPNodes = findobj(obj.Mdl.Grid,'part','VTP','-and','type','Beam');
        
        % Assume that the VTP weighs 1.4% of the MTOW
        VTPMass = 0.7/100*obj.Fame.Geometry.Weights.mtow;
        
        npoints = numel(VTPNodes);
        Inertia = zeros(3,3);
        for i = 1:numel(VTPNodes)
            Conm2Obj        = Conm2;
            Conm2Obj.id     = VTPNodes(i).id;
            Conm2Obj.grid   = VTPNodes(i).id;
            Conm2Obj.cid    = 0;
            Conm2Obj.m      = VTPMass/npoints;
            Conm2Obj.offset = [0,0,0];
            Conm2Obj.i      = Inertia;
            Conm2Obj.type   = 'Structural';
            Conm2Obj.part   = 'VTP';
            obj.Mdl.Conm2   = cat(2,obj.Mdl.Conm2,Conm2Obj);
        end
    end
end

WingNodes = findobj(obj.Mdl.Grid,'part','StbdWing','-and','type','Beam');
Wingcoord = vertcat(WingNodes.coord);

if ~isempty(obj.Fame.Geometry.Engine)
    
    % Determine the x and z locations of the wing at the engine y location
    WingBeamX = interp1(Wingcoord(:,2),Wingcoord(:,1),obj.Fame.Geometry.Engine.eta*obj.Fame.Geometry.Wing.span/2);
    WingBeamZ = interp1(Wingcoord(:,2),Wingcoord(:,3),obj.Fame.Geometry.Engine.eta*obj.Fame.Geometry.Wing.span/2);
    
    % Engine
    if isfield(obj.Fame.Geometry,'Engine')
        EngineNode_x = WingBeamX + obj.Fame.Geometry.Engine.pylon_x_loc;
        EngineNode_y = obj.Fame.Geometry.Engine.eta*obj.Fame.Geometry.Wing.span/2;
        EngineNode_z = WingBeamZ + obj.Fame.Geometry.Engine.engine_z_loc;
    end
    
    EngineNode = [EngineNode_x,EngineNode_y,EngineNode_z];
    EngineMass = obj.Fame.Geometry.Engine.engine_mass + obj.Fame.Geometry.Engine.pylon_mass;
    EngineId   = ReferenceId + 3;
    [~,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(EngineMass,5.57,1.5,0.75);
    Inertia = diag([Ixx,Iyy,Izz]);
    % Determine which node along the wing the engine should be associated
    [~,wingidx] = min(abs(Wingcoord(:,2) - EngineNode_y));
    
    EngineOffset = EngineNode - WingNodes(wingidx).coord;
    
    Conm2Obj        = Conm2;
    Conm2Obj.id     = EngineId;
    Conm2Obj.grid   = WingNodes(wingidx).id;
    Conm2Obj.cid    = 0;
    Conm2Obj.m      = EngineMass;
    Conm2Obj.offset = EngineOffset;
    Conm2Obj.i      = Inertia;
    Conm2Obj.type   = 'Engine';
    Conm2Obj.part   = 'StbdWing';
    obj.Mdl.Conm2   = cat(2,obj.Mdl.Conm2,Conm2Obj);
    
end

allconm = findobj(obj.Mdl.Conm2,'type','Engine','-or','type','Structural','-or','type','Secondary');
sum([allconm.m])

if ~isfield(obj.Fame.Geometry.Weights,'oew')
    % Assume that the fuselage mass is 25% of the MTOW
    FuselageMass = 12.5/100*obj.Fame.Geometry.Weights.mtow;
end

if isfield(obj.Fame,'TUXInput') && isfield(obj.Fame.TUXInput.airplane,'fuselage')
    
    FuselageL  = str2double(obj.Fame.TUXInput.airplane.fuselage.length.Attributes.value);
    FuselageW  = str2double(obj.Fame.TUXInput.airplane.fuselage.width.Attributes.value);
    offset = strsplit(obj.Fame.TUXInput.airplane.fuselage.AxisOffset.Attributes.value,',');
    FuselageRP = [str2double(offset{1}),str2double(offset{2}),str2double(offset{3})];
    
else
    if obj.Opts.Geom.addTail
        maxl  = max(vertcat(obj.Mdl.Grid.coord));
    else
        maxl  = 2*max(vertcat(obj.Mdl.Grid.coord));
    end
    FuselageL  = maxl(1);
    FuselageW  = obj.Fame.Geometry.UserVariables.root_wing*obj.Fame.Geometry.Wing.span;
    FuselageRP = [0,0,0];
end

% Fuselage
FuselageId   = ReferenceId + 1;
FuselageNode = FuselageRP + [0.5*FuselageL,0,0];
Inertia      = zeros(3,3);

GridObj       = Grid;
GridObj.id    = FuselageId;
GridObj.cs    = 0;
GridObj.coord = FuselageNode;
GridObj.cd    = 0;
GridObj.ps    = 0;
GridObj.seid  = 0;
GridObj.part  = 'Fuselage';
GridObj.type  = 'Fuselage';
obj.Mdl.Grid  = cat(2,obj.Mdl.Grid,GridObj);

[~,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(FuselageMass,FuselageL,0.5*FuselageW,0.75);

Inertia = diag([Ixx,Iyy,Izz]);

Conm2Obj        = Conm2;
Conm2Obj.id     = FuselageId;
Conm2Obj.grid   = FuselageId;
Conm2Obj.cid    = 0;
Conm2Obj.m      = FuselageMass;
Conm2Obj.offset = [0,0,0];
Conm2Obj.i      = Inertia;
Conm2Obj.type   = 'Structural';
Conm2Obj.part   = 'Fuselage';
obj.Mdl.Conm2   = cat(2,obj.Mdl.Conm2,Conm2Obj);

% Define the SPC at the Fuselage node
SpcObj       = Spc;
SpcObj.id    = FuselageId;
SpcObj.grids = FuselageId;
SpcObj.dof   = 123456;
SpcObj.d     = [];
SpcObj.part  = 'Fuselage';
obj.Mdl.Spc   = cat(2,obj.Mdl.Spc,SpcObj);

% Payload
PayloadNode = [25.0,0,0];
PayloadMass = 16059/2;
PayloadId   = ReferenceId + 2;

[~,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(PayloadMass,16,1.5,0.75);

Inertia = diag([Ixx,Iyy,Izz]);

% Get the Engine offset from the VTP base
PayOffset = PayloadNode - FuselageNode;

Conm2Obj        = Conm2;
Conm2Obj.id     = PayloadId;
Conm2Obj.grid   = FuselageId;
Conm2Obj.cid    = 0;
Conm2Obj.m      = PayloadMass;
Conm2Obj.offset = PayOffset;
Conm2Obj.i      = Inertia;
Conm2Obj.type   = 'Payload';
Conm2Obj.part   = 'Fuselage';
obj.Mdl.Conm2   = cat(2,obj.Mdl.Conm2,Conm2Obj);
PayIndex = numel(obj.Mdl.Conm2);

% Connect the wing to the fuselage node
WingNodes = findobj(obj.Mdl.Grid,'part','StbdWing','-and','type','Beam');

if obj.Opts.Struct.addcarrythrough
    % Add the carrythrough
    Rbe2Obj      = Rbe2;
    Rbe2Obj.id   = 6000;
    Rbe2Obj.gn   = FuselageId;
    Rbe2Obj.dof  = 246;
    Rbe2Obj.gm   = WingNodes(1).id;
    Rbe2Obj.parent = 'Fuselage';
    Rbe2Obj.child  = 'StbdWing';
    obj.Mdl.Rbe2 = Rbe2Obj;
    
    % Extract the wing root eta value to assign the RBE2
    wing_root   = str2double(cell2mat(obj.Fame.fm4Input.USER_VARIABLES.root_wing));
    wing_span   = str2double(cell2mat(obj.Fame.fm4Input.WING.GLOBAL_GEOMETRY.SPAN))/1000;
    root_loc    = wing_root * wing_span/2;
    Ynle        = struct2mat(WingNodes,'coord',2);
    [~,rootidx] = min(abs(Ynle - root_loc));
    
    Rbe2Obj      = Rbe2;
    Rbe2Obj.id   = 6001;
    Rbe2Obj.gn   = FuselageId;
    Rbe2Obj.dof  = 135;
    Rbe2Obj.gm   = WingNodes(rootidx).id;
    Rbe2Obj.parent = 'Fuselage';
    Rbe2Obj.child  = 'StbdWing';
    obj.Mdl.Rbe2 = cat(2,obj.Mdl.Rbe2,Rbe2Obj);
    
else
    
    % Do not add the carrythrough
    Rbe2Obj      = Rbe2;
    Rbe2Obj.id   = 6000;
    Rbe2Obj.gn   = FuselageId;
    Rbe2Obj.dof  = 123456;
    Rbe2Obj.gm   = WingNodes(1).id;
    Rbe2Obj.parent = 'Fuselage';
    Rbe2Obj.child  = 'StbdWing';
    obj.Mdl.Rbe2 = Rbe2Obj;
end

if obj.Opts.Geom.addTail
    if isfield(obj.Fame.Geometry,'HTP')
        % Connect the HTP to the Fuselage
        Rbe2Obj      = Rbe2;
        Rbe2Obj.id   = 6100;
        Rbe2Obj.gn   = FuselageId;
        Rbe2Obj.dof  = 123456;
        Rbe2Obj.gm   = HTPNodes(1).id;
        Rbe2Obj.parent = 'Fuselage';
        Rbe2Obj.child  = 'StbdHTP';
        obj.Mdl.Rbe2 = cat(2,obj.Mdl.Rbe2,Rbe2Obj);
    end
    if isfield(obj.Fame.Geometry,'VTP')
        % Conect the VTP to the fuselage
        Rbe2Obj      = Rbe2;
        Rbe2Obj.id   = 6200;
        Rbe2Obj.gn   = FuselageId;
        Rbe2Obj.dof  = 123456;
        Rbe2Obj.gm   = VTPNodes(1).id;
        Rbe2Obj.parent = 'Fuselage';
        Rbe2Obj.child  = 'VTP';
        obj.Mdl.Rbe2 = cat(2,obj.Mdl.Rbe2,Rbe2Obj);
    end
end
% % Connect the Engine to the VTP base
% Rbe2Obj     = Rbe2;
% Rbe2Obj.id  = 6300;
% Rbe2Obj.gn  = VTPNodes(1).id;
% Rbe2Obj.dof = 123456;
% Rbe2Obj.gm  = EngineId;
% Rbe2Obj.parent = 'VTP';
% Rbe2Obj.child  = 'Engine';
% obj.Mdl.Rbe2 = cat(2,obj.Mdl.Rbe2,Rbe2Obj);

% % Connect the Payload to the fuselage
% Rbe2Obj      = Rbe2;
% Rbe2Obj.id   = 6400;
% Rbe2Obj.gn   = FuselageId;
% Rbe2Obj.dof  = 123456;
% Rbe2Obj.gm   = PayloadId;
% Rbe2Obj.parent = 'Fuselage';
% Rbe2Obj.child  = 'Payload';
% obj.Mdl.Rbe2 = cat(2,obj.Mdl.Rbe2,Rbe2Obj);

allconm = findobj(obj.Mdl.Conm2,'type','Engine','-or','type','Structural','-or','type','Secondary');
mzfw = str2double(obj.Fame.fm4Input.USER_VARIABLES.mzfw)/2;

obj.Mdl.Conm2(PayIndex).m = mzfw - sum([allconm.m]);

obj.Fame.Geometry.Weights.oew = 2*(mzfw - obj.Mdl.Conm2(PayIndex).m);

end