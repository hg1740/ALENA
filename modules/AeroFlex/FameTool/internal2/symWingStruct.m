%% SYMWINGSTRUCT Reflect wing structure definition around XZ-plane

function obj = symWingStruct(obj)

offsetidx = obj.Opts.Struct.StrucOffsetID;

% Lets define an array of names
Names = {'StbdWing','PortWing';'StbdHTP','PortHTP'};

%% REFLECT THE GRID NODES %%

% Let's begin by reflecting the grid nodes
nGrid = numel(obj.Mdl.Grid);

for i = 1:nGrid
    
    GridObj = Grid;
    
    if abs(obj.Mdl.Grid(i).coord(2)) < 1e-6
        continue        
    end
    
    GridObj.id    = obj.Mdl.Grid(i).id + offsetidx;
    GridObj.cs    = obj.Mdl.Grid(i).cs;
    GridObj.coord = obj.Mdl.Grid(i).coord .* [1,-1,1];
    GridObj.cd    = obj.Mdl.Grid(i).cd;
    GridObj.ps    = obj.Mdl.Grid(i).ps;
    GridObj.seid  = obj.Mdl.Grid(i).seid;
    GridObj.type  = obj.Mdl.Grid(i).type;
    
    % Lets assign the port name to the grid nodes
    [~,ic]        = ismember(obj.Mdl.Grid(i).part,Names);
    if ic~=0
        GridObj.part  = Names{ic,2};
    else
        GridObj.part  = obj.Mdl.Grid(i).part;
    end
    % Add the new grid objects to the model
    obj.Mdl.Grid  = cat(2,obj.Mdl.Grid,GridObj);
end


% We have to check for the stbdwing and portwing because the center node
% still corresponds only to the stbdWing.
GridObjWing     = findobj(obj.Mdl.Grid,'part','StbdWing','-or','part','PortWing');
GridObj         = findobj(GridObjWing,'type','Beam');
GridObjCoord    = vertcat(GridObj.coord);
GridObjId       = [GridObj.id]';

% Now lets' isolate the nodes that belong to the port wing beams
LWingId      = GridObjId(GridObjCoord(:,2)<=0);

%% ADD GRIDS TO THE PART OBJECTS %%

% Add the IDs to the the Part Object.
PartObj      = PartId;
PartObj.id   = 1;
PartObj.part = 'PortWing';
PartObj.type = 'GRID';
PartObj.data = LWingId';

% Add to the rest of the part objects
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

% Now lets deal with the HTP grids:
GridObjHTP = findobj(obj.Mdl.Grid,'part','StbdHTP','-or','part','PortHTP');
GridObj    = findobj(GridObjHTP,'type','Beam');

if ~isempty(GridObj)

    GridObjCoord = vertcat(GridObj.coord);
    GridObjId    = [GridObj.id]';
    
    % Now lets' isolate the nodes that belong to the port HTP beams
    LHTPId       = GridObjId(GridObjCoord(:,2)<=0);
    
    % Add the ids to the part object
    PartObj      = PartId;
    PartObj.id   = 1;
    PartObj.part = 'PortHTP';
    PartObj.type = 'GRID';
    PartObj.data = LHTPId';
end

% Add to the rest of the part objects:
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

%% REFLECT THE PBARS %%

% Let's remind ourselves of how many cbars currently exist in the workspace
nPbar = numel(obj.Mdl.Pbar);
for i = 1:nPbar
    if ~(strcmpi(obj.Mdl.Pbar(i).part,'VTP'))
        
        PbarObj(i)     = Pbar;
        PbarObj(i).id  = obj.Mdl.Pbar(i).id + offsetidx;
        PbarObj(i).mat = obj.Mdl.Pbar(i).mat;
        PbarObj(i).a   = obj.Mdl.Pbar(i).a;
        PbarObj(i).i   = obj.Mdl.Pbar(i).i.*[1,1,-1];
        PbarObj(i).j   = obj.Mdl.Pbar(i).j;
        [~,ic]        = ismember(obj.Mdl.Pbar(i).part,Names);
        if ic~=0
            PbarObj(i).part  = Names{ic,2};
        else
            PbarObj(i).part  = obj.Mdl.Pbar(i).part;
        end
    end
end

% Add the new pbar objects to the model
obj.Mdl.Pbar  = cat(2,obj.Mdl.Pbar,PbarObj);

%% REFLECT THE CBARS %%

% Let's remind ourselves of how many cbars currently exist in the workspace
nCbar = numel(obj.Mdl.Cbar);

% Get the ids and the y coordinates of the grid nodes
gridPts    = struct2mat(obj.Mdl.Grid,'id');
gridCoordY = struct2mat(obj.Mdl.Grid,'coord',2);

for i = 1:nCbar
    if ~(strcmpi(obj.Mdl.Cbar(i).part,'VTP') || strcmpi(obj.Mdl.Cbar(i).part,'fuselage'))
        
        CbarObj = Cbar;
        CbarObj.id  = obj.Mdl.Cbar(i).id + offsetidx;
        CbarObj.pid = obj.Mdl.Cbar(i).pid + offsetidx;
        
        % Lets find the start node of the cbar        
        idxFirstNode  = find(obj.Mdl.Cbar(i).conn(1) == gridPts);
        
        % If the start node is close enough to y == 0 the node has not been
        % reflected
        if abs(gridCoordY(idxFirstNode) < 1e-6)
            firstNodeId  = obj.Mdl.Cbar(i).conn(1);
        else
            firstNodeId  = obj.Mdl.Cbar(i).conn(1) + offsetidx;
        end
        
        % Lets find the end node of the cbar 
        idxSecondNode = find(obj.Mdl.Cbar(i).conn(2) == gridPts);
        
        % If the start node is close enough to y == 0 the node has not been
        % reflected
        if abs(gridCoordY(idxSecondNode) < 1e-6)
            secondNodeId  = obj.Mdl.Cbar(i).conn(2);
        else
            secondNodeId  = obj.Mdl.Cbar(i).conn(2) + offsetidx;
        end
        
        % Let's store the connection nodes
        CbarObj.conn   = [firstNodeId secondNodeId];
        CbarObj.barg0  = obj.Mdl.Cbar(i).barg0;
        
        % Let's account for the change in orientation when reflecting
        CbarObj.orient = obj.Mdl.Cbar(i).orient .* [1,-1,1];
        
        % Account for some numerical issues with orientation vector
        idx = abs(CbarObj.orient) < 1e-5;
        CbarObj.orient(idx) = 0;
        
        CbarObj.offset = obj.Mdl.Cbar(i).offset;
        
        [~,ic] = ismember(obj.Mdl.Cbar(i).part,Names);
        if ic~=0
            CbarObj.part  = Names{ic,2};
        else
            CbarObj.part  = obj.Mdl.Cbar(i).part;
        end
        obj.Mdl.Cbar = cat(2,obj.Mdl.Cbar,CbarObj);
    end
end

%% ADD CBARS TO THE PART OBJECTS %%

% Isolate the Wing grids:
CbarObj = findobj(obj.Mdl.Cbar,'part','PortWing');

% Isolate the LHS Wing
CbarObjNode = vertcat(CbarObj.conn);

[~,NodeIdx] = intersect(CbarObjNode(:,2),LWingId);

CbarObjId   = [CbarObj.id]';

LWingId     = CbarObjId(NodeIdx);

% Add the ids to the part object
PartObj      = PartId;
PartObj.id   = 1;
PartObj.part = 'PortWing';
PartObj.type = 'CBAR';
PartObj.data = LWingId';

% Add to the objects:
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

if obj.Opts.Geom.addTail
    % Isolate the HTP grids:
    CbarObj = findobj(obj.Mdl.Cbar,'part','PortHTP');
    
    % Isolate the LHS Wing
    CbarObjNode = vertcat(CbarObj.conn);
    
    [~,NodeIdx] = intersect(CbarObjNode(:,2),LHTPId);
    
    CbarObjId   = [CbarObj.id]';
    
    LHTPId      = CbarObjId(NodeIdx);
    
    % Add the ids to the part object
    PartObj      = PartId;
    PartObj.id   = 1;
    PartObj.part = 'PortHTP';
    PartObj.type = 'CBAR';
    PartObj.data = LHTPId';
end

% Add to the objects:
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

%% REFLECT THE LUMPED MASS %%

% Reflect the masses
nConm2 = numel(obj.Mdl.Conm2);

for i = 1:nConm2
    Conm2Obj = Conm2;
    
    idxStart = find(obj.Mdl.Conm2(i).grid == gridPts);
    
    if abs(gridCoordY(idxStart) < 1e-6)
        Conm2Obj.grid   = obj.Mdl.Conm2(i).grid;
    else
        Conm2Obj.grid   = obj.Mdl.Conm2(i).grid + offsetidx;
    end
    Conm2Obj.id     = obj.Mdl.Conm2(i).id + offsetidx;
    Conm2Obj.cid    = obj.Mdl.Conm2(i).cid;
    Conm2Obj.m      = obj.Mdl.Conm2(i).m;
    Conm2Obj.offset = obj.Mdl.Conm2(i).offset .* [1,-1,1];
    Conm2Obj.i      = obj.Mdl.Conm2(i).i  .* [1,-1,1;-1,1,-1;1,-1,1];
    Conm2Obj.type   = obj.Mdl.Conm2(i).type;
    Conm2Obj.file   = obj.Mdl.Conm2(i).file;
    [~,ic]        = ismember(obj.Mdl.Conm2(i).part,Names);
    if ic~=0
        Conm2Obj.part  = Names{ic,2};
    else
        Conm2Obj.part  = obj.Mdl.Conm2(i).part;
    end
    obj.Mdl.Conm2 = cat(2,obj.Mdl.Conm2,Conm2Obj);

end
%% REFLECT THE RBE0s, AERODYNAMIC NODE CONNECTIONS %%

% Reflect the RBE0s
nRbe0 = numel(obj.Mdl.Rbe0);

for i = 1:nRbe0
    if ~(strcmpi(obj.Mdl.Rbe0(i).part,'VTP') || strcmpi(obj.Mdl.Rbe0(i).part,'fuselage'))
        idxStart = find(obj.Mdl.Rbe0(i).data(1) == gridPts);
        idxEnd   = find(obj.Mdl.Rbe0(i).data(2) == gridPts);
        
        if abs(gridCoordY(idxStart) < 1e-6) || abs(gridCoordY(idxEnd) < 1e-6)
            continue
        end
        
        Rbe0Obj = Rbe0;
        
        firstNodeId  = obj.Mdl.Rbe0(i).data(1) + offsetidx;
        secondNodeId = obj.Mdl.Rbe0(i).data(2) + offsetidx;
        
        Rbe0Obj.id   = obj.Mdl.Rbe0(i).id + offsetidx;
        Rbe0Obj.grid = obj.Mdl.Rbe0(i).grid + offsetidx;
        Rbe0Obj.data = [firstNodeId;secondNodeId];
        Rbe0Obj.part = obj.Mdl.Rbe0(i).part;
        obj.Mdl.Rbe0 = cat(2,obj.Mdl.Rbe0,Rbe0Obj);
    end
end

%% REFLECT THE RBE2s, TRADITIONAL RIGID NODE CONNECTIONS %%

% Reflect the engine RBE2
% If the master is the fuselage or the VTP do not reflect it
% If the slave is the fuselage or VTP do not reflect the slave
nRbe2 = numel(obj.Mdl.Rbe2);

count = 0;
for i = 1:nRbe2
   
    
    OffsetMaster = offsetidx;
    OffsetSlave  = offsetidx;
    
    idxMaster = find(obj.Mdl.Rbe2(i).gn == gridPts);
    idxSlave  = find(obj.Mdl.Rbe2(i).gm == gridPts);
    
    if abs(gridCoordY(idxMaster) < 1e-6)
        OffsetMaster = 0;
    end
    
    if abs(gridCoordY(idxSlave) < 1e-6)
        OffsetSlave = 0;
    end
    
    if OffsetMaster + OffsetSlave >0
        Rbe2Obj = Rbe2;
        count = count + 1;
        Rbe2Obj.id  = obj.Mdl.Rbe2(i).id + offsetidx;
        Rbe2Obj.gn  = obj.Mdl.Rbe2(i).gn + OffsetMaster;
        Rbe2Obj.dof = obj.Mdl.Rbe2(i).dof;
        Rbe2Obj.gm  = obj.Mdl.Rbe2(i).gm + OffsetSlave;
        Rbe2Obj.parent = obj.Mdl.Rbe2(i).parent;
        Rbe2Obj.child  = obj.Mdl.Rbe2(i).child;
        obj.Mdl.Rbe2(nRbe2 + count) = Rbe2Obj;
    end
end


end
