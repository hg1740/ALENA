function obj = generateEmpStructure(obj)

%HTP - THIS HARDCODED FOR THE MINUTE!
if isfield(obj.Fame.Geometry,'HTP')
    startNumber = obj.Opts.Struct.HTPStartID;
    Outline = obj.Fame.Geometry.HTP.Outline;
    
    numsec = size(Outline,1)/2;
    HTPChord = flipud(Outline(numsec+1:end,1))- Outline(1:numsec,1);
    
    % Place it at 40%
    obj.Fame.Geometry.HTP.fspar         = [0;1];
    obj.Fame.Geometry.HTP.fspar(:,2)    = [0.15;0.15];
    obj.Fame.Geometry.HTP.aspar         = [0;1];
    obj.Fame.Geometry.HTP.aspar(:,2)    = [0.65;0.65];
    
    HTPNodes = zeros(numsec,3);
    HTPNodes(:,1) = Outline(1:numsec,1) + 0.4*HTPChord;
    HTPNodes(:,2) = Outline(1:numsec,2);
    HTPNodes(:,3) = Outline(1:numsec,3);
    
    ngrid = size(HTPNodes,1);
    GridObj = Grid;
    for i = 1:ngrid
        GridObj(i)       = Grid;
        GridObj(i).id    = i + startNumber;
        GridObj(i).cs    = 0;
        GridObj(i).coord = HTPNodes(i,:);
        GridObj(i).cd    = 0;
        GridObj(i).ps    = 0;
        GridObj(i).seid  = 0;
        GridObj(i).part = 'StbdHTP';
        GridObj(i).type = 'Beam';
    end
    
    % Add the ids to the part object
    PartObj      = PartId;
    CbarObj      = Cbar;
    PbarObj      = Pbar;
    PartObj.id   = 1;
    PartObj.part = 'StbdHTP';
    PartObj.type = 'GRID';
    PartObj.data = [GridObj.id];
    
    for i = 1:ngrid-1
        
        PbarObj(i)     = Pbar;
        PbarObj(i).id  = i + startNumber;
        PbarObj(i).mat = 1;
        PbarObj(i).a   = 1e14;
        PbarObj(i).i   = [1e14,1e14,0];
        PbarObj(i).j   = 1e14;
        PbarObj(i).part = 'StbdHTP';
        
        CbarObj(i)        = Cbar;
        CbarObj(i).id     = i + startNumber;
        CbarObj(i).pid    = i + startNumber;
        CbarObj(i).conn   = [GridObj(i).id,GridObj(i + 1).id];
        CbarObj(i).orient = [0,0,1];
        CbarObj(i).offset = 'GGG';
        CbarObj(i).part = 'StbdHTP';
    end
    
    % Add the ids to the part object
    PartObj(2).id   = 1;
    PartObj(2).part = 'StbdHTP';
    PartObj(2).type = 'CBAR';
    PartObj(2).data = [CbarObj.id];
    
    % Assign to model object
    obj.Mdl.Grid    = cat(2,obj.Mdl.Grid,GridObj);
    obj.Mdl.Cbar    = cat(2,obj.Mdl.Cbar,CbarObj);
    obj.Mdl.Pbar    = cat(2,obj.Mdl.Pbar,PbarObj);
    obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);
end
if isfield(obj.Fame.Geometry,'VTP')
    
    %VTP - THIS HARDCODED FOR THE MINUTE!
    startNumber = obj.Opts.Struct.VTPStartID;
    Outline = obj.Fame.Geometry.VTP.Outline;
    numsec = size(Outline,1)/2;
    VTPChord = flipud(Outline(numsec+1:end,1))- Outline(1:numsec,1);
    % Place it at 40%
    obj.Fame.Geometry.VTP.fspar         = [0;1];
    obj.Fame.Geometry.VTP.fspar(:,2)    = [0.15;0.15];
    obj.Fame.Geometry.VTP.aspar         = [0;1];
    obj.Fame.Geometry.VTP.aspar(:,2)    = [0.65;0.65];
    
    VTPNodes = zeros(numsec,3);
    VTPNodes(:,1) = Outline(1:numsec,1) + 0.4*VTPChord;
    VTPNodes(:,2) = Outline(1:numsec,2);
    VTPNodes(:,3) = Outline(1:numsec,3);
    
    ngrid = size(VTPNodes,1);
    GridObj = Grid;
    for i = 1:ngrid
        GridObj(i)       = Grid;
        GridObj(i).id    = i + startNumber;
        GridObj(i).cs    = 0;
        GridObj(i).coord = VTPNodes(i,:);
        GridObj(i).cd    = 0;
        GridObj(i).ps    = 0;
        GridObj(i).seid  = 0;
        GridObj(i).part = 'VTP';
        GridObj(i).type = 'Beam';
    end
    
    % Add the ids to the part object
    PartObj      = PartId;
    PbarObj      = Pbar;
    CbarObj      = Cbar;
    PartObj.id   = 1;
    PartObj.part = 'VTP';
    PartObj.type = 'GRID';
    PartObj.data = [GridObj.id];
    
    for i = 1:ngrid-1
        
        PbarObj(i)     = Pbar;
        PbarObj(i).id  = i + startNumber;
        PbarObj(i).mat = 1;
        PbarObj(i).a   = 1e14;
        PbarObj(i).i   = [1e14,1e14,0];
        PbarObj(i).j   = 1e14;
        PbarObj(i).part = 'VTP';
        
        CbarObj(i)        = Cbar;
        CbarObj(i).id     = i + startNumber;
        CbarObj(i).pid    = i + startNumber;
        CbarObj(i).conn   = [GridObj(i).id,GridObj(i + 1).id];
        CbarObj(i).orient = [0,-1,0];
        CbarObj(i).offset = 'GGG';
        CbarObj(i).part = 'VTP';
    end
    
    % Add the ids to the part object
    PartObj(2).id   = 1;
    PartObj(2).part = 'VTP';
    PartObj(2).type = 'CBAR';
    PartObj(2).data = [CbarObj.id];
    
    % Assign to model object
    obj.Mdl.Grid    = cat(2,obj.Mdl.Grid,GridObj);
    obj.Mdl.Cbar    = cat(2,obj.Mdl.Cbar,CbarObj);
    obj.Mdl.Pbar    = cat(2,obj.Mdl.Pbar,PbarObj);
    obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);
    
end
end