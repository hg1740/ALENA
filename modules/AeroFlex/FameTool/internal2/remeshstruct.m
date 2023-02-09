function obj = remeshstruct(obj)

eta = obj.Opts.Struct.remeshStruct.eta;
elm = obj.Opts.Struct.remeshStruct.elementNumbers;

% Build vectors of all the grid and element data
x        = struct2mat(obj.Mdl.Grid,'coord',1);
y        = struct2mat(obj.Mdl.Grid,'coord',2);
z        = struct2mat(obj.Mdl.Grid,'coord',3);
gridId   = [obj.Mdl.Grid.id]';
cbarConn = cell2mat({obj.Mdl.Cbar.conn}');
pbarId   = cell2mat({obj.Mdl.Cbar.pid}');

halfSpan = max(y);

etaOld = y / halfSpan;

% Determine new eta values
etaNew = [];
tol    = 0;

for i = 1:numel(elm)
    ds     = (eta(i + 1) - eta(i)) / elm(i);
    etaNew = cat(2,etaNew,eta(i) + tol:ds:eta(i + 1) - tol);
end

% Rebuild grid elements
etaNew = unique(etaNew);
startNumber = 110000;

for i = 1:numel(etaNew)
    xi = interp1(etaOld,x,etaNew(i),'linear','extrap');
    yi = interp1(etaOld,y,etaNew(i),'linear','extrap');
    zi = interp1(etaOld,z,etaNew(i),'linear','extrap');
    
    GridObj(i)       = Grid;
    GridObj(i).id    = i + startNumber;
    GridObj(i).cs    = 0;
    GridObj(i).coord = [xi, yi, zi];
    GridObj(i).cd    = 0;
    GridObj(i).ps    = 0;
    GridObj(i).seid  = 0;
    GridObj(i).part  = 'Wing';
    GridObj(i).type  = 'beam';
end

% Rebuild bar and property elements
for i = 1:numel(etaNew) - 1
    idxHigher = find(etaOld > etaNew(i));
    idxLower  = find(etaOld < etaNew(i + 1));
    
    [idxGrid,~,~] = intersect(idxHigher,idxLower);
    
    gridIdSel = gridId(idxGrid);
    
    [~,~,idxCbarInboard]  = intersect(gridIdSel,cbarConn(:,1));
    [~,~,idxCbarOutboard] = intersect(gridIdSel,cbarConn(:,2));
    
    idxCbar     = unique([idxCbarInboard;idxCbarOutboard]);
    iMom        = [0,0,0];
    jTor        = 0;
    area        = 0;
    orient      = [0,0,0];
    numElements = numel(idxCbar);
    
    for j = 1:numElements
        idxPid = find(obj.Mdl.Cbar(idxCbar(j)).pid == pbarId);
        iMom   = obj.Mdl.Pbar(idxPid(1)).i + iMom;
        jTor   = obj.Mdl.Pbar(idxPid(1)).j + jTor;
        area   = obj.Mdl.Pbar(idxPid(1)).a + area;
        orient = obj.Mdl.Cbar(idxPid(1)).orient + orient;
    end
    
    PbarObj(i)     = Pbar;
    PbarObj(i).id  = i + startNumber;
    PbarObj(i).mat = 1;
    PbarObj(i).a   = area / numElements;
    PbarObj(i).i   = iMom / numElements;
    PbarObj(i).j   = jTor / numElements;
    PbarObj(i).part  = 'Wing';
    
    CbarObj(i)        = Cbar;
    CbarObj(i).id     = i + startNumber;
    CbarObj(i).pid    = i + startNumber;
    CbarObj(i).conn   = [GridObj(i).id,GridObj(i + 1).id];
    CbarObj(i).orient = orient / numElements;
    CbarObj(i).offset = 'GGG';
    CbarObj(i).part  = 'Wing';
end

% Rebuild spc elements
Xn = struct2mat(GridObj,'coord',1);
Yn = struct2mat(GridObj,'coord',2);
Zn = struct2mat(GridObj,'coord',3);

for i = 1:numel(obj.Mdl.Spc)
    
    idxSpc = find(obj.Mdl.Spc(i).grids == gridId);
    
    Xo = repmat(x(idxSpc(1)),size(Xn));
    Yo = repmat(y(idxSpc(1)),size(Yn));
    Zo = repmat(z(idxSpc(1)),size(Zn));
    
    distance = sqrt((Xo - Xn).^2 + (Yo - Yn).^2 + (Zo - Zn).^2);
    
    [~,idxMinDist] = min(distance);
    
    obj.Mdl.Spc(i).grids = GridObj(idxMinDist).id;
end


% Assign to model object
obj.Mdl.Grid = GridObj;
obj.Mdl.Cbar = CbarObj;
obj.Mdl.Pbar = PbarObj;

end

