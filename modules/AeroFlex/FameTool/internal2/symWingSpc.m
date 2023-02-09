%% SYMWINGSPC Reflect spc definition around XZ-plane
function obj = symWingSpc(obj)

Offsetidx = obj.Opts.Struct.StrucOffsetID;

nSpc = numel(obj.Mdl.Spc);

gridPts = struct2mat(obj.Mdl.Grid,'id');
gridCoordY = struct2mat(obj.Mdl.Grid,'coord',2);

for i = 1:nSpc
    SpcObj = Spc;
    
    idx = find(obj.Mdl.Spc(i).grids == gridPts);

    if abs(gridCoordY(idx) < 1e-6)
        continue
    end  
    
    SpcObj.id    = obj.Mdl.Spc(i).id;
    SpcObj.grids = obj.Mdl.Spc(i).grids + Offsetidx;
    SpcObj.dof   = obj.Mdl.Spc(i).dof;
    SpcObj.d     = obj.Mdl.Spc(i).d;
    
    obj.Mdl.Spc = cat(2,obj.Mdl.Spc,SpcObj);
end

end

