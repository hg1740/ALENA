%% SYMWINGTHRUST Reflect thrust definition around XZ-plane
function obj = symWingThrust( obj )

Offsetidx = obj.Opts.Struct.StrucOffsetID;

nThrust = numel(obj.Mdl.Thrust);

for i = 1:nThrust
    ThrustObj = Thrust;
    
    ThrustObj.lid = obj.Mdl.Thrust(i).lid + Offsetidx;
    ThrustObj.g   = obj.Mdl.Thrust(i).g + Offsetidx;
    ThrustObj.cid = obj.Mdl.Thrust(i).cid;
    ThrustObj.cx  = obj.Mdl.Thrust(i).cx;
    ThrustObj.cy  = obj.Mdl.Thrust(i).cy;
    ThrustObj.cz  = obj.Mdl.Thrust(i).cz;
    ThrustObj.ox  = obj.Mdl.Thrust(i).ox;
    ThrustObj.oy  = obj.Mdl.Thrust(i).oy;
    ThrustObj.oz  = obj.Mdl.Thrust(i).oz;
    
    obj.Mdl.Thrust(nThrust + i) = ThrustObj;
end

end
