function obj = writeSpline(obj,part)

RBF             = 2; % Radial basis function spline method 2 = Infinite Plate Spline
PartId          = [part,'SetId'];
SetId           = obj.Opts.Struct.(PartId);

if strcmpi(part,'vtp')
    PartName = part;
else
    PartName = ['Stbd',part];
end

WingCAERO   = findobj(obj.Mdl.Caero,'part',PartName);

WingNODES   = findobj(obj.Mdl.Grid,'part',PartName,'-and','type','Beam');

WingAeroNODES   = findobj(obj.Mdl.Grid,'part',PartName,'-and','type','Aerodynamic');

wingnodeidx     = [WingNODES.id];
wingnodeAeroidx = [WingAeroNODES.id];

wingnodeidx = [wingnodeidx,wingnodeAeroidx];

Se              = Sets;
Se.id           = SetId;
Se.data         = wingnodeidx';
Se.part         = part;

if isempty(obj.Mdl.Sets(1).id)
    obj.Mdl.Sets = Se;
else
    obj.Mdl.Sets = cat(2,obj.Mdl.Sets,Se);
end
count           = 0;

for i = 1:numel(WingCAERO)
    count = count +1;
    Sp(count) = Spline;
    
    Sp(count).id   = WingCAERO(i).id;
    Sp(count).aid  = WingCAERO(i).id;
    Sp(count).p1   = 1;
    Sp(count).p2   = WingCAERO(i).ny * WingCAERO(i).nx;
    Sp(count).w    = RBF;
    Sp(count).rmx  = 1;
    Sp(count).cond = 1e14;
    Sp(count).set  = Se.id;
    Sp(count).part = PartName;
    
    if ~isempty(WingCAERO(i).csType)
        Sp(count).p2 = Sp(count).p2 + WingCAERO(i).csData(4) * WingCAERO(i).ny;
    end
    
end

if isempty(obj.Mdl.Spline(1).id)
    obj.Mdl.Spline = setFields(obj.Mdl.Spline,Sp);
else
    obj.Mdl.Spline = cat(2,obj.Mdl.Spline,Sp);
end


end

