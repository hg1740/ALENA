function obj = writeSplineVTP(obj)

NODE  = obj.Mdl.Grid;
CAERO = obj.Mdl.Caero;

% nSets = 0;
discsets = 1; % 1 whole wing - 2 carrythrough / wing - 3 each panel

% fuseHalfWidth(1) = WingCAERO(1).Le(2,2) - WingCAERO(1).Le(1,2);
% fuseHalfWidth(2) = WingCAERO(1).Te(2,2) - WingCAERO(1).Le(1,2);
% fuseHalfWidth  = min(fuseHalfWidth);
%
% if fuseHalfWidth < 1.5 || fuseHalfWidth > 3.6
%  error('fame2mat:Conversion','Fuselage half - width being calculated from planform outline looks odd.\n')
% end
%
% idxStruct = find([NODE.id] < 130000)';
% idxFuse = find(struct2mat(NODE(idxStruct),'coord',2) <= fuseHalfWidth);
% idxWing = idxStruct;
% idxWing(idxFuse) = [];

% ct = (1:51);
% wg = (51:101);
RBF = 2; % Radial basis function spline method 2 = Infinite Plate Spline


% Find the wing nodes

% Find the htp nodes

% Find the vtp nodes
SetId    = 530001;

WingCAERO   = findobj(obj.Mdl.Caero,'part','VTP');

WingNODES   = findobj(obj.Mdl.Grid,'part','VTP','-and','type','Beam');

WingAeroNODES   = findobj(obj.Mdl.Grid,'part','VTP','-and','type','Aerodynamic');

wingnodeidx     = [WingNODES.id];
wingnodeAeroidx = [WingAeroNODES.id];

wingnodeidx = [wingnodeidx,wingnodeAeroidx];

if discsets == 1
    Se = Sets;
    Se.id = SetId;
    Se.data = wingnodeidx';
    Se.part = 'VTP';
    obj.Mdl.Sets = cat(2,obj.Mdl.Sets,Se);
    count = 0;
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
        Sp(count).part = 'VTP';
        
        if ~isempty(WingCAERO(i).csType)
            
            Sp(count).p2 = Sp(count).p2 + WingCAERO(i).csData(4) * WingCAERO(i).ny;
        end
        
    end
    obj.Mdl.Spline = cat(2,obj.Mdl.Spline,Sp);
    
elseif discsets == 2
    error('FAME2MAT:SplineOption','Splining option 2 for carrythrough / wing not implemented yet.')
    %  nSets = 2;
    %  Interp.sets.id(1,1) = 810001;
    %  Interp.sets.id(2,1) = 810002;
    %  Interp.sets.data{1,1} = [NODE.id(1,ct)...
    %   NODE.id(1,INFO.nMgrid + ct(1):INFO.nMgrid + ct(end))...
    %   NODE.id(1,2 * INFO.nMgrid + ct(1):2 * INFO.nMgrid + ct(end))];
    %  Interp.sets.data{2,1} = [NODE.id(1,wg)...
    %   NODE.id(1,INFO.nMgrid + wg(1):INFO.nMgrid + wg(end))...
    %   NODE.id(1,2 * INFO.nMgrid + wg(1):2 * INFO.nMgrid + wg(end))];
    %
    %  for i = 1:INFO.nWingCAERO
    %
    %   [Ywg] = find(NODE.coord(wg(2:end),2)>WingCAERO.geo.startY(i,1) & ...
    %    NODE.coord(wg(2:end),2)<(WingCAERO.geo.startY(i,1) + WingCAERO.geo.b(i,1)));
    %
    %   Interp.spline.id(i,1) = WingCAERO.id(i,1);
    %   Interp.spline.aid(i,1) = WingCAERO.id(i,1);
    %   Interp.spline.p1(i,1) = 1;
    %   Interp.spline.p2(i,1) = WingCAERO.geo.ny(i) * WingCAERO.geo.nx(i);
    %   if isempty(Ywg)
    %    Interp.spline.set(i,1) = Interp.sets.id(1,1);
    %   else
    %    Interp.spline.set(i,1) = Interp.sets.id(2,1);
    %   end
    %   Interp.spline.w(i,1) = RBF;
    %   Interp.spline.rmx(i,1) = 1;
    %   Interp.spline.cond(i,1) = 1e14;
    %  end
elseif discsets == 3
    error('FAME2MAT:SplineOption','Splining option 3 for indivudual panels not implemented yet.')
    %  for i = 1:INFO.nWingCAERO
    %   nSets = nSets + 1;
    %   [y] = find(NODE.coord(1:INFO.nMgrid,2)>WingCAERO.geo.startY(i,1) & ...
    %    NODE.coord(1:INFO.nMgrid,2)<(WingCAERO.geo.startY(i,1) + WingCAERO.geo.b(i,1)));
    %   if max(y) + 1 > INFO.nMgrid
    %    YY = [min(y) - 1;y;...
    %     NODE.id(1,min(y) - 1 + INFO.nMgrid);NODE.id(1,y + INFO.nMgrid)';...
    %     NODE.id(1,min(y) - 1 + 2 * INFO.nMgrid);NODE.id(1,y + 2 * INFO.nMgrid)']; % beam nodes
    %
    %   else
    %    YY = [min(y) - 1;y;max(y) + 1;...
    %     NODE.id(1,min(y) - 1 + INFO.nMgrid);NODE.id(1,y + INFO.nMgrid)';NODE.id(1,max(y) + 1 + INFO.nMgrid);...
    %     NODE.id(1,min(y) - 1 + 2 * INFO.nMgrid);NODE.id(1,y + 2 * INFO.nMgrid)';NODE.id(1,max(y) + 1 + 2 * INFO.nMgrid)];
    %   end
    %   Interp.sets.id(i,1) = nSets;
    %   Interp.sets.data{i,1} = YY';
    %   Interp.spline.id(i,1) = WingCAERO.id(i,1);
    %   Interp.spline.aid(i,1) = WingCAERO.id(i,1);
    %   Interp.spline.p1(i,1) = 1;
    %   Interp.spline.p2(i,1) = WingCAERO.geo.ny(i) * WingCAERO.geo.nx(i);
    %   Interp.spline.set(i,1) = Interp.sets.id(i,1);
    %   Interp.spline.w(i,1) = RBF;
    %   Interp.spline.rmx(i,1) = 1;
    %   Interp.spline.cond(i,1) = 1e14;
    %  end
end

end

