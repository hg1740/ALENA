%% CONTROLSURFACEGENERATOR   Generates aerodynamic panel defintion of control
%                            surfaces for NEOCASS and NASTRAN.
%
%  CONTROLSURFACEGENERATOR(CLASSNAME) adds control surface defintions into
%  existing |Caero| objects.
%
%  CLASSNAME is a Fame2mat object
%
%  This is a private function and should be accessed from the |Fame2mat|
%  object.
%
%   Example:
%      obj = fame2mat;
%      obj = setControlSurfaces(obj);   % calls this function
%
%   Restrictions:
%       NOTE: Only defines ailerons for now. Other control surfaces need to
%             be added later.
%
%   See also |setControlSurfaces| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function obj = controlSurfaceGeneratorVTP(obj)

CAERO     = findobj(obj.Mdl.Caero,'part','VTP');
WingNodes = findobj(obj.Mdl.Grid,'part','VTP');
wingaeroidx = find(ismember([obj.Mdl.Caero.id],[CAERO.id]));

Fame     = obj.Fame;
maxeta = max(struct2mat(WingNodes,'coord',3));
mineta = min(struct2mat(WingNodes,'coord',3));

halfSpan = maxeta - mineta;

% Specify number of cs panels
% TODO: This hard code needs to be removed so that panels set from object definition
[~,InboardCutLineIndex]  = min(abs(struct2mat(CAERO,'startZ') - (halfSpan * Fame.Geometry.VTP.ControlSurfaces.rudder.eta_beg_le) + mineta));
[~,OutboardCutLineIndex] = min(abs(struct2mat(CAERO,'startZ') + struct2mat(CAERO,'b') - (halfSpan * Fame.Geometry.VTP.ControlSurfaces.rudder.eta_end_le)));
count = 0;
aerocount  = 0;
ref_point = [obj.Mdl.Caero(wingaeroidx(1)).startX,obj.Mdl.Caero(wingaeroidx(1)).startY,obj.Mdl.Caero(wingaeroidx(1)).startZ];

for i = wingaeroidx(1:end-1)
    aerocount = aerocount + 1;
    if aerocount >= InboardCutLineIndex && aerocount <= OutboardCutLineIndex
        count = count + 1;
        eta   = [Fame.Geometry.VTP.ControlSurfaces.rudder.eta_beg_le,Fame.Geometry.VTP.ControlSurfaces.rudder.eta_end_le];
        chord = [Fame.Geometry.VTP.ControlSurfaces.rudder.chrd_beg_le,Fame.Geometry.VTP.ControlSurfaces.rudder.chrd_end_le];
        eib   = (obj.Mdl.Caero(i).startZ - ref_point(3)) / halfSpan;
        eob   = (obj.Mdl.Caero(i + 1).startZ - ref_point(3)) / halfSpan;
        cib   = 1 - interp1(eta,chord,eib,'linear','extrap');
        cob   = 1 - interp1(eta,chord,eob,'linear','extrap');
        
        obj.Mdl.Caero(i).csType = 'Rudder1';
        
        if count == 1
            obj.Mdl.Caero(i).csId   = 'rudr1r';
        else
            obj.Mdl.Caero(i).csId   = sprintf('rudr%1.0fr',count);
        end
        
        % Redistribute the number of panels
        nchord = obj.Mdl.Caero(i).nx;
        nconchord = cob*nchord;
        nconchord = round(nconchord);       
        nchord = nchord - nconchord;
        obj.Mdl.Caero(i).nx = nchord;

        obj.Mdl.Caero(i).csData = [1 cob cib nconchord];
    end
end

end
