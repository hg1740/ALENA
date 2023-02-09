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
function obj = controlSurfaceGenerator(obj,part)

% Find the corresponding CAERO cards for that part
CAEROparts  = {obj.Mdl.Caero.part};
caeropartid = ~cellfun(@isempty,regexpi(CAEROparts,part));
CAERO       = obj.Mdl.Caero(caeropartid);

AeroIdx     = find(caeropartid);

% Find the corresponding GRID cards for that part
GRIDparts   = {obj.Mdl.Grid.part};
gridpartid  = ~cellfun(@isempty,regexpi(GRIDparts,part));
PartNodes   = obj.Mdl.Grid(gridpartid);

% CAERO     = findobj(obj.Mdl.Caero,'part','StbdWing');
% WingNodes = findobj(obj.Mdl.Grid,'part','StbdWing');

StartNumber = CAERO(end).id + (CAERO(end).nx*CAERO(end).ny);
XtraPanel = 0;

Fame     = obj.Fame;

% What is the half span of the lifting surface
if strcmpi(part,'vtp')
    maxeta = max(struct2mat(PartNodes,'coord',3));
    mineta = min(struct2mat(PartNodes,'coord',3));
    halfSpan = maxeta - mineta;
    ref_point = [obj.Mdl.Caero(AeroIdx(1)).startX,...
        obj.Mdl.Caero(AeroIdx(1)).startY,...
        obj.Mdl.Caero(AeroIdx(1)).startZ];

else
    halfSpan = max(struct2mat(PartNodes,'coord',2));
end

% List the available control surfaces
CsNames = fieldnames(obj.Fame.Geometry.(part).ControlSurfaces);

IgnoreCsNames ={'slat','flap'};

% Let's identify which control surface indices we can create in the beam
% model
CsIdx   = (1:numel(CsNames))';

for i = 1:numel(IgnoreCsNames)
    CsIdx = intersect(CsIdx,find(cellfun(@isempty,regexpi(CsNames,IgnoreCsNames{i}))));
end

% Iterate through the the valid control surfaces
for k = 1:numel(CsIdx)
    
    CsPart = CsNames{CsIdx(k)};
    
    % Iterate through each sub control surface
    for j = 1:length(Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_beg_le)
        
        if strcmpi(part,'vtp')
            [~,InboardCutLineIndex]  = min(abs(struct2mat(CAERO,'startZ') - (halfSpan * Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_beg_le(j) + mineta)));
            [~,OutboardCutLineIndex] = min(abs(struct2mat(CAERO,'startZ') + struct2mat(CAERO,'b') - (halfSpan * Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_end_le(j) + mineta)));
        else
            [~,InboardCutLineIndex]  = min(abs(struct2mat(CAERO,'startY') - (halfSpan * Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_beg_le(j))));
            [~,OutboardCutLineIndex] = min(abs(struct2mat(CAERO,'startY') + struct2mat(CAERO,'b') - (halfSpan * Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_end_le(j))));
        end
        
        count = 0;
        aerocount  = 0;
        
        for i = AeroIdx
            aerocount = aerocount + 1;
            if aerocount >= InboardCutLineIndex && aerocount <= OutboardCutLineIndex
                count = count + 1;
                eta   = [Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_beg_le(j),Fame.Geometry.(part).ControlSurfaces.(CsPart).eta_end_le(j)];
                chord = [Fame.Geometry.(part).ControlSurfaces.(CsPart).chrd_beg_le(j),Fame.Geometry.(part).ControlSurfaces.(CsPart).chrd_end_le(j)];
                
                if strcmpi(part,'vtp')
                    eib   = (obj.Mdl.Caero(i).startZ - ref_point(3)) / halfSpan;
                    eob   = (obj.Mdl.Caero(i).startZ + obj.Mdl.Caero(i).b - ref_point(3)) / halfSpan;
                else
                    eib   = obj.Mdl.Caero(i).startY / halfSpan;
                    eob   = (obj.Mdl.Caero(i).startY + obj.Mdl.Caero(i).b) / halfSpan;
                    
                end
                
                cib   = 1 - interp1(eta,chord,eib,'linear','extrap');
                cob   = 1 - interp1(eta,chord,eob,'linear','extrap');
                
                obj.Mdl.Caero(i).csType = [CsNames{CsIdx(k)} num2str(j)];
                if count == 1
                    obj.Mdl.Caero(i).csId   = [CsPart(1:3) num2str(j)];
                else
                    obj.Mdl.Caero(i).csId   = [CsPart(1:3) num2str(j) '_' num2str(count-1)];
                end
                % Redistribute the number of panels
                nchord = obj.Mdl.Caero(i).nx;
                nconchord = cob*nchord;
                nconchord = round(nconchord);
                nchord = nchord - nconchord;
                obj.Mdl.Caero(i).nx = nchord;
                
                if nconchord == 0
                    nconchord = 1;
                end
                
                obj.Mdl.Caero(i).csData = [1 cib cob nconchord];
                obj.Mdl.Caero(i).fc     = [cib cob];
                
                obj.Mdl.Caero(i).fnx = nconchord;
                obj.Mdl.Caero(i).flapped = 1;
                
            end
        end
    end
end
% 
% for j = 1:length(Fame.Geometry.(part).ControlSurfaces.aileron.eta_beg_le)
%     
%     [~,aileronInboardCutLineIndex] = min(abs(struct2mat(CAERO,'startY') - (halfSpan * Fame.Geometry.Wing.ControlSurfaces.aileron.eta_beg_le(j))));
%     [~,aileronOutboardCutLineIndex] = min(abs(struct2mat(CAERO,'startY') + struct2mat(CAERO,'b') - (halfSpan * Fame.Geometry.Wing.ControlSurfaces.aileron.eta_end_le(j))));
%     count = 0;
%     
%     for i = 1:numel(obj.Mdl.Caero)
%         
%         if i >= aileronInboardCutLineIndex && i <= aileronOutboardCutLineIndex
%             count = count + 1;
%             eta   = [Fame.Geometry.Wing.ControlSurfaces.aileron.eta_beg_le(j),Fame.Geometry.Wing.ControlSurfaces.aileron.eta_end_le(j)];
%             chord_le = [Fame.Geometry.Wing.ControlSurfaces.aileron.chrd_beg_le(j),Fame.Geometry.Wing.ControlSurfaces.aileron.chrd_end_le(j)];
%             eib   = obj.Mdl.Caero(i).startY / halfSpan;
%             eob   = obj.Mdl.Caero(i + 1).startY / halfSpan;
%             cib   = 1 - interp1(eta,chord_le,eib,'linear','extrap');
%             cob_le   = 1 - interp1(eta,chord_le,eob,'linear','extrap');
%             
%             obj.Mdl.Caero(i).csType = ['Aileron' num2str(j)];
%             
%             if count == 1
%                 obj.Mdl.Caero(i).csId   = ['ail' num2str(1 + 100*j) 'r'];
%             else
%                 obj.Mdl.Caero(i).csId   = sprintf('ail%1.0fr',count + 100*j);
%             end
%             
%             % Redistribute the number of panels
%             nchord = obj.Mdl.Caero(i).nx;
%             nconchord = cob_le*nchord;
%             nconchord = round(nconchord);
%             nchord = nchord - nconchord;
%             obj.Mdl.Caero(i).nx = nchord;
%             
%             obj.Mdl.Caero(i).csData = [1 cob_le cib nconchord];
%         end
%     end
% end

% % TODO: Additional code for the spoiler
% for j = 1:length(Fame.Geometry.Wing.ControlSurfaces.spoiler.eta_beg_le)
%     [~,spoilerInboardCutLineIndex] = min(abs(struct2mat(CAERO,'startY') - (halfSpan * Fame.Geometry.Wing.ControlSurfaces.spoiler.eta_beg_le(j))));
%     
%     [~,spoilerOutboardCutLineIndex] = min(abs(struct2mat(CAERO,'startY') + struct2mat(CAERO,'b') - (halfSpan * Fame.Geometry.Wing.ControlSurfaces.spoiler.eta_end_le(j))));
%     count = 0;
%     
%     for i = 1:numel(obj.Mdl.Caero)
%         
%         if i >= spoilerInboardCutLineIndex && i <= spoilerOutboardCutLineIndex
%             count = count + 1;
%             eta   = [Fame.Geometry.Wing.ControlSurfaces.spoiler.eta_beg_le(j),Fame.Geometry.Wing.ControlSurfaces.spoiler.eta_end_le(j)];
%             
%             chord_le = [Fame.Geometry.Wing.ControlSurfaces.spoiler.chrd_beg_le(j),Fame.Geometry.Wing.ControlSurfaces.spoiler.chrd_end_le(j)];
%             chord_te = [Fame.Geometry.Wing.ControlSurfaces.spoiler.chrd_beg_te(j),Fame.Geometry.Wing.ControlSurfaces.spoiler.chrd_end_te(j)];
%             
%             eib   = obj.Mdl.Caero(i).startY / halfSpan;
%             eob   = obj.Mdl.Caero(i + 1).startY / halfSpan;
%             
%             cib_le   = 1 - interp1(eta,chord_le,eib,'linear','extrap');
%             cob_le   = 1 - interp1(eta,chord_le,eob,'linear','extrap');
%             
%             cib_te   = 1 - interp1(eta,chord_te,eib,'linear','extrap');
%             cob_te   = 1 - interp1(eta,chord_te,eob,'linear','extrap');
%             
%             % Determine how much of the CAERO Card needs to be kept
%             WingFac_ib_te = 1 - cib_te;
%             WingFac_ob_te = 1 - cob_te;
%             
%             WingFac_ib_le = 1 - cib_le;
%             WingFac_ob_le = 1 - cob_le;
%             
%             %% Alter the existing CAERO card
%             
%             % Fetch the outline.
%             
%             OldVert_1 = obj.Mdl.Caero(i).outline(1,:);
%             OldVert_2 = obj.Mdl.Caero(i).outline(2,:);
%             OldVert_3 = obj.Mdl.Caero(i).outline(3,:);
%             OldVert_4 = obj.Mdl.Caero(i).outline(4,:);
%             
% 
%             
%             %%
%             obj.Mdl.Caero(i).csType = ['Spoiler' num2str(j)];
%             
%             if count == 1
%                 obj.Mdl.Caero(i).csId   = ['Spl' num2str(1 + 100*j) 'r'];
%             else
%                 obj.Mdl.Caero(i).csId   = sprintf('Spl%1.0fr',count + 100*j);
%             end
%             
%             % Redistribute the number of panels
%             nchord = obj.Mdl.Caero(i).nx;
%             nconchord = cob_le*nchord;
%             nconchord = round(nconchord);
%             nchord = nchord - nconchord;
%             obj.Mdl.Caero(i).nx = nchord;
%             
%             obj.Mdl.Caero(i).csData = [1 cib_le cob_le nconchord];
% 
%             
%             % The CAERO card needs to be altered to account for the control
%             % surface and then an additional CAERO card needs to be added
%             % to make up for the lost trailing edge.
%             
%             if ~obj.Opts.Aero.extendspoiler
%                 
%                 % Truncate the existing chord
%                 old_c = obj.Mdl.Caero(i).c;
%                 obj.Mdl.Caero(i).c = WingFac_ib_te*obj.Mdl.Caero(i).c;
%                 
%                 % Correctly change the taper ratio ...  just in case theres
%                 % some additional sweep to the control surface
%                 old_t = obj.Mdl.Caero(i).t;
%                 obj.Mdl.Caero(i).t = (WingFac_ib_te/WingFac_ob_te)*obj.Mdl.Caero(i).t;
%                 
%                 % Recalculate the sweep angle as the chord has been truncated!
%                 obj.Mdl.Caero(i).sw = 180*atan(((OldVert_2(1)+0.25*obj.Mdl.Caero(i).t*obj.Mdl.Caero(i).c) - (OldVert_1(1)+0.25*obj.Mdl.Caero(i).c))/...
%                     obj.Mdl.Caero(i).b)/pi;
%                 
%                 obj.Mdl.Caero(i).csData = [1 (WingFac_ib_te-WingFac_ib_le)/WingFac_ib_te (WingFac_ob_te-WingFac_ob_le)/WingFac_ob_te nconchord];
%                 
%                 % Reduce the chord of the existing section and apply the
%                 % control surface to that. Then add another CAERO card with the
%                 % left over chord
%                 %% Add another CAERO Card
%                 CaeroObj        = Caero;
%                 CaeroObj.id     = StartNumber + XtraPanel;
%                 CaeroObj.startX = obj.Mdl.Caero(i).startX + WingFac_ib_te*old_c;
%                 CaeroObj.startY = obj.Mdl.Caero(i).startY;
%                 CaeroObj.startZ = obj.Mdl.Caero(i).startZ;
%                 CaeroObj.twist1 = obj.Mdl.Caero(i).twist1;
%                 CaeroObj.twist2 = obj.Mdl.Caero(i).twist2;
%                 CaeroObj.c      = (1-WingFac_ib_te)*old_c;
%                 CaeroObj.t      = ((1-WingFac_ib_te)/(1-WingFac_ob_te))*old_t;
%                 CaeroObj.b      = obj.Mdl.Caero(i).b;
%                 CaeroObj.sw     = 180*atan(((obj.Mdl.Caero(i).outline(3,1)+0.25*CaeroObj.t*CaeroObj.c) - (obj.Mdl.Caero(i).outline(4,1)+0.25*CaeroObj.c))/...
%                     CaeroObj.b)/pi;
%                 CaeroObj.rootAirfoil  = obj.Mdl.Caero(i).rootAirfoil;
%                 CaeroObj.tipAirfoil   = obj.Mdl.Caero(i).tipAirfoil;
%                 CaeroObj.dih    = obj.Mdl.Caero(i).dih;
%                 CaeroObj.meshType     = obj.Mdl.Caero(i).meshType;
%                 CaeroObj.cid    = obj.Mdl.Caero(i).cid;
%                 CaeroObj.ny     = obj.Mdl.Caero(i).ny;
%                 CaeroObj.nx     = 1;
%                 CaeroObj.part   = 'StbdWing';
%                 
%                 obj.Mdl.Caero   = cat(2,obj.Mdl.Caero,CaeroObj);
%                 
%                 XtraPanel = XtraPanel + CaeroObj.ny;
%             end
%         end
%     end
% end
end