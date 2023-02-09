%% MERGESTRUCTMASS   Create concentrated masses from structural data from
%                    FAME at structural node points
%
% This function merges all the fuel mass data onto the existing grid
% points. Thousands of small concentrated masses are used in FAME to define
% the fuel data. We want to concentrate these masses on the structural grid
% points. The user has the option to slice the data parallel to the
% XZ-plane or perpendicular to the elastic axis of the wing. A slice is
% made halfway between two grid points. A polygon is then drawn along the
% slice plane up to the point where the plane intersects the leading and
% trailing edge of the planform. All the points inside the polygon are
% identified and concentrated at the grid point within this polygon.
%
%   Example:
%      obj = fame2mat;
%      obj.Opts.Struct.massCutPlane = 'EA'; % Perpendicular to elastic axis
%      obj = run(obj);
%
%   See also |getMassData| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.

function obj = mergeStructMass(obj, logfcn)

if nargin < 2
    logfcn = @(s) fprintf('%s\n', s);
end

fmPath = obj.Dirs.fameFuelFolder;

fmPath = fullfile(fmPath);
fmDir  = dir(fullfile(fmPath,'*.205'));
fmFile = {fmDir.name};

if isdir(fmPath)
    rootpath = fmPath;
elseif isdir(obj.Dirs.fameFolder)
    rootpath = obj.Dirs.fameFolder;
else
    rootpath = pwd;
end

if isempty(fmFile) && obj.Has205
    
    logfcn(' * * WARNING * * - the ''massdist_items__sorted.205'' cannot be found.');
    logfcn(' Please tell me where it is. If it doesn'' exist then click cancel');
    
    [fmFile,fmPath,filterindex] = uigetfile([rootpath filesep '*.205'],'Pick the .205 mass distribution file');
    
    %[fmFile,fmPath,massidx] = uigetfile([rootpath filesep '*right*'],'Find the fuel distribution files for the RIGHT HAND WING','MultiSelect', 'on');
    
    if filterindex == 0
        fmPath = [];
        fmFile = [];
    end
    
    if ~iscell(fmFile)
        fmFile = {fmFile};
    end
end

if isempty(fmFile{1})
    logfcn(' * * WARNING * * - the ''*.205'' file has not been selected.');
else
    
    % Check to see which file is in the directory. The *.205 file may be sorted
    % or unsorted.
    sortidx    = find(~cellfun(@isempty,regexp(fmFile,'sorted')));
    
    if sortidx
        sortedFile = 1;
        wbdFile    = fullfile(fmPath,fmFile{sortidx});
    else
        sortedFile = 0;
        wbdFile    = fullfile(fmPath,fmFile{1});
    end

    % Save file names
    INPUTFILES.wbdFiles = wbdFile;
    obj.Inp = setFields(obj.Inp,INPUTFILES);
    
    % Read mass data
    logfcn(' Parsing FAME mass distribution file...');
    
    if sortedFile == 1
        wbd = readWbd2Sorted(wbdFile);
    else
        wbd = readWbd2(wbdFile);
    end
    
    % Correct mass data if needed. Some points default to reference position if
    % no CG position calculated or given by weight function.
    obj.Opts.Struct.massCorrection = 0;
    if obj.Opts.Struct.massCorrection == 1
        wbdCorrFile = obj.Inp.wbdCorrFiles;
        wbd         = correctWbd2(wbd,wbdCorrFile);
    end
    
    if sortedFile == 1
        % Identify the structural and secondary masses
        PartNames = {wbd.Name};
        
        %% Structural Mass
        
        % Skin
        tags = {'Lower_Skin','Upper_Skin','Spanwise_Doubler','Chordwise_Doubler','Rodfittings','Man_Holes'};
        idx_skin = [];
        for i = 1:numel(tags)
            idx_skin = [idx_skin,find(~cellfun(@isempty,regexp(PartNames,tags{i})))];
        end
        idx_skin = unique(idx_skin);
        
        % Stringers
        tags = {'Lower_Stringers','Upper_Stringers'};
        idx_string = [];
        for i = 1:numel(tags)
            idx_string = [idx_string,find(~cellfun(@isempty,regexp(PartNames,tags{i})))];
        end
        idx_string = unique(idx_string);
        
        % Spars
        tags = {'Spar_Webs','Flange','Columns','Bolting_Skin_to_Spar_','GE:slat_track_cans'};
        idx_spar = [];
        for i = 1:numel(tags)
            idx_spar = [idx_spar,find(~cellfun(@isempty,regexp(PartNames,tags{i})))];
        end
        idx_spar = unique(idx_spar);
        
        % Ribs
        tags = {'Rib_No','Doubler_Front_Spar','Doubler_Rear_Spar','Cruciform_Part_Upper_Shell',...
            'Triform_Part_Lower_Shell','TRIFORM_v02','CRUCIFORM_v02','DOUBLER_FSP','DOUBLER_RSP'};
        idx_rib = [];
        for i = 1:numel(tags)
            idx_rib = [idx_rib,find(~cellfun(@isempty,regexp(PartNames,tags{i})))];
        end
        idx_rib = unique(idx_rib);
        
        idx_struc = unique([idx_skin,idx_string,idx_spar,idx_rib]);
        
        idx_secon = 1:numel(wbd);
        idx_secon(idx_struc) = [];
        
    end
    
    % Build coordinate and mass vectors of all the points
    logfcn(' Merging structural mass data...');
    
    if sortedFile == 1
        
        wbdCoord = [];
        wbdMass  = [];
        
        wbdCoord_sec = [];
        wbdMass_sec  = [];
        
        for i = idx_struc
            wbdCoord = cat(1,wbdCoord,[wbd(i).Pts.x,wbd(i).Pts.y,wbd(i).Pts.z]);
            wbdMass  = cat(1,wbdMass,wbd(i).Pts.Mass);
        end
        
        for i = idx_secon
            wbdCoord_sec = cat(1,wbdCoord_sec,[wbd(i).Pts.x,wbd(i).Pts.y,wbd(i).Pts.z]);
            wbdMass_sec  = cat(1,wbdMass_sec,wbd(i).Pts.Mass);
        end
        
        StructuralMass = wbd(idx_struc);
        SecondaryMass  = wbd(idx_secon);
    else
        wbdCoord = [wbd.x,wbd.y,wbd.z];
        wbdMass  = wbd.Mass;
        
        StructuralMass = wbd;
        SecondaryMass  = [];
    end
    
    [~,idx]  = sort(wbdCoord(:,2),'ascend');
    wbdCoord = wbdCoord(idx,:);
    wbdMass  = wbdMass(idx);
    
    % Ignore any null mass entries
    idx = find(wbdMass == 0);
    if ~isempty(idx)
        wbdMass(idx) = [];
        wbdCoord(idx,:) = [];
    end
    
    % Check to see if some cg's are incorrect
    idx = find(wbdCoord(:,1) == 0);
    
    if ~isempty(idx)
        logfcn(sprintf([' * * WARNING * * - Some of the mass entries in the *.205 have [0,0,0] ', ...
            'location. These entries have been removed and comprise of %8.2f Kg'],sum(wbdMass(idx))));
        wbdMass(idx) = [];
        wbdCoord(idx,:) = [];
    end
    
    if sortedFile == 1
        
        [~,idx]  = sort(wbdCoord_sec(:,2),'ascend');
        wbdCoord_sec = wbdCoord_sec(idx,:);
        wbdMass_sec  = wbdMass_sec(idx);
        
        % Ignore any null mass entries
        idx = find(wbdMass_sec == 0);
        if ~isempty(idx)
            wbdMass_sec(idx) = [];
            wbdCoord_sec(idx,:) = [];
        end
        
        % Check to see if some cg's are incorrect
        idx = find(wbdCoord_sec(:,1) == 0);
        
        if ~isempty(idx)
            logfcn(sprintf([' * * WARNING * * - Some of the mass entries in the *.205 have [0,0,0] ', ...
                'location. These entries have been removed and comprise of %8.2f Kg'],sum(wbdMass(idx))));
            wbdMass_sec(idx) = [];
            wbdCoord_sec(idx,:) = [];
        end
    end
    
%     % Plot the masses
%     figure;h_struc = plot3(wbdCoord(:,1),wbdCoord(:,2),wbdCoord(:,3),'MarkerFaceColor','r','Marker','o','MarkerEdgeColor','k','LineStyle','none','MarkerSize',6);axis equal;
%     if sortedFile == 1
%         hold on;h_sec = plot3(wbdCoord_sec(:,1),wbdCoord_sec(:,2),wbdCoord_sec(:,3),'MarkerFaceColor','b','Marker','o','MarkerEdgeColor','k','LineStyle','none','MarkerSize',6);axis equal;
%     end
%     xlabel('x (m)');
%     ylabel('y (m)');
%     zlabel('z (m)');
%     
%     legend([h_struc,h_sec],'Structural Mass','Secondary Mass','Location','Best');
    
    % Reconsolidate the mass and coordinate points. Deal with this later
    
    wbdMass = [wbdMass;wbdMass_sec];
    wbdCoord = [wbdCoord;wbdCoord_sec];
    
    % Find structural reference points
    
    gridPts = findobj(obj.Mdl.Grid,'part','StbdWing','-and','type','Beam');
    
    % idx     = find([obj.Mdl.Grid.id] < obj.Opts.Struct.HTPStartID);
    % gridPts = obj.Mdl.Grid(idx);
    
    minY = min(struct2mat(gridPts,'coord',2)) - 1;
    maxY = max(struct2mat(gridPts,'coord',2)) + 1;
    
    gridVec = cell2mat({gridPts.coord}');
    incrY   = unique([minY;gridVec(:,2);maxY]);
    
    % Define points halfway between grid points
    r = numel(incrY) - 1;
    gridHalf = zeros(r,3);
    
    for i = 1:r
        yHalf = (incrY(i) + incrY(i + 1)) / 2;
        xHalf = interp1(gridVec(:,2),gridVec(:,1),yHalf,'linear','extrap');
        zHalf = interp1(gridVec(:,2),gridVec(:,3),yHalf,'linear','extrap');
        gridHalf(i,:) = [xHalf,yHalf,zHalf];
    end
    
    bounds  = [min(gridHalf(:,2)),max(gridHalf(:,2))];
    
    % Define outline of planform for when cut plane is perpindicular with
    % elastic axis, and remove points outside of planform
    outline = obj.Fame.Geometry.Wing.Outline;
    outline = defoutline(outline,bounds);
    outline = cat(1,outline,outline(1,:));
    
    [idxin,idxon] = inpolygon(wbdCoord(:,1),wbdCoord(:,2),outline(:,1),outline(:,2));
    
    idxin = find(idxin);
    idxon = find(idxon);
    
    idx      = [idxin;idxon];
    wbdCoord = wbdCoord(idx,:);
    wbdMass  = wbdMass(idx);
    
    j = 1;
    SumMass = 0;
    SumMassCG = 0;
    for i = 1:numel(gridPts)
        
        % Cut plane can either be spanwise parallel to XZ plane or
        % perpendicular to ealstic axis. a box is drawn halfway between grid
        % points and then the points inside the box are selected.
        if strcmp(obj.Opts.Struct.massCutPlane,'XZ')
            uplim      = gridHalf(i + 1,2);
            lowlim     = gridHalf(i,2);
            idxSmaller = find(wbdCoord(:,2) <=  uplim);
            idxLarger  = find(wbdCoord(:,2) > lowlim);
            
            idx = intersect(idxSmaller,idxLarger);
        else
            [xv,yv]       = defpoly(i,gridHalf,gridPts,outline,bounds);
            [idxin,idxon] = inpolygon(wbdCoord(:,1),wbdCoord(:,2),xv,yv);
            
            idxin = find(idxin);
            idxon = find(idxon);
            
            idx = cat(1,idxin,idxon);
        end
        
        Mass = sum(wbdMass(idx));
        
        if ~isempty(wbdMass) && Mass > 0
            
            % Calculate the centre of mass
            Cgx = (wbdCoord(idx,1)'* wbdMass(idx))/sum(wbdMass(idx));
            Cgy = (wbdCoord(idx,2)'* wbdMass(idx))/sum(wbdMass(idx));
            Cgz = (wbdCoord(idx,3)'* wbdMass(idx))/sum(wbdMass(idx));
            
            CG = [Cgx,Cgy,Cgz];
            
            SumMassCG = SumMassCG + Mass*CG;
            SumMass = SumMass + Mass;
            
            % Calculate Inertia around centre of mass
            dX = wbdCoord(idx,1) - ones(numel(idx),1) .*  Cgx;
            dY = wbdCoord(idx,2) - ones(numel(idx),1) .*  Cgy;
            dZ = wbdCoord(idx,3) - ones(numel(idx),1) .*  Cgz;
            
            Ixx = sum((dY.^2 + dZ.^2) .*  wbdMass(idx));
            Iyy = sum((dX.^2 + dZ.^2) .*  wbdMass(idx));
            Izz = sum((dX.^2 + dY.^2) .*  wbdMass(idx));
            Ixy = sum(dX .*  dY .*  wbdMass(idx));
            Ixz = sum(dX .*  dZ .*  wbdMass(idx));
            Iyz = sum(dY .*  dZ .*  wbdMass(idx));
            
            Inertia = [ Ixx, Ixy, Ixz;
                Ixy, Iyy, Iyz;
                Ixz, Iyz, Izz];
            
            offset = CG - gridPts(i).coord(:)';
            
            Conm2Obj(j)        = Conm2;
            Conm2Obj(j).id     = gridPts(i).id;
            Conm2Obj(j).grid   = gridPts(i).id;
            Conm2Obj(j).cid    = 0;
            Conm2Obj(j).m      = Mass;
            Conm2Obj(j).offset = offset;
            Conm2Obj(j).i      = Inertia;
            Conm2Obj(j).type   = 'Structural';
            Conm2Obj(j).part   = 'StbdWing';
            
            j = j + 1;
            
            wbdCoord(idx,:) = [];
            wbdMass(idx,:)  = [];
        end
        
    end
    
    
    obj.Mdl.Conm2 = Conm2Obj;
    
    Inertia = cell2mat({Conm2Obj.i}');
    Mass    = cell2mat({Conm2Obj.m}');
    offset  = cell2mat({Conm2Obj.offset}');
    
    obj.Fame.BeamModel.Wing.Mass     = Mass;
    obj.Fame.BeamModel.Wing.Mass_I11 = Inertia(1:3:end,1);
    obj.Fame.BeamModel.Wing.Mass_I21 = Inertia(1:3:end,2);
    obj.Fame.BeamModel.Wing.Mass_I22 = Inertia(2:3:end,2);
    obj.Fame.BeamModel.Wing.Mass_I31 = Inertia(1:3:end,3);
    obj.Fame.BeamModel.Wing.Mass_I32 = Inertia(2:3:end,3);
    obj.Fame.BeamModel.Wing.Mass_I33 = Inertia(3:3:end,3);
    obj.Fame.BeamModel.Wing.Mass_offset = offset;
    
    obj.Fame.Geometry.Weights.Wing.StructuralMass = StructuralMass;
    obj.Fame.Geometry.Weights.Wing.SecondaryMass  = SecondaryMass;
end
end
