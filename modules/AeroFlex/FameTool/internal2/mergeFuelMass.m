%% MERGEFUELMASS   Create concentrated masses from fuel data provided by
%                  FAME at structural node points
%
% This function merges all the structural mass data that is provided by
% FAME onto the existing strcutural grid points. Thousands of small
% concentrated masses are used in FAME to define the mass data. We want to
% concentrate these masses on the structural grid points. The user has the
% option to slice the data parallel to the XZ-plane or perpendicular to the
% elastic axis of the wing. A slice is made halfway between two grid
% points. A polygon is then drawn along the slice plane up to the point
% where the plane intersects the leading and trailing edge of the planform.
% All the points inside the polygon are identified and concentrated at the
% grid point within this polygon.
%
%   Example:
%      obj = fame2mat;
%      obj.Opts.Struct.massCutPlane = 'EA'; % Perpendicular to elastic axis
%      obj = run(obj);
%
%   See also |getMassData| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
%
% Edits:
%   10/02/2020 - C.Szczyglowski. Updated conversion of raw string data into
%                numeric. Previous call to 'str2num' was returning empty.

function obj = mergeFuelMass(obj, logfcn)

if nargin < 2
    logfcn = @(s) fprintf('%s\n', s);
end

fmPath = obj.Dirs.fameFuelFolder;

fmPath = fullfile(fmPath);
fmDir  = dir(fullfile(fmPath,'right*'));
fmFile = {fmDir.name};

if isdir(fmPath)
    rootpath = fmPath;
elseif isdir(obj.Dirs.fameFolder)
    rootpath = obj.Dirs.fameFolder;
else
    rootpath = pwd;
end

if isempty(fmFile) && obj.HasFuel
    
    logfcn(' * * WARNING * * - the Fuel Files were not found. ');
    logfcn(' Please tell me where it is. If it doesn''t exist then click cancel.');
    [fmFile,fmPath,massidx] = uigetfile([rootpath filesep '*right*'],'Find the fuel distribution files for the RIGHT HAND WING','MultiSelect', 'on');
        
    if massidx == 0 %Account for user clicking "Cancel"
        fmPath = [];
        fmFile = [];
    end
    if ~iscell(fmFile) %Force a cell output
        fmFile = {fmFile};
    end
end

if isempty(fmFile{1}) %Need to ask isempty on the cell as we have just forced a cell output! - C.Szczyglowski 29/06/2018
    warning('FAME2MAT:FileAccessError','File %s does not exist.',fmFile{:});
else
    
    % Save file names    
    for i = 1:numel(fmFile)
        INPUTFILES.fuelMassFiles{i} = fullfile(fmPath,fmFile{i});
    end
    
    obj.Inp = setFields(obj.Inp,INPUTFILES);
    
    % Read mass data
    logfcn(' Parsing FAME fuel mass distribution files...');
    
    % Create conm2 objects with data just retrieved
    logfcn(' Merging FUEL mass data...');
    
    % idx     = find([obj.Mdl.Grid.id] < obj.Opts.Struct.HTPStartID);
    % gridPts = obj.Mdl.Grid(idx);
    
    gridPts = findobj(obj.Mdl.Grid,'part','StbdWing','-and','type','Beam');
    
    minY = min(struct2mat(gridPts,'coord',2)) - 1;
    maxY = max(struct2mat(gridPts,'coord',2)) + 1;
    
    gridVec = cell2mat({gridPts.coord}');
    incrY = unique([minY;gridVec(:,2);maxY]);
    
    % Define points halfway between grid points
    r = numel(incrY) - 1;
    gridHalf = zeros(r,3);
    
    for i = 1:r
        yHalf = (incrY(i) + incrY(i + 1)) / 2;
        xHalf = interp1(gridVec(:,2),gridVec(:,1),yHalf,'linear','extrap');
        zHalf = interp1(gridVec(:,2),gridVec(:,3),yHalf,'linear','extrap');
        gridHalf(i,:) = [xHalf,yHalf,zHalf];
    end
    
    % Define outline of planform for when cut plane is perpendicular with
    % elastic axis
    bounds  = [min(gridHalf(:,2)),max(gridHalf(:,2))];
    outline = obj.Mdl.Caero(1).outline;
    outline = defoutline(outline,bounds);
    outline = cat(1,outline,outline(1,:));
    
    Colours = ['m','g','b','r','c','y'];
    
%     h = gcf;
    
    % Extract mass data for the different fuel cases
    for k = 1:numel(fmFile)
        
        fmFileName = fullfile(fmPath,fmFile{k});
        fid = fopen(fmFileName, 'r');
        wbd = textscan(fid, '%s', 'Delimiter', '\n');
        wbd = wbd{1};
        fclose(fid);
        
        %Remove comments
        idx = startsWith(wbd, {'c', '!'}, 'IgnoreCase',true);
        wbd(idx, :) = [];

        %Remove labels        
        wbd = char(wbd);
        ind = strfind(wbd(1, :), '!');
        if ~isempty(ind)
            wbd = wbd(:, 1 : ind - 1);
        end
        
        data = str2num(wbd);
        if isempty(data)
            %EDIT - ''str2num'' is not working as of 02 Feb 2020. Try
            %different approach instead...

            %Check format - expecting 15 columns of data before the
            %comment
            temp = strsplit(wbd(1, 1 : end - 1), ' ');
            temp(cellfun(@isempty, temp)) = [];
            assert(numel(temp) == 15, ['Exepected the .205 file to ' , ...
                'have 15 columns of data for each mass entry. Check ', ...
                'file and update code.'])
            fmt  = repmat('%f ', [1, 15]);
            data = arrayfun(@(ii) sscanf(wbd(ii, :), fmt), 1 : size(wbd, 1), 'Unif', false);
            data = horzcat(data{:})';
        end
        
        wbdCoord   = data(:, 3 : 5);
        wbdMass    = data(:, 6);
        wbdInertia = data(:, 7 : 12);
        
%         %wbdCoord(:,2) = wbdCoord(:,2);
%         h_leg = findobj(h,'Type','Legend');
%         h_line = findobj(h,'Type','Line');
%         hold on;
%         h_new = plot3(wbdCoord(:,1),wbdCoord(:,2),wbdCoord(:,3),'o','MarkerEdgeColor',Colours(k));axis equal;
%         h_leg.String{end} = fmFile(k).name;
        %legend([h_line;h_new],[h_leg.String,fmFile(k).name]);
        
        %Sort in spanwise direction
        [~,idx]    = sort(wbdCoord(:,2),'ascend');
        wbdCoord   = wbdCoord(idx,:);
        wbdMass    = wbdMass(idx);
        wbdInertia = wbdInertia(idx, :);
        
        j = 1;
         
        for i = 1:numel(gridPts)
            
            % Cut plane can either be spanwise parallel to XZ plane or
            % perpendicular to ealstic axis. A box is drawn halfway between grid
            % points and then the points inside the box are selected.
            if strcmp(obj.Opts.Struct.massCutPlane,'XZ')
                uplim      = gridHalf(i + 1,2);
                lowlim     = gridHalf(i,2);
                idxSmaller = find(wbdCoord(:,2) <=  uplim);
                idxLarger  = find(wbdCoord(:,2) > lowlim);
                
                idx = intersect(idxSmaller,idxLarger);
            else
                [xv,yv] = defpoly(i,gridHalf,gridPts,outline,bounds);
                [idxin,idxon] = inpolygon(wbdCoord(:,1),wbdCoord(:,2),xv,yv);
                
                idxin = find(idxin);
                idxon = find(idxon);
                
                idx = cat(1,idxin,idxon);
            end
            
            Mass = sum(wbdMass(idx));
            
%               if Mass > 0
%                figure;
%                plot(xv,yv);
%                hold on;
%                plot(outline(:,1),outline(:,2));
%                plot(wbdCoord(:,1),wbdCoord(:,2),'b + ');
%                plot(wbdCoord(idx,1),wbdCoord(idx,2),'r + ');
%                axis equal
%                pause
%               end
            
            if Mass > 0
                
                % Calculate Inertia around grid point of beam model
                dX = wbdCoord(idx,1) - ones(numel(idx),1) .* gridPts(i).coord(1);
                dY = wbdCoord(idx,2) - ones(numel(idx),1) .* gridPts(i).coord(2);
                dZ = wbdCoord(idx,3) - ones(numel(idx),1) .* gridPts(i).coord(3);
                
                Ixx = wbdInertia(idx,1);
                Iyy = wbdInertia(idx,2);
                Izz = wbdInertia(idx,3);
                Ixy = wbdInertia(idx,4);
                Ixz = wbdInertia(idx,6);
                Iyz = wbdInertia(idx,5);
                
                Ixx = sum((dY.^2 + dZ.^2) .* wbdMass(idx) + Ixx);
                Iyy = sum((dX.^2 + dZ.^2) .* wbdMass(idx) + Iyy);
                Izz = sum((dX.^2 + dY.^2) .* wbdMass(idx) + Izz);
                Ixy = sum(dX .* dY .* wbdMass(idx) + Ixy);
                Ixz = sum(dX .* dZ .* wbdMass(idx) + Ixz);
                Iyz = sum(dY .* dZ .* wbdMass(idx) + Iyz);
                
                Inertia = [ Ixx, Ixy, Ixz;
                    Ixy, Iyy, Iyz;
                    Ixz, Iyz, Izz];
                
                if j == 1
                    Conm2Obj       = Conm2;
                else
                    Conm2Obj(j)    = Conm2;
                end
                
                Conm2Obj(j).id     = gridPts(i).id + 500000;%100000 * k;
                Conm2Obj(j).grid   = gridPts(i).id;
                Conm2Obj(j).cid    = 0;
                Conm2Obj(j).m      = Mass;
                Conm2Obj(j).offset = [0,0,0];
                Conm2Obj(j).i      = Inertia;
                Conm2Obj(j).file   = fmFile{k};
                Conm2Obj(j).type   = 'Fuel';
                Conm2Obj(j).part   = 'StbdWing';
                
                j = j + 1;
                
                wbdCoord(idx,:) = [];
                wbdMass(idx,:) = [];
            end
            
        end
        
        obj.Mdl.Conm2 = cat(2,obj.Mdl.Conm2,Conm2Obj);
        
        
    end
    
end

% if ~exist([obj.Dirs.writeFolder filesep 'ModelFigures'],'dir')
%     mkdir([obj.Dirs.writeFolder filesep 'ModelFigures'])
% end

% saveas(h,[obj.Dirs.writeFolder filesep 'ModelFigures' filesep 'WingMasses.fig']);


end







