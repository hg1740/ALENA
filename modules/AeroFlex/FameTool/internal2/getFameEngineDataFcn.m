%% GETFAMEENGINEDATAFCN Obtain engine data from FAME
%

%   Copyright 2016 University of Bristol
%   Private function.
function obj = getFameEngineDataFcn(obj)

folderName = obj.Dirs.fameFolder;
NODE       = obj.Mdl.Grid;
CONM2      = obj.Mdl.Conm2;
Fame       = obj.Fame;

dotstr   = findstr('.fame',folderName);
slashstr = findstr(filesep,folderName);

pathName = '__mirror';
fileName = dir(fullfile(folderName,pathName,'*.fm4'));
INPUTFILES.fameEngineFiles{1} = fullfile(folderName,pathName,fileName(1).name);
neta = 0;

idxStruct = find([NODE.id] < 130000)';

startIdOffsetNumber = 200000;

if ~exist(INPUTFILES.fameEngineFiles{1},'file')
    error('Unable to find file %s.\n', INPUTFILES.fameEngineFiles{1});
else
    fp = fopen(INPUTFILES.fameEngineFiles{1}, 'r');
    skipLine = false;
    
    while ~feof(fp)
        
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
        stringdata = strfind(tline,'< ENGINE_WING_ >');
        if ~isempty(stringdata)
            neta  = neta + 1;
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            
            eta   = sscanf(tline,' < ETA >                     [%f') * max(struct2mat(NODE(idxStruct),'coord',2));
            
            fgetl(fp);
            tline = fgetl(fp);
            
            mass  = sscanf(tline,' < ENGINE_MASS >             [%f');
            
            fgetl(fp);
            tline = fgetl(fp);
            
            xLoc  = sscanf(tline,' < POD_PYLON_CG_X_LOC >      [%f') / 1000;
            
            tline = fgetl(fp);
            
            zLoc  = sscanf(tline,' < POD_PYLON_CG_Z_LOC >      [%f') / 1000;
            [~,iEngine] = min(abs(struct2mat(NODE(idxStruct),'coord',2) - eta));
            yLoc  = eta - NODE(iEngine).coord(2);
            
            % Create a new node to coincide with the Pylon - Engine c.g.
            GridObj       = Grid;
            GridObj.id    = 310001;
            GridObj.cs    = 0;
            GridObj.coord = NODE(iEngine).coord + [ - 2.5,yLoc,zLoc];
            GridObj.cd    = 0;
            GridObj.ps    = 0;
            GridObj.seid  = 0;
            NODE(end + 1) = GridObj;
            
            % Place a lumped mass at the c.g. point
            Conm2Obj        = Conm2;
            Conm2Obj.id     = GridObj.id;
            Conm2Obj.grid   = GridObj.id;
            Conm2Obj.cid    = 0;
            Conm2Obj.m      = mass;
            Conm2Obj.offset = [0,0,0];
            Conm2Obj.i      = zeros(3,3);
            Conm2Obj.type   = 'Engine';
            CONM2(end + 1)  = Conm2Obj;
            
            % Create a new RBE2 between the beam and pylon / engine
            Rbe2Obj     = Rbe2;
            Rbe2Obj.id  = GridObj.id;
            Rbe2Obj.gn  = NODE(iEngine).id;
            Rbe2Obj.dof = 123456;
            Rbe2Obj.gm  = GridObj.id;
            
            % Create a force card to mimick thrust
            ThrustObj     = Thrust;
            ThrustObj.lid = 310001;
            ThrustObj.g   = GridObj.id;
            ThrustObj.cid = 0;
            ThrustObj.cx  = Fame.LoadCases(1,8);
            ThrustObj.cy  = 0;
            ThrustObj.cz  = 0;
            ThrustObj.ox  = 0;
            ThrustObj.oy  = 0;
            ThrustObj.oz  = 0;
            
        end
        
    end % end of while
    fclose(fp);
    
end

if neta >0
    obj.Mdl.Grid   = NODE;
    obj.Mdl.Conm2  = CONM2;
    obj.Mdl.Rbe2   = Rbe2Obj;
    obj.Mdl.Thrust = ThrustObj;
    obj.Inp        = setFields(obj.Inp,INPUTFILES);
end
end

