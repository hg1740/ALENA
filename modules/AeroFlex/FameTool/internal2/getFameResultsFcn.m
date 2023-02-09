%% GETFAMERESULTSFCN Obtain Results from FAME
%
% Source FAME file : '1_structure\1_10_wing\l3_fame-w\flexible\loadcase_info.txt'
%                                                              'wng_box_properties.fmwpp01'
%                                                              'wng_addgeodata.fmwpp01'

%   Copyright 2016 University of Bristol
%   Private function.
function [Results,Inputfiles] = getFameResultsFcn(folderName)

% The pathname is essentially hardcoded 
pathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible'];

%% FETCH MORE LOAD CASES
fileName = 'loadcase_info.txt';
Inputfiles.loadCaseFiles{1} = fullfile(folderName,pathName,fileName);

spoilcount    = 0;
ailcount      = 0;
nTrim         = 0;
spoilerfactor = 0.5;

if ~exist(Inputfiles.loadCaseFiles{1},'file')
    error('Unable to find file %s.\n', Inputfiles.loadCaseFiles{1});
else
    fp = fopen(Inputfiles.loadCaseFiles{1}, 'r');
    
    skipLine = false;
    %fprintf(fid,'Reading %s file...\n', Inputfiles.loadCaseFiles{1});
    while ~feof(fp)
        
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
        % Fetch control surface deflections for each load case
        stringData = strfind(tline,'Aileron deflection [degree]:');
        
        if ~isempty(stringData)
            stringData      = textscan(tline,'Aileron deflection [degree]: %f');
            ailcount        = ailcount + 1;
            ailDef(ailcount,:) = stringData{1,1};
        end
        
        stringData = strfind(tline,'Spoiler deflection [degree]:');
        
        if ~isempty(stringData)
            stringData     = textscan(tline,'Spoiler deflection [degree]: %f %f %f %f %f %f');
            spoilcount     = spoilcount + 1;
            for i = 1:length(stringData)
                spoilDef(spoilcount,i) = stringData{1,i};
            end
        end
        
        % Fetch flight conditions
        stringData = strfind(tline,'LC   AC_mass');
        
        if ~isempty(stringData)
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                nTrim = nTrim + 1;
                
                LoadCases(nTrim,1) = stringData{1,1}(1);
                LoadCases(nTrim,2) = stringData{1,2}(1);
                LoadCases(nTrim,3) = stringData{1,3}(1);
                LoadCases(nTrim,4) = stringData{1,4}(1);
                LoadCases(nTrim,5) = stringData{1,5}(1);
                LoadCases(nTrim,6) = stringData{1,6}(1);
                LoadCases(nTrim,7) = stringData{1,7}(1);
                LoadCases(nTrim,8) = stringData{1,8}(1);
                
                tline = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            end
        end
        
        stringData = strfind(tline,'LC   Machzahl');
        
        if ~isempty(stringData)
            count = 0;
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                count = count + 1;
                
                LoadCases(count,9)  = stringData{1,2}(1);
                LoadCases(count,10) = stringData{1,3}(1);
                LoadCases(count,11) = stringData{1,4}(1);
                LoadCases(count,12) = stringData{1,5}(1);
                LoadCases(count,13) = stringData{1,6}(1);
                LoadCases(count,14) = stringData{1,7}(1);
                LoadCases(count,15) = stringData{1,8}(1);
                
                tline = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            end
        end
        
        stringData = strfind(tline,'LC   AC-Lift');
        
        if ~isempty(stringData)
            count = 0;
            fgetl(fp);
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                count = count + 1;
                
                LoadCases(count,16) = stringData{1,2}(1);
                LoadCases(count,17) = stringData{1,3}(1);
                LoadCases(count,18) = stringData{1,4}(1);
                LoadCases(count,19) = stringData{1,5}(1);
                LoadCases(count,20) = stringData{1,6}(1);
                LoadCases(count,21) = stringData{1,7}(1);
                LoadCases(count,22) = stringData{1,8}(1);
                LoadCases(count,23) = stringData{1,9}(1);
                LoadCases(count,24) = stringData{1,10}(1);
                
                tline  = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f %f %f');
            end
        end
    end % end of while
    fclose(fp);
end

%% FETCH BOX PROPERTIES
pathName = ['export',filesep,'fmwpp01',filesep,'flexible'];
fileName = 'wng_box_properties.fmwpp01';
fid      = 1;
Inputfiles.boxPropertiesFiles{1} = fullfile(folderName,pathName,fileName);

neta = 0;

if ~exist(Inputfiles.boxPropertiesFiles{1},'file')
    error('Unable to find file %s.\n', Inputfiles.boxPropertiesFiles{1});
else
    fp = fopen(Inputfiles.boxPropertiesFiles{1}, 'r');
    skipLine = false;
    
    while ~feof(fp)
        
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
        stringData = textscan(tline,'%f');
        
        if ~isempty(stringData{1,1})
            neta = neta + 1;
            BoxProperties(neta,:) = stringData{1,1}';
        end
        
    end % end of while
    fclose(fp);
    %
    
end

%% FETCH THICKNESS AND CHORD PROPERTIES
fileName = 'wng_addgeodata.fmwpp01';
fid      = 1;
Inputfiles.sectionProfileFiles{1} = fullfile(folderName,pathName,fileName);

neta = 0;

if ~exist(Inputfiles.sectionProfileFiles{1},'file')
    error('Unable to find file %s.\n', Inputfiles.sectionProfileFiles{1});
else
    
    fp = fopen(Inputfiles.sectionProfileFiles{1}, 'r');
    skipLine = false;
    
    while ~feof(fp)
        
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
        stringData = textscan(tline,'%f');
        if ~isempty(stringData{1,1})
            neta = neta + 1;
            AirfoilProperties(neta,:) = stringData{1,1}';
        end
        
    end
    fclose(fp);
end

%% FETCH AERODYNAMICS
list = dir(fullfile(folderName,pathName,'wng_ads_lc*.fmwpp01'));
% Scroll through the different trim cases and store data to array
for i = 1:length(list)
    fileName = list(i).name;
    Inputfiles.aeroCasesFiles{1} = fullfile(folderName,pathName,fileName);
    
    neta = 0;
    
    if ~exist(Inputfiles.aeroCasesFiles{1},'file')
        error('Unable to find file %s.\n', Inputfiles.aeroCasesFiles{1});
    else
        
        fp = fopen(Inputfiles.aeroCasesFiles{1}, 'r');
        skipLine = false;
        
        while ~feof(fp)
            if ~skipLine
                tline = fgetl(fp);
            else
                skipLine = false;
            end
            
            stringData = textscan(tline,'%f');
            
            if ~isempty(stringData{1,1})
                neta = neta + 1;
                Aerodynamic(neta,:) = stringData{1,1}';
            end
            
        end % end of while
        fclose(fp);
        
    end
    Results.Aerodynamic{i} = Aerodynamic;
end

list = dir(fullfile(folderName,pathName,'wng_deform_lc*.fmwpp01'));

for i = 1:length(list)
    %% FETCH DEFORMATIONS
    fileName = list(i).name;
    Inputfiles.wingDeformationsFiles{i} = fullfile(folderName,pathName,fileName);
    
    neta = 0;
    
    if ~exist(Inputfiles.wingDeformationsFiles{i},'file')
        error('Unable to find file %s.\n', Inputfiles.wingDeformationsFiles{i});
    else
        fp = fopen(Inputfiles.wingDeformationsFiles{i}, 'r');
        skipLine = false;
        
        while ~feof(fp)
            if ~skipLine
                tline = fgetl(fp);
            else
                skipLine = false;
            end
            
            stringData = textscan(tline,'%f');
            if ~isempty(stringData{1,1})
                neta = neta + 1;
                Deformation(neta,:) = stringData{1,1}';
            end
            
        end % end of while
        fclose(fp);
        
    end
    Results.Deformation{i} = Deformation;
end

list = dir(fullfile(folderName,pathName,'wng_cld_lc*.fmwpp01'));

for i = 1:length(list)
    %% FETCH INTERNAL LOADS
    fileName = list(i).name;
    Inputfiles.internalLoadsFiles{i} = fullfile(folderName,pathName,fileName);
    
    neta = 0;
    
    if ~exist(Inputfiles.internalLoadsFiles{i},'file')
        error('Unable to find file %s.\n', Inputfiles.internalLoadsFiles{i});
    else
        fp = fopen(Inputfiles.internalLoadsFiles{i}, 'r');
        skipLine = false;
        
        while ~feof(fp)
            if ~skipLine
                tline = fgetl(fp);
            else
                skipLine = false;
            end
            
            stringData = textscan(tline,'%f');
            if ~isempty(stringData{1,1})
                neta = neta + 1;
                InternalLoads(neta,:) = stringData{1,1}';
            end
            
        end
        fclose(fp);
    end
    Results.InternalLoads{i} = InternalLoads;
end

%% FETCH CONTROL SURFACE LOCATION
% IT MAY BE MORE APPOPRIATE TO USE SWITCH CARDS HERE INSTEAD - ORDER MAY
% NOT BE CONSISTENT
dotstr   = findstr('.fame',folderName);
slashstr = findstr(filesep,folderName);

pathName = '__mirror';
fileName = dir(fullfile(folderName,pathName,'*.fm4'));
Aileron  = [];
Flap     = [];
Slat     = [];
Spoiler  = [];
dih      = [];
Inputfiles.controlsLocationsFiles{1} = fullfile(folderName,pathName,fileName(1).name);

count = 0;

if ~exist(Inputfiles.controlsLocationsFiles{1},'file')
    error('Unable to find file %s\n', Inputfiles.controlsLocationsFiles{1});
else
    fp = fopen(Inputfiles.controlsLocationsFiles{1}, 'r');
    skipLine = false;
    
    while ~feof(fp)
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
        tline = tline(tline ~= ' ');
        
        % Search for the wing dihedral
        stringdata=strfind(tline,'<<<WING_>>>');
        
        if ~isempty(stringdata)
            
            while isempty(strfind(tline,'<ROOT_ETA>'))
                
                tline=fgetl(fp); 
                tline=tline(tline~=' ');
            end
            
            stringdata=textscan(tline,'<ROOT_ETA>[%f');
            
            root_wing=stringdata{1,1};
            
            while isempty(strfind(tline,'<DIHEDRAL_DISTRIBUTION_>'))
                
                tline=fgetl(fp);
                tline=tline(tline~=' ');
                
            end
            
            fgetl(fp);tline=fgetl(fp);
            
            count=0;
            
            while isempty(strfind(tline,'<_DIHEDRAL_DISTRIBUTION>'))
                
                count=count+1;
                
                stringdata=textscan(tline,'%f');
                
                if isempty(stringdata{1,1})
                    
                    stringdata=textscan(tline,' &root_wing%f');
                    
                    if ~isempty(stringdata{1,1})
                        
                        dih(count,1)=root_wing;
                        dih(count,2)=stringdata{1,1}';
                        
                    else
                        
                        break
                        
                    end
                else
                    
                    dih(count,:)=stringdata{1,1}';
                end
                tline=fgetl(fp);
            end
        end
        
        stringData = strfind(tline,'<SLAT_GEO_>');
        
        if ~isempty(stringData)
            fgetl(fp);
            tline = fgetl(fp);
            count = 0;
            
            while isempty(strfind(tline,'< _SLAT_GEO'))
                count             = count + 1;
                stringData        = textscan(tline,'%f');
                Slat.geo(count,:) = stringData{1,1}';
                tline             = fgetl(fp);
            end
        end
        
        stringData = strfind(tline,'<FLAP_GEO_>');
        
        if ~isempty(stringData)
            fgetl(fp);
            tline = fgetl(fp);
            count = 0;
            
            while isempty(strfind(tline,'< _FLAP_GEO'))
                count             = count + 1;
                stringData        = textscan(tline,'%f');
                Flap.geo(count,:) = stringData{1,1}';
                tline             = fgetl(fp);
            end
        end
        
        stringData = strfind(tline,'<SPOILER_GEO_>');
        
        if ~isempty(stringData)
            fgetl(fp);
            tline = fgetl(fp);
            count = 0;
            
            while isempty(strfind(tline,'< _SPOILER_GEO'))
                count                = count + 1;
                stringData           = textscan(tline,'%f');
                Spoiler.geo(count,:) = stringData{1,1}';
                tline                = fgetl(fp);
            end
        end
        
        stringData = strfind(tline,'<AILERON_GEO_>');
        
        if ~isempty(stringData)
            fgetl(fp);
            tline = fgetl(fp);
            count = 0;
            
            while isempty(strfind(tline,'< _AILERON_GEO'))
                count                = count + 1;
                stringData           = textscan(tline,'%f');
                Aileron.geo(count,:) = stringData{1,1}';
                tline                = fgetl(fp);
            end
        end
    end
    
    fclose(fp);
end

Aileron.def               = ailDef;
Spoiler.def               = spoilDef * spoilerfactor;
Results.Wing.dih          = dih; % TODO Added this Dario -> Etienne
Results.Control.Aileron   = Aileron;
Results.Control.Flap      = Flap;
Results.Control.Slat      = Slat;
Results.Control.Spoiler   = Spoiler;
Results.LoadCases         = LoadCases;
Results.AirfoilProperties = AirfoilProperties;
Results.BoxProperties     = BoxProperties;

end
