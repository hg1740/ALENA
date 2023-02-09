%GETFAMELOADCASES: A function that extracts the loadcase information
%
% Date of Creation: 12/2015
%
% Authors:  Dario Calderon - University of Bristol
%           Etienne Coetzee - Airbus (FPO)
%
%
function [results,fName] = getFameLoadCases(fName,OEW)

% The pathname is essentially hardcoded
% pathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible'];
% 
% %% FETCH MORE LOAD CASES
% fileName = 'loadcase_info.txt';
% 
% Inputfiles.loadCaseFiles{1} = fullfile(folderName,pathName,fileName);

if nargin == 1
    OEW = 0;
end
    
spoilcount      = 0;
ailcount        = 0;
flapcount       = 0;
fuelcount       = 0;

spoilerReductionfactor = 1;

if ~exist(fName,'file')
    error('Unable to find file %s.\n', fName);
else
    fp = fopen(fName, 'r');
    
    skipLine = false;
    
    while ~feof(fp)
        
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        
%         % Fetch the LoadCase label:
%         
%         stringData = strfind(tline,'Loadcase: %f Loadcasetype: %s');
%         if ~isempty(stringData)
%             stringData      = textscan(tline,'Loadcase: %i Loadcasetype: %s');
%             
%         end
%         
        % Fetch the fuel mass:
        stringData = strfind(tline,' Betankungszustand:      !!! Both Sides !!!');
        if ~isempty(stringData)
            tline = fgetl(fp);
            tline = fgetl(fp);
            tline = fgetl(fp);
            tline = fgetl(fp);
            tline = fgetl(fp);
            tline = fgetl(fp);
            tline = fgetl(fp);
            stringData      = textscan(tline,' Tankmenge(SOLL), gesamt [kg] :  %f');
            fuelcount        = fuelcount + 1;
            FuelDef(fuelcount,:) = stringData{1,1};
            
        end
        
        % Fetch control surface deflections for each load case
        stringData = strfind(tline,'Aileron deflection [degree]:');
        
        if ~isempty(stringData)
            %stringData      = textscan(tline,'Aileron deflection [degree]: %f');
            stringData      = regexp(tline,'-?\d*[.]\d*','Match');
            ailcount        = ailcount + 1;
            %ailDef(ailcount,:) = stringData{1,1};
            ailDef(ailcount,:) = str2double(stringData{1,1});
            
        end
        
        stringData = strfind(tline,'Flap deflection    [degree]:');
        
        if ~isempty(stringData)
            stringData     = regexp(tline,'-?\d*[.]\d*','Match');
            %stringData     = textscan(tline,'Flap deflection    [degree]: %f %f');
            flapcount      = flapcount + 1;
            for i = 1:length(stringData)
                %FlapDef(flapcount,i) = stringData{1,i};
                FlapDef(flapcount,i) = str2double(stringData{1,i});
            end
            
            
        end
        
        stringData = strfind(tline,'Spoiler deflection [degree]:');
        
        if ~isempty(stringData)
            stringData     = regexp(tline,'-?\d*[.]\d*','Match');
            %stringData     = textscan(tline,'Spoiler deflection [degree]: %f %f %f %f %f %f');
            spoilcount     = spoilcount + 1;
            for i = 1:length(stringData)
                %spoilDef(spoilcount,i) = stringData{1,i};
                spoilDef(spoilcount,i) = str2double(stringData{1,i});
            end
            
            
        end
        
        % Fetch flight conditions
        stringData = strfind(tline,'LC   AC_mass');
        
        if ~isempty(stringData)
            nTrim = 0;
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                nTrim = nTrim + 1;
                
                LoadCases(nTrim).LC             = stringData{1,1}(1);
                LoadCases(nTrim).AC_mass        = stringData{1,2}(1);
                LoadCases(nTrim).Loadfactor 	= stringData{1,3}(1);
                LoadCases(nTrim).Safetyfactor   = stringData{1,4}(1);
                LoadCases(nTrim).AC_cg          = stringData{1,5}(1);
                LoadCases(nTrim).PitchAngle     = stringData{1,6}(1);
                LoadCases(nTrim).Rollangle      = stringData{1,7}(1);
                if ~isempty(stringData{1,8})
                    LoadCases(nTrim).EngineThrust   = stringData{1,8}(1);
                else
                    LoadCases(nTrim).EngineThrust   = [];
                end
                
                tline = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            end
        end
        
        stringData = strfind(tline,'LC   Machzahl');
        
        if ~isempty(stringData)
            nTrim = 0;
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                nTrim = nTrim + 1;
                
                LoadCases(nTrim).Mach           = stringData{1,2}(1);
                LoadCases(nTrim).Velocity       = stringData{1,3}(1)*1000/(60*60);
                LoadCases(nTrim).Altitude       = stringData{1,4}(1);
                LoadCases(nTrim).Temperature    = stringData{1,5}(1);
                LoadCases(nTrim).AirDensity     = stringData{1,6}(1);
                LoadCases(nTrim).Staticpressure = stringData{1,7}(1);
                LoadCases(nTrim).Dynamicpressure= stringData{1,8}(1);
                
                tline = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f');
            end
        end
        
        stringData = strfind(tline,'LC   AC-Lift');
        
        if ~isempty(stringData)
            nTrim = 0;
            fgetl(fp);
            fgetl(fp);
            fgetl(fp);
            tline = fgetl(fp);
            stringData = textscan(tline,'%f %f %f %f %f %f %f %f %f %f');
            
            while ~isempty(stringData{1,1})
                nTrim = nTrim + 1;
                
                LoadCases(nTrim).AC_lift    = stringData{1,2}(1);
                LoadCases(nTrim).Wing_lift  = stringData{1,3}(1);
                LoadCases(nTrim).Tail_lift  = stringData{1,4}(1);
                
                tline  = fgetl(fp);
                stringData = textscan(tline,'%f %f %f %f %f %f %f %f %f %f');
            end
        end
    end % end of while
    fclose(fp);
end

for i = 1:size(ailDef,1)
    
    LoadCases(i).aileron = ailDef(i,:);
    LoadCases(i).spoiler = spoilDef(i,:)*spoilerReductionfactor; % FIX THIS Spoler reduction factor from fame-w.fmc
    LoadCases(i).flap    = FlapDef(i,:);
    
    % Add the rudder loadcase for future reference
    LoadCases(i).rudder  = 0;
    
    % Add the fuel mass for each loadcase
    LoadCases(i).fuel    = FuelDef(i);
    
    % Payload Mass??
    if OEW == 0
        LoadCases(i).PayloadMass = [];
    else
        LoadCases(i).PayloadMass = LoadCases(i).AC_mass - LoadCases(i).fuel - OEW;
    end
    
end


results = LoadCases;

end