%% INPUTFILES Input files used from FAME
%
% This object contains all the files that were used from FAME. Only one of
% the properties can be set by the user - |wbdCorrFiles|. This property
% contains the name of the Excel correction file that is used to correct
% certain mass values and positions. 

%   Copyright 2016 University of Bristol
%   Private function.
classdef Inputfiles
    
    properties
        loadCaseFiles
        boxPropertiesFiles
        boxStiffnessFiles
        boxAreaFiles
        sectionProfileFiles
        aeroCasesFiles
        wingDeformationsFiles
        internalLoadsFiles
        controlsLocationsFiles
        nastranFiles
        adamsFiles
        wbdFiles
        wbdCorrFiles
        fuelMassFiles
        aeroPlanformFiles
        fameEngineFiles
        thicknessFiles
        fm4File
        TUXFile
        airfoilFiles
        NeofuelMassFiles
    end
    
    methods
        function obj = Inputfiles(~)
            
        end
        
        function obj = setFields(obj,INPUTFILES)
            f = fieldnames(INPUTFILES);
            
            for i = 1:numel(f)
                obj.(f{i}) = INPUTFILES.(f{i});
            end
        end
        
        function display(obj)
            f = fieldnames(obj);
            l = max(cellfun(@numel,f));
            
            if ~isempty(obj.wbdCorrFiles)
                fprintf('\n')
                fprintf(' USER DEFINED INPUT FILES :\n')
                fprintf('\n')
                fprintf('%s %s\n',padString('wbdCorrFiles',l),obj.wbdCorrFiles);
            end
            
                fprintf('\n\n')            
            
                fprintf(' FAME GENERATED INPUT FILES :\n')
                fprintf('\n')
            
            for i = 1:numel(f)
                
                if strcmp(f{i},'wbdCorrFiles')
                    continue
                end
                
                if iscell(obj.(f{i}))
                    nCells = numel(obj.(f{i}));
                    
                    for j = 1:nCells
                        if j == 1
                            fprintf('%s %s\n',padString(f{i},l),obj.(f{i}){1});
                        else
                            fprintf('%s %s\n',padString('',l),obj.(f{i}){j});
                        end
                    end
                else
                    fprintf('%s %s\n',padString(f{i},l),obj.(f{i}));
                end
            end
        end
    end
    
end

function str = padString(str,l)

nBlanks = l - numel(str) + 7;
str = sprintf('%s%s :',blanks(nBlanks),str);

end
