classdef (ConstructOnLoad) FameResult < awi.model.ResultSet
    %
    % Represents content imported from a FAME results file
    
    properties
        
        %Import from where ?
        ImportFolder;
        InternalLoads;
        Aerodynamic;
        Thickness;
        Deformation;
        Stiffness;
        BoxAreas;
        LoadCases
        
    end
    
    properties (Dependent)
        
        Fx
        Fy
        Fz
        Mx
        My
        Mz
        Eta
        
        Legend
        
    end
    
    methods % construction / destruction
        
        function obj = FameResult(varargin)
            %
            % Construct a SizeResult object OBJ, based on the Trim Analysis PAR
                
            %Pass it on
            obj@awi.model.ResultSet(varargin{:});
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, 'Import...', 'import')
            end
            
        end
        
        function import(obj, dn, varargin)
            
            %For FAME4 files, starting point is a directory not a file
            if nargin < 2 || isempty(dn)
                
                %Ask the user
                dn = uigetdir(char(obj.ImportFolder),'Pick a folder ...');
                
                if isempty(dn) || isequal(dn,0)
                    return;
                end
            end
            
            % Read INTERNAL LOADS results
            obj.InternalLoads  = getFameInternalLoads(dn);
            
            % Read AERODYNAMIC results
            obj.Aerodynamic    = getFameAerodynamics(dn);
            
            % Read BOX PROPERTIES results
            obj.Thickness      = getFameThickness(dn);
            
            % Read DEFORMATION results
            obj.Deformation   = getFameDeformations(dn);
            
            % Read BOX STIFFNESS properties
            obj.Stiffness     = getFameStiffness(dn);
            
            % Read BOX AREAS
            obj.BoxAreas      = getFameBoxAreasFcn(dn);
            
            % Read LOADCASE
            pathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible'];
            fileName = 'loadcase_info.txt';
            obj.LoadCases      = getFameLoadCases(fullfile(dn,pathName,fileName));
            
        end
        
    end
    
    methods %get / set
        
        function val = get.Legend(obj)
           
            val = arrayfun(@(x) ['Fame LC ' num2str(x)], [obj.LoadCases.LC],'UniformOutput',false);
            
        end
        
        function val = get.Eta(obj)
            
            val = [obj.InternalLoads.eta];
            
        end
        
        function val = get.Fx(obj)
            
            val = [obj.InternalLoads.Fx];
            
        end
        
        function val = get.Fy(obj)
            
            val = [obj.InternalLoads.Fy];
            
        end
        
        function val = get.Fz(obj)
            
            val = [obj.InternalLoads.Fz];
            
        end
        
        function val = get.Mx(obj)
            
            val = [obj.InternalLoads.Mx];
            
        end
        
        function val = get.My(obj)
            
            val = [obj.InternalLoads.My];
            
        end
        
        function val = get.Mz(obj)
            
            val = [obj.InternalLoads.Mz];
            
        end
        
    end
    
end
