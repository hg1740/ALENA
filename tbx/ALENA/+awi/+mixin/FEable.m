classdef (ConstructOnLoad) FEable < matlab.mixin.SetGet
    %FEable Denotes a class that can be converted to a Finite Element (FE)
    %model.
    %
    % Provides a method 'convertToFE' which invokes the method
    % 'convertThisToFE' on the object before recursing through all the 
    % children.
            
    properties
        %Controls what type of elements will be used to generate the model.
        %   - '1D' : Generates a beam-stick version of the aircraft
        %   - '2D' : Generates a shell-element model of the aircraft
        %            including ribs and spars, etc.
        ModelType = '1D';
    end
    
    %Structural mesh properties
    properties %(AbortSet, SetObservable)
        %Length of a beam element
        BeamLength = 0.8;
        %Number of beam elements along the beam
        NumBeamElem
        %Width of a shell element
        ShellWidth
        %Aspect-Ratio (AR) of a shell element
        ShellAR = 1;
    end
    
    %Aero mesh properties
    properties
        %Length of an aerodynamic panel in the global X-axis
        AeroPanelLength = 0.3;
        %Number of aerodynamic panels along the chord
        NumAeroPanel        
        %Aspect-Ratio (AR) of an aerodynamic panel
        AeroPanelAR = 1;
        %Flag to model control surfaces
        ModelControlSurf = false;
        %Flag to model control surface structural connections
        ModelControlSurfStructure = false;
    end
    
    %Controlling the FE conversion process
    properties
        %Logical flag indicating whether the beam properties (I11, I22, GJ,
        %etc.) or beam cross-section should be used to define the sectional
        %properties.
        UseBeamCrossSection = false;
        %Logical flag for modelling the aerodynamic panels 
        ModelAero = true;
        %Logical flag for modelling aerodynamic control surfaces
        ModelAeroControls = true;
        %Logical flag for generating the lifting surface planform nodes
        GenerateLiftSurfPlanform = true;
        %Logical flag for generating the cross-section nodes using rigid 
        %links
        GenerateCrossSectionNodes = false;
        %Specifies whether RBE2 (small-angle) or RBE3 elements should be
        %used for splining the cross-section nodes to the beam axis.
        CrossSectionRigidBarType = 'RBE2';
    end
    
    properties (Dependent)
        %Logical flag indicating whether the object can be converted to a
        %Finite Element representation
        CanConvertToFE
    end
    
    properties (SetAccess = protected)
        %Logical flag which denotes whether the object is a primary
        %component of the Finite Element model
        IsFEComponent = false;
    end
    
    methods % set / get 
        function set.BeamLength(obj, val)               %set.BeamLength      
            %set.BeamLength Set method for the property 'BeamLength'.
            
%             validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
%                 'nonnan', 'finite', 'real', 'nonempty'}, class(obj), ...
%                 'BeamLength')
            obj.BeamLength = val;
        end
        function set.NumBeamElem(obj, val)              %set.NumBeamElem     
            %set.NumBeamElem Set method for the property 'NumBeamElem'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'integer', 'nonnan', 'finite', 'real'}, ...
                class(obj), 'BeamLength')
            obj.NumBeamElem = val;
        end
        function set.AeroPanelLength(obj, val)          %set.AeroPanelLength 
            %set.BeamLength Set method for the property 'AeroPanelLength'.
            
%             validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
%                 'nonnan', 'finite', 'real'}, class(obj), ...
%                 'AeroPanelLength')
            obj.AeroPanelLength = val;
        end
        function set.NumAeroPanel(obj, val)             %set.NumAeroPanel
            %set.NumAeroPanel Set method for the property 'NumAeroPanel'.

            %             validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
            %                 'integer', 'nonnan', 'finite', 'real'}, ...
            %                 class(obj), 'NumAeroPanel')
            validateattributes(val, {'numeric'}, {'positive','integer'}, ...
                class(obj), 'NumAeroPanel')
            obj.NumAeroPanel = val;
        end
        function set.AeroPanelAR(obj, val)              %set.AeroPanelAR     
            %set.AeroPanelAR Set method for the property 'AeroPanelAR'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'nonnan', 'finite', 'real', 'nonempty'}, class(obj), ...
                'AeroPanelAR')
            obj.AeroPanelAR = val;
        end
        function set.ModelAero(obj, val)                %set.ModelAero       
            %set.ModelAero Set method for the property 'ModelAero'.
            %
            % 'ModelAero' must be a scalar logical
            
            validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, ...
                class(obj), 'ModelAero');
            obj.ModelAero = val;
        end
        function set.ModelAeroControls(obj, val)        %set.ModelAeroControls       
            %set.ModelAeroControls Set method for the property
            %'ModelAeroControls'.
            %
            % 'ModelAeroControls' must be a scalar logical
            
            validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, ...
                class(obj), 'ModelAeroControls');
            obj.ModelAeroControls = val;
        end
        function set.GenerateLiftSurfPlanform(obj, val) %set.GenerateLiftSurfPlanform
            validateattributes(val, {'logical'}, {'scalar', 'nonempty'}, ...
                class(obj), 'GenerateLiftSurfPlanform');
            obj.GenerateLiftSurfPlanform = val;
        end
        function set.CrossSectionRigidBarType(obj, val) %set.CrossSectionRigidBarType
            val = validatestring(val, {'RBE2', 'RBE3'}, class(obj), 'CrossSectionRigidBarType');
            obj.CrossSectionRigidBarType = val;
        end
        function val = get.CanConvertToFE(obj)          %get.CanConvertToFE  
            %get.CanConvertToFE Get method for the dependent property
            %'CanConvertToFE'.
            
            %Pass it on
            val = canConvertToFE(obj);
        end
    end
    
    methods % converting to FE model       
        function FEModel = convertToFE(obj, FEModel, varargin)
            %convertToFE Call the class specific method 'convertThisToFE'
            %for this object before passing the method onto the objects
            %children.
            
            if nargin < 2  || isempty(FEModel) %Create a default model
                %Returns a blank/empty instance of 'awi.model.FEModel'
                %which will be screened in the parent conversion process.
                FEModel = obj.blankFEModel;
                FEModel.IsComponent = obj.IsFEComponent;
            end
            if obj.IsFEComponent %Extend the collection for a new component
                FEModel(end + 1)  = obj.blankFEModel;
                FEModel(end).Name = obj.Name;
                FEModel(end).IsComponent = obj.IsFEComponent;
                %Store a reference to the geometry object
                FEModel(end).GeometryObject = obj;
            end
            
            %Parse the extra arguments
            [bRecurse, varargin] = parseTokens(obj, varargin);
            function [bRecurse, argIn] = parseTokens(obj, argIn)
                
                %Find special tokens
                idxChar = cellfun(@ischar, argIn);
                indChar = find(idxChar);
                tok     = argIn(idxChar);
                index   = indChar(cellfun(@(x) isequal(x(1), '-'), tok));

                %Recursion?
                bRecurse = ismember('-UnderRecursion', tok);
                
                %Remove tokens from the list
                argIn(index) = [];
                
                %Log function?
                idxLog = ismember('LogFcn', tok);
%                 if ~any(idxLog)
% %                     if isa(obj, 'mvc.mixin.UiTools')
% %                         logFcn = progressdlg(obj, 'FE conversion process');
% %                     else
% %                         logFcn = @(str) fprintf('%s\n', str);
% %                     end
%                     logFcn = @(str) strcat(str, ' '); 
%                     argIn = [argIn, {'LogFcn', logFcn}];
%                 end
                
            end            
            
            %Has the parent FEModel been supplied 
            ind = find(ismember(varargin(1 : 2 : end), 'Parent'), 1);
            if ~isempty(ind)
                par = varargin{ind + 1};
                %Only add the parent if it has been populated with data
                if par.HasFEData
                    FEModel(end).Parent = par;
                end
            end
            
            %Convert this object...
            FEModel(end) = convertThisToFE(obj, FEModel(end), varargin{:});            
            
            %...and any children...
            idx = arrayfun(@(x) isa(x, 'awi.mixin.FEable'), obj.Children);
            ch = arrayfun(@(c) convertToFE(c, [], '-UnderRecursion', ...
                'Parent', FEModel(end), varargin{:}), obj.Children(idx), 'Unif', false);
                       
            %...make sure children are connected to their parents... 
            %   - But remove empty children before we do!
            ch = [ch{:}];
            if ~isempty(ch)
                idx = [ch.HasFEData];
                FEModel(end).Children = ch(idx);
            end

            %...remove any empty FEModels
            %   - If we are at the top of the tree and there is no data in
            %     the FE model then step down a level in the hierachy and 
            %     use that FE model... TODO : Add this in            
            hasData = [FEModel.HasFEData];
            FEModel = FEModel(hasData);

            %Anything?
            if isempty(FEModel)
                %Force empty matrix
                FEModel = [];
            else           
                %...if we are at the top level then we need to establish 
                %the connectivity between components in the model and 
                %update the ID numbers for all objects in the 
                %component(s)/model
                if ~bRecurse                
                    %Ensure the parent and child are connected somehow...
                    connectParentAndChild(FEModel);
                    %Resolve 'Connector' objects
                    updateConnectors(FEModel);
                    %Assign unique ID numbers to every FE object
                    assignIDnumbers(flatlist(FEModel));
                    %Close the logging function
                end
            end
           
        end        
        function FEModel = convertThisToFE(obj, FEModel, varargin)
           %convertThisToFE Converts the object to a collection of 
           %Finite Element (FE) data objects. 
           %
           % This method is a placeholder method. The actual conversion to
           %a FE model will be done at the subclass level. 
           
           assert(obj.CanConvertToFE, sprintf(['Unable to convert the ', ...
               'object ''%s'' to a Finite Element representation.']    , ...
               obj.Name));
           
           %Custom log function?
           p = inputParser;
           p.KeepUnmatched = true;
           addParameter(p, 'LogFcn', @(s) fprintf('%s\n', s), @(x)isa(x, 'function_handle'));
           parse(p, varargin{:});
           
           if nargin < 2 %Create a default model
               FEModel = awi.fe.FEModel;
               FEModel.IsComponent = obj.IsFEComponent;
           end  
           
           FEModel.LogFcn = p.Results.LogFcn;
           
        end
        function tf = canConvertToFE(~)
            %canConvertToFE Returns a logical (true/false) stating whether
            %the conversion process to an FE representation can take place.
            %
            % More comprehensive checks can be initiated at the subclass
            % level. For now, just allow code to continue.
            
            tf = true;
        end
    end
    
    methods(Static)
        function FEModel = blankFEModel()
           FEModel = awi.fe.FEModelAircraft; 
        end
    end
        
end

