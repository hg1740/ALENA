classdef (ConstructOnLoad) Aircraft < awi.model.Component

    properties (AbortSet, SetObservable)
        
        %Properties related to mass cases
        MTOM = 0;
        MZFM = 0;
        MLM = 0;
        OEM = 0;
        
        %Properties related to Trim of the aircraft
        Trim = '-none-';
        
    end
    
    %Handles to specific components of the 'Aircraft'
    properties 
        %Handle to the 'awi.model.LiftingSurface' objects that are defined
        %as the 'wings'
        Wings
        %Handle to the 'awi.model.Fuselage' objects that are defined as
        %'fuselages'
        Fuselage
    end
    
    %Reference quantities
    properties 
        %Reference wing area (both wings) - Used for MSC.Nastran 
        RefArea
        %Reference chord - Used by MSC.Nastran 
        RefChord
        %Reference wing span (full span) - Used by MSC.Nastran
        RefSpan
    end
    
    properties (Dependent)
        
        ControlSurfaces;
        TrimResults;
        
    end
    
    methods % set / get        
        function set.RefArea(obj, val)  %set.RefArea  
            %set.RefArea Set method for the property 'RefArea'.
            %
            % 'RefArea' must be a positive scalar numeric.            
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'RefArea');
            obj.RefArea = val;
        end        
        function set.RefChord(obj, val) %set.RefChord 
            %set.RefChord Set method for the property 'RefChord'.
            %
            % 'RefChord' must be a positive scalar numeric.            
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'RefChord');
            obj.RefChord = val;
        end
        function set.RefSpan(obj, val)  %set.RefSpan  
            %set.RefSpan Set method for the property 'RefSpan'.
            %
            % 'RefSpan' must be a positive scalar numeric.            
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'RefSpan');
            obj.RefSpan = val;
        end
        
        
        function val = get.ControlSurfaces(obj)
        
            %Simple search
            %   - This won't work because there are subclasses of
            %   'ControlSurface'
            %val = findall(obj, 'Type', 'ControlSurface');
            Collection = flatlist(obj);
            idx = arrayfun(@(o) isa(o, 'awi.model.ControlSurface'), Collection);
            val = Collection(idx);
        end
        
        function val = get.TrimResults(obj)
            %
            %Return a list of applicable trim results by NAME,
            % rather than the actual instance of the trim.
            
            %Start with '-none-' (TODO: what if there was an actual trim result with this name ?)
            val = {'-none-'};
            
            %Need a root object
            if ~isempty(obj.Root)
                
                %Simple search
                trm = findall(obj.Root, 'Type', 'TrimResult');
                
                %Anything ?
                if ~isempty(trm)

                    %Extend list of names
                    val = [val, {trm.Name}];

                end
            
            end
            
        end
        
        function set.Trim(obj, val)
            
            %Make a note - trim set by NAME, rather than an actual TrimResult object
            if isa(val, 'awi.model.TrimResult')
                val = val.Name;
            end
            assert(ischar(val), 'Trim must be specified by name');
            
            %Make a note
            obj.Trim = val;
            
            %Equivalent rotation ?
            rot = [0, 0, 0];
            
            %Need a root object
            if ~isempty(obj.Root)
            
                %Look for corresponding object
                trm = findall(obj.Root, 'Type', 'TrimResult', 'Name', val);
                
                %If found
                if ~isempty(trm)
                    
                    %Set rotation accordingly
                    [r1, r2, r3] = dcm2angle(trm.RotMat);
                    rot = rad2deg([r1, r2, r3]);
                    
                end
                
            end
            
            %Apply
            obj.Rotation = rot;
            
        end
        
    end
    
    methods % construction / destruction
    
        function obj = Aircraft(varargin)
        
            %Pass it on
            obj@awi.model.Component(varargin{:});
            
            %Aircraft does not actually draw itself, just its content
            obj.Visible = false;
            
            %Configure property groups
            obj.addPropertyGroup('Mass Cases', ...
                'MTOM', 'MTOM', 'Maximum take-off mass', [], ...
                'MZFM', 'MZFM', 'Maximum zero-fueled mass', [], ...
                'MLM' , 'MLM',  'Maximum landing mass', [], ...
                'OEM' , 'OEM',  'Operating empty mass', []);
            
            %Configure property groups
            obj.addPropertyGroup('Trim', ...
                'Trim', 'Trim', 'Trim result applied to Aircraft', @obj.TrimResults);
            
            %An aircraft describes a collection of finite element
            %components
            obj.IsFEComponent = true;
            
        end
                
    end
    
    methods % Analyses
    
        function M = mass(obj, varargin)
            
            %Two alternative approaches
            if true
                
                %COULD use the TotalMass property (which traverses the model)
                M = obj.TotalMass;
                
            else
                
                %OR do it the hard way
                %For mass, easier to work with a flat list
                L = flatlist(obj);
                
                %Total mass of all components in model
                M = sum([0, L.Mass]);
                
            end
            
            %If caller did not want it back
            if nargout == 0
                
                %Tell us about it
                msgbox(obj, ['Total mass = ', num2str(M)]);
                
            end
            
        end
        
    end
    
    methods % converting to FE model
        function FEModel = convertThisToFE(obj, FEModel, varargin)
           %convertThisToFE Converts the object to a collection of 
           %Finite Element (FE) data objects. 
           %
           % The 'awi.model.Aircraft' just parents a set of
           % 'awi.mixin.FEable' objects therefore there is no need to
           % generate any FE data for this object.
           
           %Pass to superclass
           FEModel = convertThisToFE@awi.mixin.FEable(obj, FEModel, varargin{:});
           
        end 
    end
    
    methods (Static)
        function FEModel = blankFEModel()
           FEModel = awi.fe.FEModelAircraft; 
        end
        function Aircraft = defaultAircraft(aircraftType)
            %defaultAircraftModel Returns an instance of the
            %'awi.model.Aircraft' class which describes a default version
            %of an aircraft. If the aircraft type is not specified the
            %default is an A320.
            %
            % Available aircraft types:
            %   - A320 -> 'A320 (default)'
            
            validAircraftType = {'A320 (default)'};
            
            if nargin < 1
               aircraftType = 'A320 (default)';                
            end
            
            Aircraft = awi.model.Aircraft('Name', aircraftType);
            
            %Assign the relevant data
            switch aircraftType
                case 'A320 (default)'
                    Aircraft = awi.model.Aircraft.makeA320(Aircraft);
                otherwise
                    validatestring(aircraftType, validAircraftType, ...
                        'awi.model.Aircraft.defaultAircraft', 'the aircraft type ("aircraftType")');
            end                          
            
            %Build the aircraft
            build(Aircraft);
            
        end
        function Aircraft = makeA320(Aircraft)
            %makeA320 Generates a version of the A320 aircraft.
            %
            %The dimensions for the A320 are taken from the Wikipeida page
            % (https://en.wikipedia.org/wiki/Airbus_A320_family#Overview)
            % or have been estimated from the GA.
            
            %% Make a default fuselage
            F = awi.model.BluffBody('Name', 'Fuselage');
            
            %Use the radius and eta values from Dario Calderon's default
            %aircraft fuselage
            F.Eta = [0;0.005268;0.010536;0.015805;0.021073;...
                0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
                0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
                0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
                0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';
            
            F.Radius = [0;0.3844030;0.565081;0.707928;0.830682;0.940375;...
                1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
                1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
                2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
                1.75;1.6;1.4;1.2;1.0;0]';
            
            %Place it on the origin
            F.Origin = [0, 0, 0];
            
            %Length?
            F.Length = 37.57;
            
            %% Wings
            
            StbdWing = awi.model.LiftingSurface('Name', 'Starboard Wing');            
            StbdWing.ActiveSet   = 'sSet';
            
            %Position along the fuselage
            StbdWing.SOffset = 0.45;
            StbdWing.ZOffset = -1.0;
            
            %Wing dimensions
            StbdWing.SpanVector  = 'Y';
            StbdWing.Span        = 35.8/2;
            StbdWing.LESweep     = [25, 25];
            StbdWing.LESweep_eta = [0, 1];
            StbdWing.TESweep     = [0, 10, 10];
            StbdWing.TESweep_eta = [0, 0.3, 1];
            StbdWing.RootChord   = 5;
            
            PortWing = awi.model.LiftingSurface('Name', 'Port Wing', ...
                '-Mirror', StbdWing);
            
            %% Add the parts to the model
            Aircraft.add(F);
            F.add(StbdWing);
        end
    end

end
