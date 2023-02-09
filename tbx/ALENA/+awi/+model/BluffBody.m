classdef (ConstructOnLoad) BluffBody < awi.model.Beam
    %BluffBody Defines a generic bluff-body.
    %
    %   - A bluff-body is assumed to be aligned with the global x-axis
    %   - The beam line of a bluff-body can vary in the global-X and
    %     global-Z directions but must be aligined with the X-Z plane.
    %
    %
    % TODO - Add engine class (subclass of bluff body) and add Cant Angle
    % and Pitch Angle to the property list (same as FAME I think).
    
    %Appearance
    properties (AbortSet, SetObservable)
        %Appearance described in terms of face colour and transparency
        FaceColor    = [0, 1, 0];
        FaceAlpha    = 0.25;
        FaceLighting = 'gouraud'
    end
    
    %Length Set
    properties (AbortSet, SetObservable)
        %Total length of the bluff-body
        Length
        %Normalised position along the bluff-body axis
        Eta
        %Radius of the bluff-body at 'Eta'        
        Radius
        %Offset in the positive global y-axis
        DeltaZ = [0, 0];
        %Normalised position along the bluff-body axis of 'DeltaZ'
        DeltaZ_eta = [0, 1];
    end
    
    %Coordinate Set
    properties (AbortSet, SetObservable)
        %X-Coordinates of the beam in the global reference frame
        AbsX 
        %Y-Coordinates of the beam in the global reference frame
        AbsY
        %Z-Coordinates of the beam in the global reference frame
        AbsZ
    end
    
    %Analytic Properties
    properties (AbortSet, SetObservable)
        %Number of points defining the circumference of the bluff-body
        NPoints = 20;
    end
    
    %Shadow properties
    properties (Dependent)
       Eta_     %Normalised position along the bluff-body axis accounting for Radius and z-offset eta values
       Radius_  %Bluff-body radius at 'Eta_'
       DeltaZ_  %Offset in the global z-axis at 'Eta_'
    end
    
    methods % set / get 
        function val = get.Eta_(obj)
            %get.Eta_ Get method for the property 'get.Eta_'.
            %
            % 'obj.Eta_' is the finest distribution required to define the
            % position of the stick based on the length and the z-offset.
            
            eta = get(obj, {'Eta', 'DeltaZ_eta'});
            val = unique([cat(2, eta{:}), 0, 1]);
        end
        function val = get.Radius_(obj)
            %get.Radius_ Get method for the property 'get.Eta_'.
            %
            % 'obj.Radius_' is the radius of the bluff body defined at
            % 'obj.Eta_'.
            
            if isempty(obj.Eta) || isempty(obj.Radius)
                val = [0, 0];
            else
                val = interp1(obj.Eta, obj.Radius, obj.Eta_);
            end
        end
        function val = get.DeltaZ_(obj)
            %get.DeltaZ_ Get method for the property 'DeltaZ_'.
            %
            % 'obj.DeltaZ_' is the offset in the global-Z direction at the
            % normalised position 'Eta_'.
            
            val = interp1(obj.DeltaZ_eta, obj.DeltaZ, obj.Eta_);
        end
    end
    
    methods % construction / destruction
    
        function obj = BluffBody(varargin)
        
            %Pass it on
            obj@awi.model.Beam(varargin{:});  
            
            %'Radius'     , 'Radius of the bluff body at ''Eta''', ...
            %'Radius', 'Radius of the bluff body at ''XData'', ''YData'', ''ZData'''
            %Add parameter set
            obj.addParameterSet('lSet', ...
                'DisplayName' , 'Length Set', ...
                'Precedence'  , 1, ...
                'Description' , ['Defines the bluff body using the overall ', ...
                'length and the radius at normalised positions along the length.'], ...
                'Parent'     , '', ...
                'Length'     , 'Total length of the bluff-body', ...
                'Eta'        , 'Normalised position along the bluff-body axis', ...                
                'DeltaZ'    , 'Offset in the positive global Z-axis', ...
                'DeltaZ_eta', 'Normalised position along the bluffbody of ''DeltaZ''');            
            obj.addParameterSet('cSet', ...
                'DisplayName', 'Coordinate Set', ...
                'Precedence' , 2, ...
                'Description', 'Defines the coordinates of the bluff body stick', ...
                'AbsX'  , 'X-coordinate of the bluff-body stick in the global frame', ...
                'AbsY'  , 'Y-coordinate of the bluff-body stick in the global frame', ...
                'AbsZ'  , 'Z-coordinate of the bluff-body stick in the global frame');
        
            %Extend property groups
            obj.addPropertyGroup('Appearance', ...
                'FaceColor', 'Colour applied to planform and elevation', ...
                'FaceAlpha', 'Transparency of planform and elevation', ...
                'FaceLighting', '');
            
            %Bluff-body is a primary component in a FE model
            obj.IsFEComponent = true;
            
        end
        
    end
        
    methods % class building
        
        function build_lSet(obj)
            %build_cSet Build method for the parameter set 'lSet'.
            
            %Set 'XData', 'YData' & 'ZData'
            
            obj.XData = obj.Eta_ .* obj.Length;
            obj.YData = zeros(size(obj.XData));
            obj.ZData = obj.DeltaZ_;
     
        end
        
        function build_cSet(obj)
            %build_cSet Build method for the parameter set 'cSet'.
                        
            %Update the 'Origin' property
            if isempty(obj.Origin)
                obj.Origin   = [obj.AbsX(1), obj.AbsY(1), obj.AbsZ(1)];
            end

            %Invoke get method once!
            pos = obj.Position;
            
            %Set the 'XData', 'YData' &  'ZData' relative to 'obj.Position'
            %   - i.e. in the local stick frame
            obj.XData = obj.AbsX - pos(1);
            obj.YData = obj.AbsY - pos(2);
            obj.ZData = obj.AbsZ - pos(3);
        end
        
        function update_lSet(obj)
            %update_lSet Update method for the parameter set 'lSet'.
            
            obj.Length  = obj.XData(end) - obj.XData(1);
            obj.Eta     = [0, obj.XData(2:end) ./ obj.Length];
            obj.DeltaZ = [0, obj.ZData(2:end) - obj.ZData(1)];
            obj.DeltaZ_eta = obj.Eta;
        end
        
        function update_cSet(obj)
            %update_cSet Updates method for the parameter set 'cSet'.
            
            %Invoke get method once!
            pos = obj.AbsPosition;
            
            %Set the 'XData', 'YData' &  'ZData' to the global coordinate 
            %system
            obj.AbsX = obj.XData + pos(1);
            obj.AbsY = obj.YData + pos(2);
            obj.AbsZ = obj.ZData + pos(3);
            
        end
        
    end
    
    methods % converting to FE model
        function FEModel = convertThisToFE(obj, FEModel, varargin)
            %convertThisToFE Class specfic method for converting an
            %instance of 'awi.model.BluffBody' into a collection of FE data
            %objects.
            %
            % Steps:
            %   1. Determine which points along the beam/Stick should be
            %   passed to the FE-converter for the 'awi.model.Beam' class.
            %   2. Pass these points to the 'awi.model.Beam' FE-converter.
            %   3. Define the aerodynamic panels...
            %
            % TODO - Update this to handle a generic cross-section. Should
            % probably hand it down to the Beam class.
            
            %Start with base-class
            FEModel = convertThisToFE@awi.mixin.FEable(obj, FEModel);            
            
            %Points from geometry that MUST be included in the beam mesh
            xG   = obj.XData;
            yG   = obj.YData;
            zG   = obj.ZData;
            rG   = obj.RData;
            etaG = rG ./ rG(end);
                        
            %Pass to 'awi.model.Beam' method...
            FEModel = convertThisToFE@awi.model.Beam(obj, FEModel, ...
                'etaR', etaG, 'xR', xG, 'yR', yG, 'zR', zG);
            
            %Add FE objects to describe the representative aerodynamic
            %surface of the bluff-body
            if obj.ModelAero && ~any(isnan(obj.Radius))
                
                %Make the aerodynamic panels
                AeroPanel = i_makeAeroPanels(obj);
                
                %Make the nodes, rigid bars, strucutal sets, etc.
                
                
                %% Making Rigid elements (charles)
                FuselageBeamNodes = FEModel.BeamNodes;
                FuselageBeamNodes = [FuselageBeamNodes(1, :), FuselageBeamNodes(2, end)];
                FuselagebeamX     = [FuselageBeamNodes.X];
                nGrid     = size( FuselagebeamX, 2);
                
                %Preallocate
                Nodesupper = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
                Nodeslower = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
                Nodesleft = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
                Nodesright = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
                RB_vertical      = arrayfun(@(~) awi.fe.RigidBar, 1 : nGrid   , 'Unif', false);
                RB_horizontal      = arrayfun(@(~) awi.fe.RigidBar, 1 : nGrid   , 'Unif', false);
       
                %Set     = awi.fe.StructuralSet;
                Nodesupper = horzcat(Nodesupper{:});
                Nodeslower = horzcat(Nodeslower{:});
                Nodesleft = horzcat(Nodesleft{:});
                Nodesright = horzcat(Nodesright{:});
                RB_vertical = horzcat(RB_vertical{:});
                RB_horizontal = horzcat(RB_horizontal{:});
                
                Coordupper = zeros(3, nGrid);
                Coordlower = zeros(3, nGrid);
                Coordleft = zeros(3, nGrid);
                Coordright = zeros(3, nGrid);
                
                X_ = obj.Eta*obj.Length+obj.Origin(1);
                R_=obj.Radius;
                
                upper_ = [X_ ; zeros(1,length(X_)) ; obj.Origin(3) + R_];
                lower_ = [X_ ; zeros(1,length(X_)) ; obj.Origin(3) - R_];
                left_ = [X_ ; obj.Origin(2) - R_ ; zeros(1,length(X_))];
                right_ = [X_ ; obj.Origin(2) + R_ ; zeros(1,length(X_))];
                
                              
                Coordupper(3, :)      = interp1(upper_(1, :), upper_(3, :), FuselagebeamX(1, :));
                Coordupper([1, 2], :) = [FuselagebeamX(1, :) ; FuselagebeamX(2, :)];
                
                Coordlower(3, :)      = 2*obj.Origin(3) - Coordupper(3, :);
                Coordlower([1, 2], :) = [FuselagebeamX(1, :) ; FuselagebeamX(2, :)];
                
                Coordleft(2, :)      = interp1(left_(1, :), left_(2, :), FuselagebeamX(1, :));
                Coordleft([1, 3], :) = [FuselagebeamX(1, :) ; FuselagebeamX(3, :)];
                
                Coordright(2, :)      = 2*obj.Origin(2) - Coordleft(2, :);
                Coordright([1, 3], :) = [FuselagebeamX(1, :) ; FuselagebeamX(3, :)];
                               

                set(Nodesupper, {'X'}, num2cell(Coordupper , 1)');
                set(Nodeslower, {'X'}, num2cell(Coordlower, 1)');
                set(Nodesleft, {'X'}, num2cell(Coordleft, 1)');
                set(Nodesright, {'X'}, num2cell(Coordright, 1)');
                
                %Define the rigid bar independant/dependant nodes
                set(RB_vertical, {'NodesI'}, num2cell(FuselageBeamNodes')); %Conversion to double from cell is not possible.
                set(RB_vertical, {'NodesD'}, num2cell([Nodesupper ; Nodeslower ; Nodesleft ; Nodesright], 1)');
                
%                 set(RB_horizontal, {'NodesI'}, num2cell(FuselageBeamNodes'));
%                 set(RB_horizontal, {'NodesD'}, num2cell([Nodesleft ; Nodesright], 1)');
                
%                 set(RB_horizontal, 'CN', 123456);
                set(RB_vertical, 'CN', 123456);
                
                %%-----------------------------------------------------
                
                %Add objects to the FE model
%                 addFEData(FEModel, AeroPanel,Nodesupper,Nodeslower, Nodesleft, Nodesright, RB_vertical);
                
                %% temperary remove aeropanels on body 
                
                addFEData(FEModel,Nodesupper,Nodeslower, Nodesleft, Nodesright, RB_vertical);
                %% Define aeropanels(charles)
                
                  %No need to proceed if we aren't modelling the aero
                  if ~obj.ModelAero
                      return
                  end
                  
                  
                  %Create a single 'awi.fe.AeroProp' object
                  AeroProp = awi.fe.AeroProp;
                  
                  AeroPanels = [AeroPanel];
                  nAeroSeg = numel(AeroPanels);
                  
                  %Assign a reference to the aerodynamic properties
                  set(AeroPanels, 'AeroProp', AeroProp);
                  
                  %Preallocate
                  StrucSet = arrayfun(@(~) awi.fe.StructuralSet    , 1 : nAeroSeg, 'Unif', false);
                  AeroSet  = arrayfun(@(~) awi.fe.AeroPanelSet     , 1 : nAeroSeg, 'Unif', false);
                  Spline   = arrayfun(@(~) awi.fe.AeroelasticSpline, 1 : nAeroSeg, 'Unif', false);
                  
                  StrucSet = horzcat(StrucSet{:});         
                  AeroSet  = horzcat(AeroSet{:});
                  Spline   = horzcat(Spline{:});
                  
                  %Set up coordinate matrices - Quicker outside the loop
                  %   - Need absolute coordinate of beam coordinate and panel LE
                  %     coordinates.
                  aBeam = abs(FuselagebeamX);
                  eBeam = (FuselagebeamX(1, :) - FuselagebeamX(1, 1)) ./ obj.Length;
                  aX1   = abs(horzcat(AeroPanels.X1));
                  aX4   = abs(horzcat(AeroPanels.X4));
                  
                  %Must create new coordinate systems for the splines
                  %   - These coordinate systems must have their y-axis along the
                  %     beam axis.
                  [cs, eta] = createCoordSysObjects(obj);
                  SplineCoordSys = makeFECoordSys(obj, cs(1 : end - 1));
                  coordSysIndex  = 1 : numel(cs);
                  
                  for jAP = 1:4
                      Num=0.25*numel(AeroPanels);
                      
                      for iAP = (jAP-1)*Num+1 : Num*jAP
                          
                          %Index the coordinates within X1/X4 of this aero panel set
                          idx = and(aBeam(1, :) >= aX1(1, iAP), aBeam(1, :) <= aX4(1, iAP));
                          
                          %Grab the nodes and add to the structural set
                          beamNode = FuselageBeamNodes(idx);
                                                   
                          if jAP==1
                              StrucSet(iAP).Nodes = horzcat(beamNode, Nodesright(idx))';
                          elseif jAP==2
                              StrucSet(iAP).Nodes = horzcat(beamNode, Nodesupper(idx))';
                              
                          elseif jAP==3
                              StrucSet(iAP).Nodes = horzcat(beamNode, Nodesleft(idx))';
                          elseif jAP==4
                              StrucSet(iAP).Nodes = horzcat(beamNode, Nodeslower(idx))';
                          end
                                                                       
                          %Select coordinate system which defines the beam axis
                          %   - Each CAERO1 segment must be planar, therefore, each
                          %     node in the 'beamNode' variable will have the same
                          %     coordinate system defining the orientation of the
                          %     beam. -> Can just reference the output coordinate
                          %     system of the node as each node has its output
                          %     coordinate system aligned with the beam axis.
                          index = interp1(eta, coordSysIndex, eBeam(idx), 'previous');
                          Spline(iAP).CoordSys = SplineCoordSys(index(1));
                          
                      end
                      
                  end
                  
                  %Populate the aerodynamic cards
                  %   - Each CAERO1 card has an AELIST with all the ID numbers
                  %   - Each spline points to its corresponding SET1 entry
                  %   - Each Splint points to its corresponding AELIST entry
                  %     i.e. 1 x spline per CAERO1
                  set(AeroSet, {'AeroPanels'}  , num2cell(AeroPanels)');
                  set(Spline , {'AeroPanelSet'}, num2cell(AeroSet)');
                  set(Spline , {'StrucSet'}    , num2cell(StrucSet)');
                  
                  
                  %% temperary remove aeropanels 
%                   addFEData(FEModel, AeroProp, ...
%                       AeroSet, StrucSet, Spline, SplineCoordSys);
                 
                
            end
            
            function AeroPanel = i_makeAeroPanels(obj)
                %i_makeAeroPanels Makes the aerodynamic panels for the
                %BluffBody.
                %   - The BluffBody object is represented by a cruciform of
                %   aerodynamic panels.
                %   - As the BluffBody can have a varying cross-section we
                %   need to use a combination of quadralateral and
                %   triangular aero panels to best approximate the shape of
                %   the body.
                %   - Currently, the cross-section of the fuselage is
                %   assumed to be circular. Therefore the cross-section is
                %   symmetric about 2 axes and we only need to calculate
                %   the panel coordinates in one quadrant before reflecting
                %   the coordinates in the other quadrants.
                
                %Make the quad/tri aerodynamic panels for one quadrant
                QuadPanel = i_makeQuadPanels(obj);
                %charles commented out triangular panel shape
%                 TriPanel  = i_makeTriPanels(obj);
                
                %Combine the quadralateral and triangular panels into a single
                %set
%                 AeroPanel = [QuadPanel, TriPanel];
                
                %charles commented out triangular 
                AeroPanel = [QuadPanel];
                
%                 %Assign EID numbers
%                 totalPanel       = cumsum(AeroPanel.NumPanels);
%                 pID              = [FEModel.ID.AeroPanels, totalPanel + FEModel.ID.AeroPanels];
%                 AeroPanel.EID = pID(1 : end - 1);
%                 
%                 %Check that the total number of panels does not exceed the FE
%                 %object increment
%                 if ~isempty(totalPanel) && totalPanel(end) > FEModel.CardID
%                     error(['The number of aerodynamic panels in the object ', ...
%                         '''%s'' exceeds the default setting for the number ', ...
%                         'of FE objects per type. Change the property ' , ...
%                         '''CardID'' of the object ''awi.fe.FEModel''.'], obj.NameAndType);
%                 end
                
            end
            
            % charles comment out triangular methods
%             function TriPanel = i_makeTriPanels(obj)
%                 %i_makeTriPanels Creates the triangular panel sets for the
%                 %aerodynamic mesh of the BluffBody.                
%  
%                 TriPanel = [];
%                 
%                 %Get coordinates of the stick in the global coordinate system
%                 Pos = [obj.XData ; obj.YData ; obj.ZData] + obj.AbsPosition';
%                 
%                 %Determine segment length in the chordwise direction
%                 c = diff(obj.Eta_ .* obj.Length);
%                 
%                 %Determine points along the BluffBody where there is an
%                 %increase/decrease in radius.
%                 %   - If dr ~= 0 then it will be necessary to generate a
%                 %   triangular panel section to approximate the shape of the
%                 %   bluff-body
%                 r       = obj.Radius_;
%                 dr      = diff(r);
%                 idx_tri = (dr ~= 0);
%                 nTriSet = nnz(idx_tri);
%                 nZTri   = zeros(1, nTriSet);
%                 
%                 %Bail out if we don't actually need to define any panels
%                 if nTriSet ==  0
%                     return
%                 end
%                 
%                 %Create variables that describe the radius ('spanwise' length)
%                 %and the coordinates at the start and end of each segment.
%                 %   - Need this for logical indexing!
%                 r1 = r(1 : end - 1);
%                 r2 = r(2 : end);
%                 Pos1 = Pos(:, 1 : end - 1);
%                 Pos2 = Pos(:, 2 : end);                                
%                 
%                 %Spanwise length of each triangular section
%                 s = abs(r1(idx_tri) - r2(idx_tri));
%                 
%                 %Calculate aero panel size
%                 %   - 'panelLength' is dimension in the global X direction
%                 %   - 'panelWidth' is dimension in the global spanwise axis
%                 if isempty(obj.NumAeroPanel)
%                     panelLength = repmat(obj.AeroPanelLength, size(s));
%                 else
%                     panelLength = c ./ obj.NumAeroPanel;
%                 end
%                 panelWidth = panelLength .* obj.AeroPanelAR;
%                 
%                 %TODO - Check skewness of the triangular panels
%                 
%                 %Generate the objects
%                 TriPanel = arrayfun(@(x) awi.fe.AeroPanel(), 1 : nTriSet, 'Unif', false);
%                 TriPanel = horzcat(TriPanel{:});
%                 
%                 %Populate data for the triangular panel sets
%                 x1 = Pos1(:, idx_tri) + [nZTri ; r1(idx_tri) ; nZTri];
%                 x4 = Pos2(:, idx_tri) + [nZTri ; r2(idx_tri) ; nZTri];
%                 set(TriPanel, {'X1'}    , num2cell(x1, 1)');
%                 set(TriPanel, {'X4'}    , num2cell(x4, 1)');
%                 set(TriPanel, {'CHORD'} , num2cell([c(idx_tri) ; nZTri], 1)');
%                 set(TriPanel, {'NSPAN'} , num2cell(ceil(s ./ panelWidth), 1)');
% %                 %charles 
% %                 set(TriPanel, {'NSPAN'} , num2cell(ceil(s ./ panelWidth(idx_tri)), 1)');
%                 set(TriPanel, {'NCHORD'}, num2cell(ceil(c(idx_tri) ./ panelLength), 1)');   
% %                 %charles 
% %                 set(TriPanel, {'NCHORD'}, num2cell(ceil(c(idx_tri) ./ panelLength(idx_tri)), 1)');
%             end
            
            function QuadPanel = i_makeQuadPanels(obj)
                %i_makeQuadPanels Creates the quadralateral panel sets for
                %the aerodynamic mesh of the BluffBody.
                    
                QuadPanel = [];
                
                %Get coordinates of the stick in the global coordinate system
                Pos = [obj.XData ; obj.YData ; obj.ZData] + obj.AbsPosition';
                
                %Determine segment length in the chordwise direction
                c = diff(obj.Eta_ .* obj.Length);
                
                %Determine points along the BluffBody where there is an
                %increase/decrease in radius.
                %   - If dr ~= 0 then it will be necessary to generate a
                %   triangular panel section to approximate the shape of the
                %   bluff-body
                r  = obj.Radius_;
                
                %Create variables that describe the radius ('spanwise' length)
                %and the coordinates at the start and end of each segment.
                %   - Need this for logical indexing!
                r1 = r(1 : end - 1);
                r2 = r(2 : end);
                Pos1 = Pos(:, 1 : end - 1);
                Pos2 = Pos(:, 2 : end);     %#ok<NASGU>
                
                %Determine points along the BluffBody where it is necessary
                %to generate a quadralateral set of panels.
                %   - If r1 == 0 || r2 == 0 then the aerodynamic panels can be
                %   represented by a triangular set only.
                idx_quad = ~or(r1 == 0, r2 == 0);
                nQuadSet = nnz(idx_quad);
                nZQuad   = zeros(1, nQuadSet);                
                                
                %Bail out if we don't actually need to define any panels
                if nQuadSet ==  0
                    return
                end
                
                %Spanwise length of each quadrilateral section
                s = r1(idx_quad);
                
                %Calculate aero panel size
                %   - 'panelLength' is dimension in the global X direction
                %   - 'panelWidth' is dimension in the global spanwise axis
                if isempty(obj.NumAeroPanel)
                    panelLength = repmat(obj.AeroPanelLength, size(s));
                else
                    panelLength = c ./ obj.NumAeroPanel;
                end
               
                %original code
%                 panelWidth = panelLength .* obj.AeroPanelAR;
                
                %charles redefine width
                panelWidth = s ./ 5;
                
                %Define offset values for the coordinates of point 4 in all
                %quadrants                
                dX4{1} = [nZQuad ; r1(idx_quad) ; nZQuad];  %right set
                dX4{2} = [nZQuad ; nZQuad ; r1(idx_quad)];  %upper set
                dX4{3} = [nZQuad ; -r1(idx_quad) ; nZQuad]; %left set
                dX4{4} = [nZQuad ; nZQuad ; -r1(idx_quad)]; %lower set
                
                %Add offsets to the coordinates of point 1 to obtain the
                %coordinates of point 4 in all quadrants
                X4 = cellfun(@(x) Pos1(:, idx_quad) + x, dX4, 'Unif', false);
                
                %Generate the objects
                %   - One set of panel objects for each quadrant
                QuadPanel = arrayfun(@(x) awi.fe.AeroPanel(), 1 : (4 * nQuadSet), 'Unif', false);
                QuadPanel = horzcat(QuadPanel{:});
                QuadPanel = reshape(QuadPanel, [4, nQuadSet]);
                
                %Populate data for the quadralateral panel sets 
                %   - TODO : Vectorise the SET for all objects
                x1 = num2cell(Pos1(:, idx_quad), 1)';
                c_ = num2cell(repmat(c(idx_quad), [2, 1]), 1)';
                ns = num2cell(ceil(s ./ panelWidth), 1)';
%                 %charles 
%                 ns = num2cell(ceil(s ./ panelWidth(idx_quad)), 1)';
                nc = num2cell(ceil(c(idx_quad) ./ panelLength)');
%                 %charles
%                 nc = num2cell(ceil(c(idx_quad) ./ panelLength(idx_quad))');
                
                for i = 1 : 4
                    set(QuadPanel(i, :), {'X1'}    , x1);
                    set(QuadPanel(i, :), {'X4'}    , num2cell(X4{i}, 1)');
                    set(QuadPanel(i, :), {'CHORD'} , c_);
                    set(QuadPanel(i, :), {'NSPAN'} , ns);
                    set(QuadPanel(i, :), {'NCHORD'}, nc);
                end

                %Collapse the [4,nP] matrix into a vector for storing
                QuadPanel = QuadPanel(:)';
                
                
            end
          
        end
    end
    
    methods % visualisation 
        function hg = drawElement(obj, ht)
            
            %Start with base-class
            hg = drawElement@awi.model.Beam(obj, ht);
            
        end        
    end
    
    methods % analytic methods
        function [cs, eta] = createCoordSysObjects(obj)
            %createCoordSysObjects Creates the 'awi.model.CoordSys' objects
            %along the span of the 'BluffBody'.
            %
            % As we assume the bluff-body is orientated along the global
            % x-direction we will assume the local x-axis goes along the
            % beam, the local-yaxis is aligned with the global y-axis and
            % the z-axis is normal to that plane.
            
            nDat = numel(obj.XData);
            dG2  = [zeros(1, nDat) ; ones(1, nDat) ; zeros(1, nDat)];
            
            %Use the superclass method
            [cs, eta] = createCoordSysObjects@awi.model.Beam(obj, dG2);
                        
        end
        function [Crs, eta] = createCrossSectionObjects(obj, varargin)
            %createCrossSectionObjects Generates the cross-section objects
            %that describe the aerofoil profile along the beam and assigns
            %them to the beam.
            %
            % The profiles can be assigned without the model being 'built'
            % but the orientation cannot be defined as it depends on the
            % beam object being fully populated.            
            
            [Crs, eta] = createCrossSectionObjects@awi.model.Beam(obj);     
            if isempty(eta)
                eta = obj.Eta;
            end
            nCrs       = max([numel(Crs), numel(eta)]);
            
            r = num2cell(2 * obj.Radius)';
            if isempty(r)
                return
            end            
            
            %Do we need to make new objects?
            if isempty(Crs) || numel(Crs) ~= numel(r)
                %Make the 'awi.model.CrossSection' objects and assign the
                %'awi.model.CrossSection' objects to the object
                Crs  = arrayfun(@(~) awi.model.CrossSection, 1 : nCrs);                
                eta  = obj.Eta;
                assignBeamObject(obj, Crs, eta, 'replace');                
            end
            
            eta  = [Crs.BeamEta];            
            
            %Generate the 'awi.model.CrossSection' objects and assign them to the beam
            set(Crs, 'CrossSectionName', 'Ellipse');
            set(Crs, {'MajorAxis'}, r);
            set(Crs, {'MinorAxis'}, r);
            generateCoordsFromLibrary(Crs);            
            
        end
    end
    
end
