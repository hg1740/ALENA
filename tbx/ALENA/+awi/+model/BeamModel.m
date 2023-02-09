classdef (ConstructOnLoad) BeamModel < awi.model.Entity %awi.model.ResultSet
    %BeamModel Describes a mass and stiffness distribution for a given 
    %aircraft.
    
    %AWI-specific properties
    properties    
        %Extend metadata with content specific to awi
        Aircraft; 
        LoadCase;
        BM;        
    end
    
    methods % set / get
       function set.Aircraft(obj, val)    %set.Aircraft    
%             validateattributes(val, {'awi.model.Aircraft'}, {'scalar'}, ...
%                 class(obj), 'Aircraft');
            validateattributes(val, {'char'}, {'row'}, class(obj), ...
                'Aircraft');
            obj.Aircraft = val;
        end 
    end
     
    methods % construction / destruction
        
        function obj = BeamModel(varargin)
            %BeamModel Constructor for the 'BeamModel' class.
            
            %Adjust properties to allow the Beam Model to be an intermediate node in result sets
            args = {'IsLeafNode', false, 'StructureLocked', false};
            
            %Pass it on
            obj@awi.model.Entity(args{:}, varargin{:});
            
            %Configure collectabls
            obj.addCollectionSpec(@awi.model.ResultSet         , [], [], [], true, [], 'Static Results');    %Hidden, grouped
            obj.addCollectionSpec(@awi.model.TransientResultSet, [], [], [], true, [], 'Transient Results'); %Hidden, grouped
            obj.addCollectionSpec(@awi.model.Collector, 'Collector');            
            
            %Extend property groups
            obj.addPropertyGroup('General', ...
                'LoadCase', 'LoadCase'    , ...
                'Aircraft', 'Aircraft');
            
        end
        
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ha, tag)
            
            hg = [];

            if nargin < 3
                tag = obj.Name;
                if isempty(tag)
                   tag = 'Beam Model'; 
                end
            end
            if isempty(obj.BM)
                return
            end
            
            if isa(obj.BM, 'awi.fe.FEModel')
                hg = draw(obj.BM, ha);
                return
            end
            %TODO - Update the drawing the beam model so that it is more
            %robust. Currently, when a model is implemented we are calling
            %draw after each primary node (e.g. Aircraft, Load Cases) has
            %been imported. If I have a FAME beam model (FAME FEM) in the
            %model it does not have the structure 'BM.Aero.lattice' so it
            %cannot be drawn and subsequently errors.
            if ~isfield(obj.BM.Aero, 'lattice')
                return
            end
            
            % Plot the aerodynamic lattice
            [npwing, ~, ~] = size(obj.BM.Aero.lattice.XYZ);
            
            % Setup handle groups for specific properties of the plot
            h_lattice  = hggroup(ha);
            h_beam     = hggroup(ha);
            h_bnodes   = hggroup(ha);
            h_anodes   = hggroup(ha);
            h_conm2    = hggroup(ha);
            h_conm2conn= hggroup(ha);
            h_body     = hggroup(ha);
            
            % Plot aerodynamic mesh
            hg = plot3(obj.BM.Aero.lattice.XYZ(1:npwing,:,1)',...
                obj.BM.Aero.lattice.XYZ(1:npwing,:,2)',...
                obj.BM.Aero.lattice.XYZ(1:npwing,:,3)','-bo',...
                'MarkerSize', 1, 'MarkerFaceColor','b','LineWidth',2,'Parent',h_lattice);
            
            %hold(ha,'on');
            
            % Plot the beam connections
            BeamNodes = [];
            for i = 1: length(obj.BM.Beam.ID)
                [~,idx] = intersect(obj.BM.Node.ID,obj.BM.Beam.Conn(i,:)'); 
                BeamNodes = [BeamNodes;idx];
                hg(end+1) = plot3(obj.BM.Node.Coord(idx,1),obj.BM.Node.Coord(idx,2),...
                    obj.BM.Node.Coord(idx,3),'k-','LineWidth',1.5,'Parent',h_beam);
                %hold(ha,'on');
            end
            
            % Plot Beam nodes
            BeamNodes = sort(unique(BeamNodes));
            %hold(ha,'on');
            hg(end+1) = plot3(obj.BM.Node.Coord(BeamNodes,1),...
                obj.BM.Node.Coord(BeamNodes,2),...
                obj.BM.Node.Coord(BeamNodes,3),...                
                'MarkerFaceColor','b',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_bnodes);
            
            % Plot the aerodynamic nodes
            AeroNodes = setdiff((1:length(obj.BM.Node.ID))',BeamNodes);
            %hold(ha,'on');
            hg(end+1) = plot3(obj.BM.Node.Coord(AeroNodes,1),...
                obj.BM.Node.Coord(AeroNodes,2),...
                obj.BM.Node.Coord(AeroNodes,3),...
                'MarkerFaceColor','r',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_anodes);


            % Plot the lumped masses
            
            nodeidx = arrayfun(@(x) find(obj.BM.Node.ID == x), obj.BM.Conm2.Node);
            hg(end+1) = plot3(obj.BM.Node.Coord(nodeidx,1) + obj.BM.Conm2.Offset(:,1),...
                obj.BM.Node.Coord(nodeidx,2) + obj.BM.Conm2.Offset(:,2),...
                obj.BM.Node.Coord(nodeidx,3) + obj.BM.Conm2.Offset(:,3),...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_conm2);
            
            for i = 1:numel(obj.BM.Conm2.Node)
                hg(end+1) = plot3([obj.BM.Node.Coord(nodeidx(i),1), obj.BM.Node.Coord(nodeidx(i),1)+ obj.BM.Conm2.Offset(i,1)],...
                    [obj.BM.Node.Coord(nodeidx(i),2), obj.BM.Node.Coord(nodeidx(i),2)+ obj.BM.Conm2.Offset(i,2)],...
                    [obj.BM.Node.Coord(nodeidx(i),3), obj.BM.Node.Coord(nodeidx(i),3)+ obj.BM.Conm2.Offset(i,3)],'g-','Parent',h_conm2conn);
            end
            
            % Plot the bluff_bodies            
            for i = 1:numel(obj.BM.Aero.body.ID)
                BodyColour = obj.BM.Aero.body.Colour(i,:);
                hg(end+1) = patch('Faces',obj.BM.Aero.body.lattice.Elem.Conn{i},'Vertices',obj.BM.Aero.body.lattice.Elem.Node{i},'FaceColor',BodyColour,'EdgeColor','k','FaceAlpha',0.7,'Parent',h_body);
            end
           
            
            %legend(ha,[h_lattice,h_beam,h_bnodes,h_anodes,h_conm2,h_conm2conn,h_body],{'Aero Mesh','Beam Elements','Beam Nodes','Aero Nodes','Conm2','Conm2 Connection','Bluff Body'});
            
            set(hg,'Tag',tag);
            
            function hg = i_drawMSC(BeamModel, ha)
                
                %Setup handle groups for specific properties of the plot
                hNodes = hggroup(ha);                
                hBeams = hggroup(ha);   
                hRigid = hggroup(ha);
                hAero  = hggroup(ha);
                
                emptyFlag = ~[BeamModel.HasFEData];
                if nnz(emptyFlag) == numel(BeamModel)
                    hg = gobjects(1);
                    return
                end
                
                %% Draw Nodes 
                
                %Collect everything
                Nodes = [BeamModel.Nodes];
                GID   = [Nodes.GID];
                X     = [Nodes.X];
                                
                hg = plot3(X(1, :), X(2, :), X(3, :), ...
                    'LineStyle'      , 'none', ...
                    'Marker'         , 'o', ...
                    'MarkerFaceColor', 'g', ...
                    'MarkerEdgeColor', 'k', ...
                    'Tag'            , 'Nodes', ...
                    'Parent'         , hNodes);
                
                %% Draw Beams 
                
                %Collect everything
                Beams = [BeamModel.Beams];                                
                GA    = [Beams.GA];
                GB    = [Beams.GB]; 
                
                %Index the nodes
                [~, ind_A] = ismember(GA, GID); 
                [~, ind_B] = ismember(GB, GID);
                X_A = X(:, ind_A);
                X_B = X(:, ind_B);
                
                %Grab coordinates
                xb = [X_A(1, :) ; X_B(1, :)];
                yb = [X_A(2, :) ; X_B(2, :)];
                zb = [X_A(3, :) ; X_B(3, :)];
                
                %Plot
                hL = plot3(xb,yb,zb, 'k-', ...
                    'Parent', hBeams, ...
                    'Tag'   , 'Beam Elements');
                
                %% Draw Rigid Bars
                
                %Collect everything
                RigidBars = [BeamModel.RigidBars];                
                GN = arrayfun(@(s) repmat(s.GN, size(s.GM)), RigidBars, 'Unif', false);
                GN = vertcat(GN{:});
                GM = vertcat(RigidBars.GM);
                
                %Index the nodes
                [~, ind_N] = ismember(GN, GID); 
                [~, ind_M] = ismember(GM, GID);
                X_N = X(:, ind_N);
                X_M = X(:, ind_M);
                
                %Grab coordinates
                xr = [X_N(1, :) ; X_M(1, :)];
                yr = [X_N(2, :) ; X_M(2, :)];
                zr = [X_N(3, :) ; X_M(3, :)];
                
                %Plot
                hR = plot3(xr, yr, zr, ...
                    'Color'    , 'r', ...
                    'LineStyle', '-', ...
                    'LineWidth', 1  , ...
                    'Marker'   , 'o', ...
                    'Tag'      , 'Rigid Elements', ...
                    'Parent'   , hRigid);                
                
                %% Draw aerodynamic panels                
                AeroPanel = definePanels([BeamModel.AeroPanels]);
                
                hA = fill3( ...
                    AeroPanel.coord(:, :, 1), ...
                    AeroPanel.coord(:, :, 2), ...
                    AeroPanel.coord(:, :, 3), ...
                    'k'     , 'FaceColor', 'none' , ...
                    'Tag'   , 'Aerodynamic Panels', ...
                    'Parent', hAero);
                
                %% Combine graphics objects
                hg = [hg ; hL ; hR ; hA];
                
            end
        end
        
    end
    
end