classdef SizeAnalysis < awi.view.Analysis
    
    properties
        
        Linear          = true;  % true or false to invoke linear analysis
        NonLinear       = false; % true or false to invoke nonlinear analysis
        Elements        = 25;    % number of beam elements
        TrimResults     = {};
        OptimVariables  = {};
        BeamModel       = {};
        AnalysisType    = {};
        AeroOptimisation   = false;
        StructOptimisation = true;
        AeroDistribution   = 'Elliptical';
        
    end
    
    properties (Dependent)
        
        CruiseLoadCase;
        
    end
    
    properties (Transient, SetAccess = protected )
        
        hLinear;
        hNonLinear;
        hElements;
        hAeroOptim;
        hStructOptim;
        hAeroType;
        
    end
    
    methods % get/set
        
        function val = get.CruiseLoadCase(obj)
            
            %Start with all LoadCases
            val = obj.LoadCases;
            
            %Anything ?
            if ~isempty(val)
                
                %Down-select
                idx = find(val, 'LoadCaseType', 'Cruise Design Case', '-index');
                
                %Nothing ?
                if isempty(idx)
                    
                    %A Cruise loadcase is just a manouevre with a 1g LoadFactor
                    idx = find(val, 'LoadCaseType', 'manoeuvre', 'LoadFactor', 1, '-index');
                    
                end
                
                %Send back whatever we've got
                val = val(idx);
                
            end
            
        end

    end
    
    methods
        
        function obj = SizeAnalysis( model, varargin )
            
            %Call superclass constructor - do not pass the model in, better to wait
            obj@awi.view.Analysis( model, ...                'Title', 'Size Analysis', ...
                'OfferLoadCaseSelection', 'none', ...
                'PushButtons', {'Analyse...', @analyse}, ...
                varargin{:} );
            
            % Allow the user to select linear anlysis
            obj.hLinear = uicontrol('Parent',obj.hControls,...
                'Style','checkbox',...
                'String','Linear',...
                'Value',double(obj.Linear),'Callback',@obj.cbLinear);
            
            %Good height
            obj.hControls.Heights(end) = 25;
            
            % Allow the user to select nonlinear            % anlysis
            obj.hNonLinear = uicontrol('Parent',obj.hControls,...
                'Style','checkbox',...
                'String','Nonlinear',...
                'Value',double(obj.NonLinear),'Callback',@obj.cbNonLinear);
            
            %Good height
            obj.hControls.Heights(end) = 25;
            
            % Allow the user to change the number of beam elements
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', '# beam elements:','HorizontalAlignment','Left');
            
            obj.hElements = uicontrol('Parent',obj.hControls,...
                'Style','edit',...
                'String',num2str(obj.Elements),'Callback',@obj.cbElements);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;
            
            % Allow the user to select linear anlysis
            obj.hAeroOptim = uicontrol('Parent',obj.hControls,...
                'Style','checkbox',...
                'String','Optimise Jig Twist',...
                'Value',double(obj.AeroOptimisation),...
                'Callback',@obj.cbAeroOptim);
            
            %Good height
            obj.hControls.Heights(end) = 25;
            
            % Allow the user to change the number of beam elements
            obj.hAeroType = uicontrol('Parent', obj.hControls,...
                'Style', 'popup',...
                'String', {'TargetCl','Elliptical','Triangular','Bell-Shaped'},...
                'Callback',@obj.cbAeroType,...
                'Enable','off');
            
            %Good height
            obj.hControls.Heights(end) = 25;
            
            % Allow the user to select linear anlysis
            obj.hStructOptim = uicontrol('Parent',obj.hControls,...
                'Style','checkbox',...
                'String','Optimise Wing Box',...
                'Value',double(obj.StructOptimisation),'Callback',@obj.cbStructOptim);
            
            %Good height
            obj.hControls.Heights(end) = 25;
            
            %Now store the model - which will in turn trigger an update
            obj.Model = model;
            
        end
        
    end
    
    methods
        
        function R = analyse(obj, varargin)
            
            %Have we got enough to work on ?
            assert(numel(obj.LoadCases) > 1, 'Size analysis requires at least TWO loadcases, one of which must be a CRUISE case');
            assert(numel(obj.CruiseLoadCase) == 1, 'Size analysis requires one CRUISE loadcase');
            
            % Delete the current results view
            delete(obj.hPanel.Children);
            
            obj.OptimVariables  = {};
            obj.TrimResults     = {};
            obj.BeamModel       = {};
            obj.AnalysisType    = {};
            
            %Something to look at during analysis
            logfcn = progressdlg(obj.Model, 'Size analysis...');
            logfcn('Performing size analysis...');
            
            %Which interface to custom code ?
            if ~isa(obj, 'mvc.mixin.PathAdjustable') || ...
                    obj.Model.Root.CustomCodeFolderIndex == 1
                
                %Find all lifting surfaces and bluff bodies
                LSObj = findall(obj.Aircraft, 'Type', 'LiftingSurface');
                BBObj = findall(obj.Aircraft, 'Type', 'BluffBody');
                
                %BUT convertFmwk2Neo ASSUMES that all lifting surfaces have spars and control surfaces,
                % so to avoid downstream error, identify and eliminate any that do not
                %TODO: Check with Robbie what is right thing to do here
                LSObj(arrayfun(@(x)isempty(findall(x, 'Type', 'Spar')), LSObj)) = [];
                LSObj(arrayfun(@(x)isempty(x.ControlSurfaces), LSObj)) = [];
                
                % Convert to a usable Aeroflex workspace
                [beam_model,Aircraft_param,Optim_Variables,LoadCases] = ...
                    convertFmwk2Neo(LSObj,BBObj,obj.Aircraft,obj.LoadCases,obj.Elements, logfcn);
            
            else
                
                %OLD CODE INTERFACE
                
                % Find the starboard wing
                stbdWing = findall(obj.Aircraft,'Name','StbdWing');
                
                % Run an assertion to check its existence
                assert(~isempty(stbdWing),'Starboard wing not found');
                
                % Convert to a usable Aeroflex workspace
                [beam_model,Aircraft_param,Optim_Variables,LoadCases] = ...
                    convertFmwk2Neo(stbdWing,obj.Aircraft,obj.LoadCases,obj.Elements, logfcn);
                
            end
            
            % Restructure the variables so that they can be used by the
            % AircraftSizer
            GUI_Input                = obj.InitialiseGUI;
            GUI_Input.LoadCases      = LoadCases;
            GUI_Input.beam_model     = beam_model;
            GUI_Input.Aircraft_param = Aircraft_param;
            
            GUI_Input.Optim_Variables = Optim_Variables;
            
            % Initialise the panel for the view
            hp = uiextras.TabPanel('Parent', obj.hPanel, 'TabWidth', 100);
            hb = uiextras.HBox('Parent', hp);
            hp.TabTitles{end} = 'Progress';
            
            % Update the Aeroflex variable according to the
            % linear/nonlinear choice
            
            GUI_Input.Parameters.Struc_optimisation = obj.StructOptimisation;
            GUI_Input.Parameters.Aero_optimisation  = obj.AeroOptimisation;
            
            GUI_Input.Parameters.AeroDistribution = obj.AeroDistribution;
            
            if obj.Linear
                
                GUI_Input.Parameters.Linear = true;
                
                obj.AnalysisType{end+1} = 'Lin';
                
                % Run the sizing
                [obj.TrimResults{end+1},obj.OptimVariables{end+1},AeroFlexModel] = AircraftSizer(GUI_Input,hb, logfcn);
                
                % Pass the beam model here to extract what is necessary for
                % the framework
                obj.BeamModel{end + 1} = ExtractBareBeamModel(obj,AeroFlexModel);
            end
            
            if obj.NonLinear
                
                GUI_Input.Parameters.Linear = false;
                
                obj.AnalysisType{end+1} = 'NonLin';
                
                % Run the sizing
                [obj.TrimResults{end+1},obj.OptimVariables{end+1},AeroFlexModel] = AircraftSizer(GUI_Input,hb, logfcn);
                
                % Pass the beam model here to extract what is necessary for
                % the framework
                obj.BeamModel{end + 1} = ExtractBareBeamModel(obj,AeroFlexModel);
                
            end
            
            
            % Plot some important results immediately after the size
            % analysis
            
            % Plot the internal loads
            plot_intload_loadcase(obj,hp);
            hp.TabTitles{end} = 'Internal Loads';
            hp.Selection = numel(hp.Children);
            
            % Plot the discretised beam model
            plot_beam_model(obj,hp);
            hp.TabTitles{end} = 'Beam Model';
            hp.Selection = numel(hp.Children);
            
            for i = 1:numel(obj.BeamModel)
                
                %Insert an intermediate Beam Model
                B = awi.model.BeamModel( ...
                    'Aircraft', obj.Aircraft.Name, ...
                    'LoadCase', {obj.LoadCases.Name},...
                    'BM', obj.BeamModel{i},...
                    'Name','AeroFlex FEM');
                
                %Add to session
                obj.Model.add(B);
                
                % What is the number of load cases
                sizeParts = fieldnames(obj.OptimVariables{i});
                nResults = size(obj.OptimVariables{i}.(sizeParts{1}).Fx,2);
                
                % Construct StbdWingResults
                StbdWingResults = arrayfun(@(i) awi.model.AeroelasticResult, 1 : nResults, 'Unif', false);
                StbdWingResults = horzcat(StbdWingResults{:});
                
                % Construct Manoeuver results
                ManResults = arrayfun(@(i) awi.model.ResultSet, 1 : nResults, 'Unif', false);
                ManResults = horzcat(ManResults{:});

                % Assign some meaningful labels
                Labels = arrayfun(@(i) sprintf('AeroFlex Result (LC %i)', i), 1:nResults , 'Unif', false);
                set(ManResults, {'Name'}, Labels');
                
                StbdWing = findall(obj.Aircraft, 'Name', 'StbdWing');
                set(StbdWingResults, 'Beam', StbdWing);
                
                % Format the results from AEROFLEX
                %   - All loads have the same eta positions
                eta = arrayfun(@(x) obj.OptimVariables{i}.Wing.y_mbox/Aircraft_param.Wing.Planform.b_ref,1:nResults,'Unif',false)';
                Fx  = arrayfun(@(x) obj.OptimVariables{i}.Wing.Fx(:,x),1:nResults,'Unif',false)';
                Fy  = arrayfun(@(x) obj.OptimVariables{i}.Wing.Fy(:,x),1:nResults,'Unif',false)';
                Fz  = arrayfun(@(x) obj.OptimVariables{i}.Wing.Fz(:,x),1:nResults,'Unif',false)';
                Mx  = arrayfun(@(x) obj.OptimVariables{i}.Wing.Mx(:,x),1:nResults,'Unif',false)';
                My  = arrayfun(@(x) obj.OptimVariables{i}.Wing.My(:,x),1:nResults,'Unif',false)';
                Mz  = arrayfun(@(x) obj.OptimVariables{i}.Wing.Mz(:,x),1:nResults,'Unif',false)';
                
                partidx = and(ismember({obj.BeamModel{1}.PartId.Part},'Wing_R'),ismember({obj.BeamModel{1}.PartId.Type},'GRID')) ;
                
                nodeidx = obj.BeamModel{1}.PartId(partidx).index;
                eta_disp =  arrayfun(@(x) AeroFlexModel.Node.Coord(nodeidx,2)/max(AeroFlexModel.Node.Coord(nodeidx,2)),1:nResults,'Unif',false)';            
                Tx  = arrayfun(@(x) obj.TrimResults{i}.Displacements{x}(nodeidx,1),1:numel(obj.TrimResults{i}.Displacements),'Unif',false)';
                Ty  = arrayfun(@(x) obj.TrimResults{i}.Displacements{x}(nodeidx,2),1:numel(obj.TrimResults{i}.Displacements),'Unif',false)';
                Tz  = arrayfun(@(x) obj.TrimResults{i}.Displacements{x}(nodeidx,3),1:numel(obj.TrimResults{i}.Displacements),'Unif',false)';
                
                %Assign to the objects
                set(StbdWingResults, {'Fx_eta'}, eta);
                set(StbdWingResults, {'Fx'}    , Fx);  
                set(StbdWingResults, {'Fy_eta'}, eta);
                set(StbdWingResults, {'Fy'}    , Fy);  
                set(StbdWingResults, {'Fz_eta'}, eta);
                set(StbdWingResults, {'Fz'}    , Fz);
                set(StbdWingResults, {'Mx_eta'}, eta);
                set(StbdWingResults, {'Mx'}    , Mx);
                set(StbdWingResults, {'My_eta'}, eta);
                set(StbdWingResults, {'My'}    , My);
                set(StbdWingResults, {'Mz_eta'}, eta);
                set(StbdWingResults, {'Mz'}    , Mz);
                
                                %Assign to the objects
                set(StbdWingResults, {'FxFAME_eta'}, eta);
                set(StbdWingResults, {'FxFAME'}    , Fz);  
                set(StbdWingResults, {'FyFAME_eta'}, eta);
                set(StbdWingResults, {'FyFAME'}    , Fx);  
                set(StbdWingResults, {'FzFAME_eta'}, eta);
                set(StbdWingResults, {'FzFAME'}    , Fy);
                set(StbdWingResults, {'MxFAME_eta'}, eta);
                set(StbdWingResults, {'MxFAME'}    , Mz);
                set(StbdWingResults, {'MyFAME_eta'}, eta);
                set(StbdWingResults, {'MyFAME'}    , Mx);
                set(StbdWingResults, {'MzFAME_eta'}, eta);
                set(StbdWingResults, {'MzFAME'}    , My);
                
                set(StbdWingResults, {'Tx'}    , Tx);
                set(StbdWingResults, {'Tx_eta'}, eta_disp);
                set(StbdWingResults, {'Ty'}    , Ty);
                set(StbdWingResults, {'Ty_eta'}, eta_disp);
                set(StbdWingResults, {'Tz'}    , Tz);
                set(StbdWingResults, {'Tz_eta'}, eta_disp);
                
                partidx  = and(ismember({obj.BeamModel{1}.PartId.Part},'Wing_R'),ismember({obj.BeamModel{1}.PartId.Type},'CAERO')) ;
                aeroidx  = obj.BeamModel{1}.PartId(partidx).index;
                
                Cl        = arrayfun(@(x) obj.TrimResults{i}.NonDimAeroForces{x}(aeroidx,3),1:numel(obj.TrimResults{i}.NonDimAeroForces),'Unif',false)';
                eta_aero =  arrayfun(@(x) obj.TrimResults{i}.p_mid_r{1}(aeroidx)/max(AeroFlexModel.Node.Coord(nodeidx,2)),1:nResults,'Unif',false)';            
                
                set(StbdWingResults, {'Cl'}      , Cl);
                set(StbdWingResults, {'Cl_eta'}  , eta_aero);
               
                BR = num2cell(StbdWingResults);
                set(ManResults, {'BeamResults'}, BR');
                
                B.add(ManResults);
                
            end
            
            %Retain log for audit trail, if possible
            if isa(obj.Model, 'mvc.mixin.Auditable')
                
                %What happened ?
                [~, str] = logfcn();
                
                %Yes
                obj.Model.addAuditTrailEntry(str);
                
                %Done with progress view
                logfcn('close');
               
            end
            
        end
        
    end
    
    methods ( Access = protected )
        
        function val = getLoadCases(obj)
            
            %For size analysis, we do not just return all load cases in session,
            % although it is a good place to start
            val = getLoadCases@awi.view.Analysis(obj);
            
            %Anything ?
            if ~isempty(val)
                
                % Find the manouevre load cases:
                ManObj = findobj(val, 'LoadCaseType', 'Manoeuvre');
                
                % Find the cruise case:
                CruiseObj = findobj(val, 'LoadCaseType', 'Cruise Design Case');
                
                %NO too harsh to error at this point - wait till later
                % Run assert on the cruise object.
                % assert((numel(CruiseObj) == 1),'A single cruise case must be defined.')
                
                % Concatenate the loadcase objects
                val = cat(1,ManObj,CruiseObj);
            
            end
            
        end
        
        function update(obj)
            
            %Start with base class
            update@awi.view.Analysis(obj);
            
            %TODO: Display the result of the analysis (if it exists)
            
        end
        
        function cbElements(obj,~,~)
            
            obj.Elements = str2double(get(obj.hElements, 'String'));
            
        end
        
        function cbLinear(obj,~,~)
            
            obj.Linear = get(obj.hLinear,'Value')  == 1;
            
        end
        
        function cbNonLinear(obj,~,~)
            
            obj.NonLinear = get(obj.hNonLinear,'Value')  == 1;
            
        end
        
        function cbAeroOptim(obj,~,~)
            
            obj.AeroOptimisation = get(obj.hAeroOptim,'Value')  == 1;
            
            if obj.AeroOptimisation
                obj.hAeroType.Enable = 'on';
            else
                obj.hAeroType.Enable = 'off';
            end
            
        end
        
        function cbStructOptim(obj,~,~)
            
            obj.StructOptimisation = get(obj.hStructOptim,'Value')  == 1;
            
        end
        
                
        function cbAeroType(obj,~,~)
            
            AeroTypes = {'TargetCl','Elliptical','Triangular','Bell-Shaped'};
            
            obj.AeroDistribution = AeroTypes{get(obj.hAeroType,'Value')};
            
        end
        
                        
        function plot_intload_loadcase(obj,hp)
            
            hg = uiextras.Grid('Parent',hp);
            
            for i = 1:6
                hc = uicontainer('Parent',hg);
                ha(i) = axes('Parent',hc);
            end
            
            hg.Widths = [-1,-1,-1];
            
            legendstring = {};
            
            % Generate a list of line coloUrs
            col  = get(ha(1),'ColorOrder');
            
            if numel(obj.LoadCases) > size(col,1)
                col = hsv(numel(obj.LoadCases));
            end
            
            % Generate a list of linestyles to be used by the plot function
            hl = line(ha(1),nan,nan);
            sty  = set(hl,'Linestyle');%{'-','-.','--',':'};
            delete(hl);
            
            
            isty = 0;
            
            for jj = 1:numel(obj.OptimVariables)
                
                isty = isty + 1;
                if isty > numel(sty)
                    isty = 1;
                end
                
                Optim = obj.OptimVariables{jj};
                
                for ii = 1:size(Optim.Wing.Fx,2)
                    legendstring{end+1} = ['LC: ' obj.LoadCases(ii).Name ' (' obj.AnalysisType{jj} ')'];
                end
                
                icol = 0;
                
                for ii = 1:size(Optim.Wing.Fx,2)
                    
                    icol = icol + 1;
                    
                    if icol > size(col,1)
                        icol = 1;
                    end
                    
                    args = {'Color',col(icol,:),'Linestyle',sty{isty}};
                    
                    hl(1) = plot(ha(1),Optim.Wing.y_mbox,Optim.Wing.Fx(:,ii),args{:});
                    hold(ha(1),'on');
                    xlabel(ha(1),'y [m]');ylabel(ha(1),'F_x');
                    
                    hl(3) = plot(ha(3),Optim.Wing.y_mbox,Optim.Wing.Fy(:,ii),args{:});
                    hold(ha(3),'on');
                    xlabel(ha(3),'y [m]');ylabel(ha(3),'F_y');
                    
                    hl(5) = plot(ha(5),Optim.Wing.y_mbox,Optim.Wing.Fz(:,ii),args{:});
                    hold(ha(5),'on');
                    xlabel(ha(5),'y [m]');ylabel(ha(5),'F_z');
                    
                    hl(2) = plot(ha(2),Optim.Wing.y_mbox,Optim.Wing.Mx(:,ii),args{:});
                    hold(ha(2),'on');
                    xlabel(ha(2),'y [m]');ylabel(ha(2),'M_x');
                    
                    hl(4) = plot(ha(4),Optim.Wing.y_mbox,Optim.Wing.My(:,ii),args{:});
                    hold(ha(4),'on');
                    xlabel(ha(4),'y [m]');ylabel(ha(4),'M_y');
                    
                    hl(6) = plot(ha(6),Optim.Wing.y_mbox,Optim.Wing.Mz(:,ii),args{:});
                    hold(ha(6),'on');
                    xlabel(ha(6),'y [m]');ylabel(ha(6),'M_z');
                    
                end
                
            end
            
            legend(ha(1),legendstring,'Location','Best');
            
        end
        
        function plot_beam_model(obj,hp)
            
            hg = uiextras.Grid('Parent',hp);
            
            hc = uicontainer('Parent',hg);
            ha = axes('Parent',hc);
            
            % Plot the aerodynamic lattice
            [npwing, ~, ~] = size(obj.BeamModel{end}.Aero.lattice.XYZ);
            
            % Setup handle groups for specific properties of the plot
            h_lattice  = hggroup(ha);
            h_beam     = hggroup(ha);
            h_bnodes   = hggroup(ha);
            h_anodes   = hggroup(ha);
            h_conm2    = hggroup(ha);
            h_conm2conn= hggroup(ha);
            h_body     = hggroup(ha);
            
            % Plot aerodynamic mesh
            plot3(obj.BeamModel{end}.Aero.lattice.XYZ(1:npwing,:,1)',...
                obj.BeamModel{end}.Aero.lattice.XYZ(1:npwing,:,2)',...
                obj.BeamModel{end}.Aero.lattice.XYZ(1:npwing,:,3)','-bo',...
                'MarkerSize', 1, 'MarkerFaceColor','b','LineWidth',2,'Parent',h_lattice);
            
            hold(ha,'on');
            
            % Plot the beam connections
            BeamNodes = [];
            for i = 1: length(obj.BeamModel{end}.Beam.ID)
                [~,idx] = intersect(obj.BeamModel{end}.Node.ID,obj.BeamModel{end}.Beam.Conn(i,:)'); 
                BeamNodes = [BeamNodes;idx];
                plot3(obj.BeamModel{end}.Node.Coord(idx,1),obj.BeamModel{end}.Node.Coord(idx,2),...
                    obj.BeamModel{end}.Node.Coord(idx,3),'k-','LineWidth',1.5,'Parent',h_beam);
                hold(ha,'on');
            end
            
            % Plot Beam nodes
            BeamNodes = sort(unique(BeamNodes));
            hold(ha,'on');
            plot3(obj.BeamModel{end}.Node.Coord(BeamNodes,1),...
                obj.BeamModel{end}.Node.Coord(BeamNodes,2),...
                obj.BeamModel{end}.Node.Coord(BeamNodes,3),...
                'MarkerFaceColor','b',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_bnodes);
            
            % Plot the aerodynamic nodes
            AeroNodes = setdiff((1:length(obj.BeamModel{end}.Node.ID))',BeamNodes);
            hold(ha,'on');
            plot3(obj.BeamModel{end}.Node.Coord(AeroNodes,1),...
                obj.BeamModel{end}.Node.Coord(AeroNodes,2),...
                obj.BeamModel{end}.Node.Coord(AeroNodes,3),...
                'MarkerFaceColor','r',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_anodes);

            % Plot the lumped masses
            
            nodeidx = arrayfun(@(x) find(obj.BeamModel{end}.Node.ID == x), obj.BeamModel{end}.Conm2.Node);
            plot3(obj.BeamModel{end}.Node.Coord(nodeidx,1) + obj.BeamModel{end}.Conm2.Offset(:,1),...
                obj.BeamModel{end}.Node.Coord(nodeidx,2) + obj.BeamModel{end}.Conm2.Offset(:,2),...
                obj.BeamModel{end}.Node.Coord(nodeidx,3) + obj.BeamModel{end}.Conm2.Offset(:,3),...
                'MarkerFaceColor','g',...
                'MarkerEdgeColor','k',...
                'Marker','o',...
                'LineStyle','none',...
                'Parent',h_conm2);
            
            for i = 1:numel(obj.BeamModel{end}.Conm2.Node)
                plot3([obj.BeamModel{end}.Node.Coord(nodeidx(i),1), obj.BeamModel{end}.Node.Coord(nodeidx(i),1)+ obj.BeamModel{end}.Conm2.Offset(i,1)],...
                    [obj.BeamModel{end}.Node.Coord(nodeidx(i),2), obj.BeamModel{end}.Node.Coord(nodeidx(i),2)+ obj.BeamModel{end}.Conm2.Offset(i,2)],...
                    [obj.BeamModel{end}.Node.Coord(nodeidx(i),3), obj.BeamModel{end}.Node.Coord(nodeidx(i),3)+ obj.BeamModel{end}.Conm2.Offset(i,3)],'g-','Parent',h_conm2conn);
            end
            
            % Plot the bluff_bodies            
            for i = 1:numel(obj.BeamModel{end}.Aero.body.ID)
                BodyColour = obj.BeamModel{end}.Aero.body.Colour(i,:);
                patch('Faces',obj.BeamModel{end}.Aero.body.lattice.Elem.Conn{i},'Vertices',obj.BeamModel{end}.Aero.body.lattice.Elem.Node{i},'FaceColor',BodyColour,'EdgeColor','k','FaceAlpha',0.7,'Parent',h_body);
            end
            
            axis(ha,'tight');
            axis(ha,'equal');
            grid(ha,'on');
            view(ha,[-45,45]);
            
            legend(ha,[h_lattice,h_beam,h_bnodes,h_anodes,h_conm2,h_conm2conn,h_body],{'Aero Mesh','Beam Elements','Beam Nodes','Aero Nodes','Conm2','Conm2 Connection','Bluff Body'});
            
        end
        
        function Output = ExtractBareBeamModel(obj,Input)
            
            %%ExtractBareBeamModel:
            % The following function takes the beam model used by Aeroflex and
            % outputs the bare minimum needed for other analyses to run. This
            % should ideally be in objects. The script currently converts the
            % workspace to a structure which mimics an object, therefore,
            % should be very easy to transfer to objects.
            
            % Remove the midnode generated by Aeroflex
            MidNodeIdx = nonzeros(Input.Beam.Conn(:,2));
            
            % This is where the new beam model will sit
            Output = [];
            
            % Grids IDs and Coordinates
            
            Output.Node.ID     = Input.Node.ID;
            Output.Node.Coord  = Input.Node.Coord;
            
            
            Output.Node.ID(MidNodeIdx)      = [];
            Output.Node.Coord(MidNodeIdx,:) = [];
            
            % Beam Connections and Orientations
            
            Output.Beam.ID         = Input.Beam.ID;
            Output.Beam.PID        = Input.PBeam.ID(Input.Beam.PID);
            Output.Beam.Conn(:,1)  = Input.Node.ID(Input.Beam.Conn(:,1));
            Output.Beam.Conn(:,2)  = Input.Node.ID(Input.Beam.Conn(:,3));
            Output.Beam.Orient     = Input.Beam.Orient;
            Output.Beam.Offset     = Input.Beam.Offset;
            
            
            % Beam Properties
            
            Output.PBeam.ID         = Input.PBeam.ID;
            Output.PBeam.Mat        = Input.PBeam.Mat;
            Output.PBeam.A          = Input.PBeam.A;
            Output.PBeam.I          = Input.PBeam.I;
            Output.PBeam.J          = Input.PBeam.J;
            Output.PBeam.RhoNS      = Input.PBeam.RhoNS;
            Output.PBeam.Kshear     = Input.PBeam.Kshear;
            Output.PBeam.X_L        = Input.PBeam.X_L;
            
            
            % Material Properties
            
            Output.Mat.ID   = Input.Mat.ID;
            Output.Mat.E    = Input.Mat.E;
            Output.Mat.G    = Input.Mat.G;
            Output.Mat.nu   = Input.Mat.nu;
            Output.Mat.Rho  = Input.Mat.Rho;
            
            
            % Lumped Masses
            
            Output.Conm2.ID     = Input.Conm2.ID;
            Output.Conm2.CID    = Input.Conm2.CID;
            Output.Conm2.Node   = Input.Conm2.Node;
            Output.Conm2.M      = Input.Conm2.M;
            Output.Conm2.Offset = Input.Conm2.Offset;
            
            
            % Spatial Constraints
            
            Output.SPC.ID       = Input.SPC.ID;
            Output.SPC.DOF      = Input.SPC.DOF;
            Output.SPC.Nodes    = Input.SPC.Nodes;
            
            
            % Aerodynamic Properties
            Output.Aero.ID      = Input.Aero.ID;
            Output.Aero.CP      = Input.Aero.CP;
            Output.Aero.geo.ny      = Input.Aero.geo.ny;
            Output.Aero.geo.nx      = Input.Aero.geo.nx;
            Output.Aero.geo.startx  = Input.Aero.geo.startx;
            Output.Aero.geo.starty  = Input.Aero.geo.starty;
            Output.Aero.geo.startz  = Input.Aero.geo.startz;
            Output.Aero.geo.c       = Input.Aero.geo.c;
            Output.Aero.geo.b       = Input.Aero.geo.b;
            Output.Aero.geo.T       = Input.Aero.geo.T;
            Output.Aero.geo.SW      = Input.Aero.geo.SW;
            Output.Aero.geo.TW      = Input.Aero.geo.TW;
            Output.Aero.geo.dihed   = Input.Aero.geo.dihed;
            Output.Aero.lattice     = Input.Aero.lattice_vlm;
            Output.Aero.body        = Input.Aero.body;
            
            % Part IDs
            if ~isa(obj, 'mvc.mixin.Pathadjustable') || obj.Model.Root.CustomCodeFolderIndex == 1
                
                %Copy it across
                Output.PartId = Input.PartId;
                
            else
                
                %PartID not present
                
            end
            
        end
    end
    
    methods (Static)
        
        function GUI_Input = InitialiseGUI
            
            % This is where objects might be quite nice
            
            % Initialise some GUI options here (*.Sizing.*)
            GUI_Input.Parameters                = [];
            
            GUI_Input.Parameters.Sizing         = 1;
            GUI_Input.Parameters.Analysis       = 1;
            
            GUI_Input.Parameters.box_type       = 3;
            GUI_Input.Parameters.Ref_model      = 1;
            GUI_Input.Parameters.FakeVTP        = 0;
            GUI_Input.Parameters.WorkshopOutput = 0;
            GUI_Input.Parameters.WorkshopInput  = 0;
            GUI_Input.Parameters.NastranSolver  = 1;
            GUI_Input.Parameters.ARStudy        = 0;
            GUI_Input.Parameters.ARRange        = [];
            GUI_Input.Parameters.BeamElement    = 0;
            GUI_Input.Parameters.WorkshopFile   = [];
            GUI_Input.Parameters.ParallelPro    = 0;
            GUI_Input.Parameters.ResultsPath    = [];
            GUI_Input.Parameters.LoadCases      = [];
            GUI_Input.Parameters.ElementMass    = 1;
            GUI_Input.Parameters.AllowableSensitivity = 0;
            GUI_Input.Parameters.SplineType     = 'Surface';
            GUI_Input.Parameters.BeamType       = 1;
            GUI_Input.Parameters.AddSecondMass  = 1;
            GUI_Input.Parameters.AddFuel        = 1;
            GUI_Input.Parameters.Struc_optimisation = 1;
            GUI_Input.Parameters.Aero_optimisation  = 1;
            GUI_Input.Parameters.Jig_optimisation   = 0;
            
            % Initialise the beam_model structure, required for simulations
            GUI_Input.beam_model                = [];
            
            % Initialise the Aircraft_param structure, required for sizing
            GUI_Input.Aircraft_param            = [];
            
            % Initialise the Fame structure, output from FAME file
            GUI_Input.Fame                      = [];
            
            % Initialise the LoadCases structure, require for sizing
            GUI_Input.LoadCases                 = [];
            
            % Initialise the type of distribution (*.Sizing.*)
            GUI_Input.Parameters.AeroDistribution = 'Elliptical';
            
            % Type of Discipline (*.Analysis.*)
            GUI_Input.Parameters.Structural     = 0;
            GUI_Input.Parameters.Aeroelastic    = 0;
            GUI_Input.Parameters.Aerodynamic    = 0;
            
            % Activate Modal Analysis (*.Analysis.*)
            GUI_Input.Parameters.Modal          = 0;
            
            % Static or Dynamic Analysis (*.Analysis.*)
            GUI_Input.Parameters.Static         = 0;
            GUI_Input.Parameters.Dynamic        = 0;
            
            % Linear or Nonlinear Analysis (*.Analysis.*)
            GUI_Input.Parameters.Linear         = 1;
            GUI_Input.Parameters.Nonlinear      = 0;
            
            % Solver (*.Analysis.*)
            GUI_Input.Parameters.Neo            = 1;
            GUI_Input.Parameters.Nastran        = 0;
            
            % Type of Aeroelastic Analysis (*.Analysis.*)
            GUI_Input.Parameters.Flutter        = 0;
            GUI_Input.Parameters.Trim           = 0;
            GUI_Input.Parameters.FixedAngle     = 0;
            GUI_Input.Parameters.trimI          = 0;
            
            % Rigid or Flexible Analysis (*.Analysis.*)
            GUI_Input.Parameters.Rigid          = 0;
            
            % Aerodynamic Solver
            GUI_Input.Parameters.AeroSolver     = 'VLM-HS';
            
            % Structural Damping (*.Dynamic.*)
            GUI_Input.Parameters.StructDamp     = 0;
            
            % Time Integration Parameter (*.Dynamic.*)
            GUI_Input.Parameters.NewmarkParam   = 0;
            GUI_Input.Parameters.Timestep       = 0;
            GUI_Input.Parameters.NSteps         = 0;
            GUI_Input.Parameters.ODESolver      = 'Gen_alpha';
            GUI_Input.Parameters.TightAero      = 0;
            
            % Parallelisation Properties (*.UVLM.*)
            GUI_Input.Parameters.GPUAcc         = 0;
            GUI_Input.Parameters.ParallelPro    = 0;
            
            % UVLM Properties (*.UVLM.*)
            GUI_Input.Parameters.WakeRollup     = 0;
            GUI_Input.Parameters.WakeCutOff     = 0;
            
            % Gust Properties (*.UVLM.*)
            GUI_Input.Parameters.GustFamily     = 0;
            GUI_Input.Parameters.NDir           = 1;
            GUI_Input.Parameters.MinGustDir     = 0;
            GUI_Input.Parameters.MaxGustDir     = 0;
            GUI_Input.Parameters.NGusts         = 0;
            GUI_Input.Parameters.MinGust        = 0;
            GUI_Input.Parameters.MaxGust        = 0;
            GUI_Input.Parameters.Gust           = 0;
            
            % Some randeom parameters to organise
            GUI_Input.Parameters.Restart        = 0;
            GUI_Input.Parameters.IDrag          = 0;
            GUI_Input.Parameters.AeroLinear     = 0;
            GUI_Input.Parameters.RotGrav        = 0;
            
            GUI_Input.Parameters.WingAngle      = 0;
            
            GUI_Input.Parameters.ResultsPath    = [pwd filesep 'FrameworkResults'];
            
            GUI_Input.Parameters.ProfileCode    = 0;
            
            GUI_Input.Results                   = [];
            
        end
        
        
        
    end
    
end
