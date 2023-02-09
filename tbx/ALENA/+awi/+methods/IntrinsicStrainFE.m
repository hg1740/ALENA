classdef IntrinsicStrainFE < awi.methods.Analysis
    %IntrinsicStrainFE Contains the methods for analysing the nonlinear,
    %instrinsic, strain-based finite element methodology.
    
    methods % static analysis
        function StaticResults = static(obj, StaticOptions, load)
            %static Executes a static analysis using the analysis model
            %assigned to this method object and the options defined in
            %'StaticOptions'.            
                        
            %import awi.methods.isfe.* 
            
            StaticResults = [];
            
            if strcmp(StaticOptions.StructuralResponse, 'linear')
                warning(['The intrinsic, strain-based, finite element ', ...
                    'structural solver is only valid for nonlinear '   , ...
                    'structural analysis. Setting parameter '          , ...
                    '''StructuralResponse'' to ''nonlinear''.']);
                StaticOptions.StructuralResponse = 'nonlinear';
            end
            
            %Pass it on
            static@awi.methods.Analysis(obj, StaticOptions);
            
            AllFEM = flatlist(obj.AnalysisModel);
            
            %% Go via Robbie's original code
            BeamModel = convertFEM2BeamModel(AllFEM);          
                                    
            function BeamModel = convertFEM2BeamModel(AllFEM)
                %convertFEM2BeamModel Uses the awi.fe.FEModel object to
                %return a MATLAB structure matching the 'aircraft_data'
                %input structure to the code 'InitialiseRobbie'.                
                                
                BeamModel = i_initialiseBeamModel;   
                
                %awi.fe objects
                BeamNodes = [AllFEM.Beams.Nodes];
                Nodes     = [BeamNodes(1, :), BeamNodes(2, end)];
                BeamProps = [AllFEM.BeamProps];
                BeamMat   = [BeamProps.Material];
                
                nBeam =  numel(AllFEM.Beams);
                
                %Part IDs
                PartId = struct('Part', 'PortWing', 'Type', [], 'data', []);
                PartId = repmat(PartId, [1, 2]);
                %   - Beam grids
                PartId(1).Type = 'GRID';
                PartId(1).data = [Nodes.ID];
                %   - Aero panels
                PartId(2).Type = 'CAERO';
                PartId(2).data = [AllFEM.AeroPanels.ID];
                
                BeamModel.PartId = PartId;
                
                %Nodes
                BeamModel.Node.ID    = [Nodes.ID]';
                BeamModel.Node.Coord = [Nodes.X]';
                
                %Beams
                BeamModel.Beam.ID     = [AllFEM.Beams.ID];
                BeamModel.Beam.PID    = [BeamProps.ID];
                BeamModel.Beam.Conn   = [BeamNodes(1, :).ID ; BeamNodes(2, :).ID]';
                BeamModel.Beam.Orient = repmat([1, 0, 0], [numel(BeamNodes(1, :)), 1]);
                BeamModel.Beam.Offset = zeros(numel(BeamNodes(1, :)), 9);
                
                %Beam properties
                BeamModel.PBeam.ID  = [BeamProps.ID];
                BeamModel.PBeam.MID = [BeamMat.ID];
                %   - Area
                BeamModel.PBeam.A = arrayfun(@(ii) struct('data', BeamProps(ii).A), 1 : nBeam);
                BeamModel.PBeam.I = arrayfun(@(ii) struct('data', [BeamProps(ii).I11 , BeamProps(ii).I22, BeamProps(ii).I12]'), 1 : nBeam);
                BeamModel.PBeam.J = arrayfun(@(ii) struct('data', BeamProps(ii).J), 1 : nBeam);
                
                %Materials
                BeamModel.Mat.ID  = [AllFEM.Materials.ID];
                BeamModel.Mat.E   = [AllFEM.Materials.E];
                BeamModel.Mat.G   = [AllFEM.Materials.G];
                BeamModel.Mat.Rho = [AllFEM.Materials.Rho];
                
                %Aero 
                rStart = [AllFEM.AeroPanels.X1];
                chord  = [AllFEM.AeroPanels.CHORD];
                BeamModel.Aero.ID         = [AllFEM.AeroPanels.ID]';
                BeamModel.Aero.geo.c      = chord(1, :)';
                BeamModel.Aero.geo.startx = rStart(1, :)';
                BeamModel.Aero.geo.starty = rStart(2, :)';
                BeamModel.Aero.geo.startz = rStart(3, :)';
                BeamModel.Aero.geo.TW     = zeros(numel(AllFEM.AeroPanels), 1);
                
                function BeamModel = i_initialiseBeamModel 
                                       
                    Node  = struct('ID', [], 'Coord', []);                   
                    Beam  = struct('ID', [], 'PID', [], 'Conn', [], 'Orient', []);
                    PBeam = struct('ID', [], 'MID', [], 'A', [], 'I', [], 'J', []);
                    Mat   = struct('ID', [], 'E', [], 'G', []);
                    Conm2 = struct('Node', [], 'M', [], 'Offset', []);
                    Aero  = struct('ID', [], 'geo'  , struct( ...
                        'c', [], 'startx', [], 'starty', [], 'startz', [], 'TW', []));
                    
                    BeamModel = struct('Node', Node, 'Beam', Beam, ...
                        'PBeam', PBeam, 'Mat', Mat, 'Conm2', Conm2, 'Aero', Aero);
                    
                end
                
            end
            
            [Matrices ,AeroData ,MassProps] = InitialiseRobbie(BeamModel);
            
            %Tip load in positive z
            switch AllFEM.PointLoads(1).LoadBehaviour
                case 'follower'
                    Matrices.f{1}(end - 3)   = -AllFEM.PointLoads(1).Magnitude;
                case 'non-follower'
                    Matrices.f_a{1}(end - 3) = AllFEM.PointLoads(1).Magnitude;
            end
            
            
            Matrices.CGa  = eye(3);
            
            x0.x_f        = {zeros(size(Matrices.f_mp_pnt{1}))};              
            Sim           = initSystem(x0, Matrices, AeroData, StaticOptions);
            x0            = initSim(x0, Matrices, StaticOptions);
            Sim.Soln      = 1;
            Sim.aero_flag = 0;
            Sim.rb_flag   = 1;            
            AeroData.inflow_switch = 1;

            %Run static solution
            x_stat = getStatic(x0, Matrices, AeroData, Sim, StaticOptions);
            
            [coords, CaB] = strains2coords_all(x_stat.x_f, Matrices);
            
            %% New ALENA implementation             
         
            %Make the model
            AllFEM = flatlist(obj.AnalysisModel);
            ModelData = arrayfun(@(~) awi.methods.isfe.SlenderBody, 1 : numel(AllFEM));
            setup(ModelData, AllFEM, 'static');
            
            %Setup initial system 
            System = awi.methods.isfe.System(ModelData, StaticOptions);
            
            %Run analysis
            x_f = awi.methods.isfe.getStatic(ModelData, System, StaticOptions);

            [coords_new, CaB_new] = strains2coords(ModelData, x_f.x_f);
            [coords_0  , ~] = strains2coords(ModelData, zeros(size(x_f.x_f)));
            dT = coords_new - coords_0;
            
            %Generate the results sets
            BN = AllFEM.BeamNodes;
            DispResults = awi.results.NodalResults;
            DispResults.Nodes       = [BN(1, :), BN(2, end)];
            DispResults.Translation = dT(:, ModelData.MatrixIndexing.NodeInd);
            
            StaticResults = DispResults;
            
            %% error  calculation
            
            figure, hold on
            plot3(squeeze(coords{1}(1, :, :)), squeeze(coords{1}(2, :, :)), squeeze(coords{1}(3, :, :)), 'b-', 'DisplayName', 'FALCON Code');
            plot3(coords_new(1, :), coords_new(2, :), coords_new(3, :), 'r-', 'DisplayName', 'ALENA Code');
            legend
            
            rmse = sqrt(mean((abs(x_stat.x_f{1}) - abs(x_f.x_f)).^2));
            if rmse > 1
                warning('Bad results between FALCON and ALENA code.');
            end
            
            
        end
        
    end
    
end

