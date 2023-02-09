function [beam_model,Aircraft_param,Optim_Variables] = setupbox(Aircraft_param,box_type)
%
cd_current = cd;
filename = [cd filesep 'WingSizing.dat'];
fid = 1;

beam_model = InitRead(filename,fid);

% %% Store structure array
% %warning('OFF','MATLAB:structOnObject')
% beam_model.COORD1 = COORD1;
% beam_model.COORD2 = COORD2;
% beam_model.Coord  = Coord;
% beam_model.Node   = Node;
% beam_model.Bar    = Bar;
% beam_model.Beam   = Beam;
% beam_model.PBar   = PBar;
% beam_model.PBeam  = PBeam;
% beam_model.Aero   = Aero;
% beam_model.RBE2   = RBE2;
% beam_model.RBE0   = RBE0;
% beam_model.Param  = Param;
% beam_model.Conm1  = Conm1;
% beam_model.Conm2  = Conm2_;
% beam_model.Mat    = Mat_;
% beam_model.Info   = Info_;
% beam_model.SPC    = SPC;
% beam_model.F      = FORCE;
% beam_model.Gust   = GUST;
% beam_model.Aeros  = AEROS;
% beam_model.Dextload = DEXTLOAD;
% beam_model.M      = MOMENT;
% beam_model.F_FLW  = F_FLW;
% beam_model.Celas  = CELAS;
% beam_model.CBush  = CBUSH;
% beam_model.PBush  = PBUSH;
% beam_model.Aero.body  = BAero;
% beam_model.RJoint = RJOINT;
% beam_model.WB     = [];
% beam_model.RBar   = RBar;
% beam_model.BAero  = BAero;

%warning('ON','MATLAB:structOnObject')

%%%%% 2 - No Stringers, 3 - With Stringers, 4 - Smeared Stringers %%%%%
ACParts = fieldnames(Aircraft_param);
wingIdx = ~cellfun(@isempty,regexpi(ACParts,'wing'));
WingPartName = ACParts(wingIdx);

NWingBeams   = Aircraft_param.(WingPartName{1}).Planform.NBeams;
beam_model.Engine.Exist = 0;

AircraftParts = fieldnames(Aircraft_param);

PartCount = 0;
for i = 1:length(AircraftParts)
    
    if isfield(Aircraft_param.(AircraftParts{i}),'LiftingSurface')
        PartCount = PartCount + 1;
        if Aircraft_param.(AircraftParts{i}).LiftingSurface == 1
            
            % Define the starting index for each of the parts
            Aircraft_param.(AircraftParts{i}).Planform.NodeIdx    = 100000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.LENodeIdx  = 201000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.TENodeIdx  = 202000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.RBE0Idx    = 100000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.CBEAMIdx   = 100000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.PBEAMIdx   = 100000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.CONM2Idx   = 800000 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.SET1Idx    = 500001 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.SPLINE1Idx = 500001 + PartCount*10000;
            Aircraft_param.(AircraftParts{i}).Planform.CAEROIdx   = 710000 + PartCount*10000;
            
            [beam_model.Aero.(AircraftParts{i}),Optim_Variables.(AircraftParts{i}),Aircraft_param.(AircraftParts{i}).Planform] = ...
                setupWingV3(box_type,Aircraft_param.(AircraftParts{i}).Planform,Aircraft_param.(AircraftParts{i}).Box);
            
            beam_model.PartIDs.(AircraftParts{i}).CAEROIDs = beam_model.Aero.(AircraftParts{i}).ID;
            beam_model.PartIDs.(AircraftParts{i}).Mat      = PartCount*100;
            
            if isfield(Aircraft_param.(AircraftParts{i}).Planform,'CT_y_eta')
                y_ct = Aircraft_param.(AircraftParts{i}).Planform.CT_y_eta*Aircraft_param.(AircraftParts{i}).Planform.b_ref;
                
                [~,ct_idx] = min(abs(y_ct - Optim_Variables.(AircraftParts{i}).nodey));
                Aircraft_param.(AircraftParts{i}).Planform.CTNode = Aircraft_param.(AircraftParts{i}).Planform.NodeIdx + ct_idx;
            end
            
        else
            
            Aircraft_param.(AircraftParts{i}).Planform.NodeIdx     = 100000 + PartCount*10000;
             Aircraft_param.(AircraftParts{i}).Planform.CONM2Idx   = 800000 + PartCount*10000;
            % Check to see if it has a parent object
            if ~isempty(Aircraft_param.(AircraftParts{i}).Parent)
                
                Parent = Aircraft_param.(AircraftParts{i}).Parent;

                % Check if the parent object is a lifting surface or not
                
                if Aircraft_param.(Parent).LiftingSurface == 1
                    
                    % Update the offset so that it is correctly placed
                    Node = [Optim_Variables.(Parent).nodex,Optim_Variables.(Parent).nodey,Optim_Variables.(Parent).nodez];
                    
                    deltaNode = Node(2:end,:)-Node(1:end-1,:);
                    
                    % Find the length along the part to identify the dominating
                    % axis. This is important for distinguishing VTPs and wings
                    % for example.
                    PartLength = sum(deltaNode);
                    
                    if Aircraft_param.(Parent).Planform.dihedral(1).value_beg == 90
                        DomAxis = 3;
                    else
                        DomAxis = 2;
                    end
                    
                    
                    % Determine where along the parent axis the body will
                    % lie
                    if ~isempty(Aircraft_param.(AircraftParts{i}).Planform.y_eta_parent)
                        
                        LocationWithParent = Aircraft_param.(AircraftParts{i}).Planform.y_eta_parent*(Node(end,DomAxis)-Node(1,DomAxis));
                        
                    elseif ~isempty(Aircraft_param.(AircraftParts{i}).Planform.y_parent)
                        
                        LocationWithParent = Aircraft_param.(AircraftParts{i}).Planform.y_parent-Node(1,DomAxis);
                        
                    end
                    
                    % Find the difference between the location of the body
                    % and the nodes assigned to the parent body
                    Loc_diff = Node(:,DomAxis) - LocationWithParent;
                    
                    % Determine which node to attach the body to
                    [~,Parent_index] = min(abs(Loc_diff));
                    
                    Node_offset            = zeros(1,3);
                    Node_offset(1,DomAxis) = -Loc_diff(Parent_index);
                    
                    % No offset exists then use the CG position to work out
                    % the node offset and the reference point.
                    if isempty(Aircraft_param.(AircraftParts{i}).Planform.Offset)
                        
                        Aircraft_param.(AircraftParts{i}).Planform.Offset =  Aircraft_param.(AircraftParts{i}).Properties.CG - ...
                            Node(Parent_index,:) + Node_offset - [Aircraft_param.(AircraftParts{i}).Planform.Length/2,0,0];
                        
                         Aircraft_param.(AircraftParts{i}).Planform.ref_point = Aircraft_param.(AircraftParts{i}).Properties.CG ...
                            - [Aircraft_param.(AircraftParts{i}).Planform.Length/2,0,0];
                    
                    else % if an offset exists then use that to determine the CG, offset and ref_point.
                                              
                         Aircraft_param.(AircraftParts{i}).Planform.Offset = Node_offset + ...
                             Aircraft_param.(AircraftParts{i}).Planform.Offset; 
                         
                        Aircraft_param.(AircraftParts{i}).Planform.ref_point = Node(Parent_index,:) + ...
                            Aircraft_param.(AircraftParts{i}).Planform.Offset +...
                            - [Aircraft_param.(AircraftParts{i}).Planform.Length/2,0,0];
                        
                        Aircraft_param.(AircraftParts{i}).Properties.CG = Node(Parent_index,:) + ...
                            Aircraft_param.(AircraftParts{i}).Planform.Offset;
                        
                    end
                    
                    Aircraft_param.(AircraftParts{i}).Planform.ParentIdx = Parent_index;
                else
                    
                    Aircraft_param.(AircraftParts{i}).Planform.ParentIdx = 1;
                    
                    Aircraft_param.(AircraftParts{i}).Planform.ref_point = Aircraft_param.(AircraftParts{i}).Properties.CG ...
                        - [Aircraft_param.(AircraftParts{i}).Planform.Length/2,0,0];
                    %error('Setting the engine as a part of the fuselage needs to be fixed');
                end
            else
                Aircraft_param.(AircraftParts{i}).Planform.ref_point = Aircraft_param.(AircraftParts{i}).Properties.CG ...
                    - [Aircraft_param.(AircraftParts{i}).Planform.Length/2,0,0];
            end
            
            [beam_model,beam_model.Aero.body] = setupFuselage(beam_model,beam_model.Aero.body,Aircraft_param.(AircraftParts{i}).Planform);
            
            beam_model.PartIDs.(AircraftParts{i}).BodyIdx = beam_model.Info.nbaero;
            
        end
    end
end

% List lifting surfaces that do not lie on the symmetry plane

SizLift = {};
for i = 1:numel(ACParts)
    Part = ACParts{i};
    if isfield(Aircraft_param.(Part),'LiftingSurface')
        if Aircraft_param.(Part).LiftingSurface == 1 && Aircraft_param.(Part).Planform.SymPlane == 0
            SizLift = [SizLift;ACParts{i}];
        end
    end
    
end

% Combine the Aerodynamic mesh onto one
beam_model = Combine_Aero(beam_model,SizLift);

% Define the wing nodes
beam_model.WingNodes = (1: NWingBeams+1);
beam_model.WingBar = (1: 2*NWingBeams);

% Define mean bar
beam_model.WingMeanBar  = (1: NWingBeams);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% Generate the VORTEX LATTICE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[beam_model.Aero.lattice, beam_model.Aero.ref] = ...
    vlm_setup(1, beam_model.Aero.geo, beam_model.Aero.state, beam_model.Aero.ref);

beam_model.Aero.lattice.Control.Name = beam_model.Aero.Control.Name;

beam_model.Aero = rmfield(beam_model.Aero,'Control');

beam_model.Aero.lattice_vlm = beam_model.Aero.lattice;

beam_model.Info.amesh_av_vlm = 1;
beam_model.Info.ncaero = numel(beam_model.Aero.geo.ny);

beam_model.twist_corr = 1;
beam_model.camber_corr = 1;

end