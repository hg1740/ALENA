%% solve_static
% The following code is a static solver where the applied loads are
% derived from force, moment and thrust cards. A given gravitational
% constant can be given.
%
% Outputs: Displacements, Internal loads and positions.
% Restart ouputs can be requested.

function [Static,beam_model] = AeroFlexStatic(beam_model,Linear)

beam_model.Param.SUBCASE{1}.LOAD    = 1;
beam_model.Param.SUBCASE{1}.SPC     = 1001;

beam_model.Grav.ID                  = 1;
beam_model.Grav.Scale               = 29.0562;
beam_model.Grav.Orient              = [0,0,-1];

beam_model.F.ID = 1;
beam_model.F.Node = 56;
beam_model.F.Mag = 0;
beam_model.F.Orient = [0,0,1];
beam_model.Info.nf = 1;

% Run through the subcases
Static = [];

if isempty(beam_model.Param.SUBCASE{1}.ID)
    disp('No Subcases defined.');
    return
end

for i = 1:length(beam_model.Param.SUBCASE)
    
    beam_model = getSubcase(beam_model,i);
    
    % Define some useful counters
    ndof  = beam_model.Info.ndof;
    ngrid = beam_model.Info.ngrid;
    
    % Define the mass matrix
    M = ms_matrix(beam_model.Info, beam_model.Node.DOF, beam_model.Node.R, beam_model.ConM, beam_model.Bar, beam_model.Beam);
    
    % Define the stiffness matrix here for a Linear Analysis
    if Linear
        K = st_lin_matrix(beam_model.Info, beam_model.Node.DOF, beam_model.Node.R, ...
            beam_model.Node.Coord, beam_model.Bar, beam_model.Beam, beam_model.Celas,beam_model.CBush);
        if ~isempty(beam_model.RBE2)
            K = RBE2Assembly(beam_model.RBE2,K);
        end
        
        % Account for suport cards
        [Kll,ldof] = get_free_stiff(K, beam_model.Node, beam_model.Param.SUPORT, beam_model.Param.EPS);
        
    end
    
    % Define the acceleration vector
    ACC = zeros(6,1);
    ACC(1:3) = beam_model.Param.GRAV';    
    
    % Define inertial load
    Fi = gf_iner_nodal(ndof, beam_model.Node.DOF, M, ACC);
    
    % Define applied loads
    if ~isempty(beam_model.Param.LOAD)
        Fapplied = gf_lin_nodal(beam_model.Param.LOAD, beam_model.Info, beam_model.F, beam_model.M, beam_model.Node.DOF);
    else
        Fapplied = zeros(ndof, 1);
    end
    
%     % Process the pre-strain in the wing: TODO FIX THIS for the stiffness
%     F_pre = ...
%         get_precurvature_load(beam_model.Info.nbeam, beam_model.Beam, beam_model.Node, ...
%         zeros(ngrid,6), beam_model.Info.ngrid,beam_model.Beam.PreStrain);
%     
%     F_curv = get_internal_forces(ndof, beam_model.Node.DOF, F_pre);
    
    % Combine loads
    Fext = Fapplied + Fi;% - F_curv;
    
    % NONLINEAR ANALYSIS STATIC
    if Linear
        fprintf('\nRunning Linear static analysis ...')
        [NODEPOS,beam_model,BarForces] = linear_static_solver(Fext,beam_model,Kll,ldof,0);
        
    else
        fprintf('\nRunning Nonlinear static analysis ...')
        [NODEPOS,beam_model,BarForces] = nonlinear_static_solver(Fext,beam_model,40,1);
    end
    
    % Restructure and calculate rotational displacements
    ex_node = zeros(ngrid,3);
    ey_node = zeros(ngrid,3);
    ez_node = zeros(ngrid,3);
    
    for k = 1:ngrid
        ex_node(k,1)=beam_model.Res.NRd(1,1,k);
        ey_node(k,1)=beam_model.Res.NRd(1,2,k);
        ez_node(k,1)=beam_model.Res.NRd(1,3,k);
        ex_node(k,2)=beam_model.Res.NRd(2,1,k);
        ey_node(k,2)=beam_model.Res.NRd(2,2,k);
        ez_node(k,2)=beam_model.Res.NRd(2,3,k);
        ex_node(k,3)=beam_model.Res.NRd(3,1,k);
        ey_node(k,3)=beam_model.Res.NRd(3,2,k);
        ez_node(k,3)=beam_model.Res.NRd(3,3,k);
    end
    
    Displacements = zeros(ngrid,6);
    Displacements(:,1:6)    = beam_model.Res.NDispl(:,1:6);
    if Linear == 0
        Displacements(:,4)      = atan(ey_node(:,3)./ez_node(:,3));
        Displacements(:,5)      = atan(-ex_node(:,3)./sqrt(ey_node(:,3).^2 + ez_node(:,3).^2));
        Displacements(:,6)      = atan(ex_node(:,2)./ex_node(:,1));
    end
    
    Static.Displacements{i}     = Displacements;
    Static.BarForces{i}         = BarForces;
    Static.Positions{i}         = NODEPOS;
    Static.NR{i}                = beam_model.Res.NRd;
    Static.BARR{i}              = beam_model.Res.Bar.R;
    Static.BEAMR{i}             = beam_model.Res.Beam.R;
end


end
