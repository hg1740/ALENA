%%SOLVE_NONLINEAR_TRIM_CUSTOMISED
% 
% 
%   Author: Dario Calderon 

function [Static,beam_model] = solve_nonlinear_trim_customised(Static,beam_model,TRIM_INDEX, hp, logfcn)

NSTEP = 20;

PLOT_CONV = 1;
PLOT_ITER = 0;

if PLOT_CONV == 1
    if nargin < 4 || isempty(hp)
                hp = figure;
                ha(1) = subplot(1,2,1,'Parent',hp);
                ha(2) = subplot(1,2,2,'Parent',hp);
        %PLOT_CONV = 0;
        logfcn = @disp;
    elseif isa(hp,'uiextras.HBox')
        delete(hp.Children);
        ha(1) = axes('Parent',hp);
        ha(2) = axes('Parent',hp);
    else
        error('Subplot issues - FIX ME!!!!!!!!!!!!!');
    end
end

WngNodes = beam_model.WingNodes;
GRAV = beam_model.Param.G;
ndof = beam_model.Info.ndof;
nc = 1;

ngrid = beam_model.Info.ngrid;


Fa = zeros(ndof, 2);
Fi = zeros(ndof, 2);
M = ms_matrix(beam_model.Info, beam_model.Node.DOF, beam_model.Node.R, beam_model.ConM, beam_model.Bar, beam_model.Beam);
%F = gf_lin_nodal(beam_model.Param.LOAD, beam_model.Info, beam_model.F, beam_model.M, beam_model.Node.DOF);

%     beam_model.Aero.geo.ref_point = beam_model.WB.CG;
%     beam_model.Aero.geo.CG = beam_model.WB.CG;

%   Total external forces (Applied loads + Inertial + Follower)
SOL = zeros(ndof, 1);

ACC = zeros(6,1);
ACC(1) = 0;
ACC(3) = - (GRAV);
Fi1 = gf_iner_nodal(ndof, beam_model.Node.DOF, M, ACC);

% AERO COUPLED ITERATIVE SOLVER

%fprintf(['\nPerforming nonlinear trim ... for load case ' num2str(i) '. #' num2str(trimcount) ' out of ' num2str(length(TRIM_INDEX))]);
for N = 2 : NSTEP+1
    
    if N>2
        beam_model.Aero.lattice = update_vlm_mesh(beam_model.Node, node_pos, ...
            beam_model.Res.NRd, beam_model.Aero,beam_model);
        if nc
            beam_model.Aero.lattice = update_hinge_lines(beam_model.Aero.geo, ...
                beam_model.Aero.lattice);
        end
        %restart_rigid_trim(beam_model.Res.Aero.x_hist(:,end), 1);
        beam_model.Aero.AIC = [];
        %beam_model = rigid_vlm_trim(beam_model);
        beam_model.Res = rigid_trim(beam_model,beam_model.Res,beam_model.Aero);
        
        
    end
%     logfcn(sprintf('Coupled Iteration %g, Alpha: %f deg', N, rad2deg(beam_model.Res.state.alpha)));
%     %fprintf('... Alpha: %f deg.',   rad2deg(beam_model.Res.state.alpha));
%     if beam_model.Aero.geo.nc
%         logfcn('Control Surfaces:');
%         [~,J] = find(beam_model.Res.CS.Fixed == 0);
%         for n_nc = J
%             logfcn(sprintf(' %s: %g deg.', beam_model.Aero.lattice.Control.Name{1,n_nc},180*beam_model.Res.CS.Value(n_nc)/pi));
%         end
%     end
    beam_model.Res.SOL = 'Static nonlinear aeroelastic';
    
    %%
    
    Fi(:,2) = Fi1;
    
    Fa(:,2) = gf_trim_aero_nodal(beam_model.Info, beam_model.Node.DOF, ...
        beam_model.Node, beam_model.Aero, beam_model.Res.Aero,1);
    
    AeroWeight = 0.8;
    
    Fa(:,1) = AeroWeight* Fa(:,2)+(1-AeroWeight)*Fa(:,1);
    Fi(:,1) = AeroWeight* Fi(:,2)+(1-AeroWeight)*Fi(:,1);
    
    Fext = Fa(:,1) + Fi(:,1);
    
    %% PERFORM NONLINEAR STATIC ANALYSIS
    
    %[node_pos,node_rot,~,~,beam_model,~] = nonlinear_static_solver(Fext,beam_model,15,0);
    [node_pos,beam_model,~] = nonlinear_static_solver(Fext,beam_model,15,0);

    DISP_RMS(N)   = sqrt(sum(sum(beam_model.Res.NDispl(:,1:3) .* beam_model.Res.NDispl(:,1:3),2),1));
    DS_RESOVERALL = abs(DISP_RMS(N) - DISP_RMS(N-1));
    RESIDUAL_AERO(N)=DS_RESOVERALL;
    TipDisplacement(N)=beam_model.Res.NDispl(WngNodes(end),3);
    
    if PLOT_CONV ==1
        if N >= 1
            
            semilogy(ha(1),(2:N), RESIDUAL_AERO(2:N), '-ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'LineWidth', 1);
            grid(ha(1),'on'); xlabel(ha(1),'Coupled iteration'); ylabel(ha(1),'Displacements variation norm');
            
            plot(ha(2),(1:N), TipDisplacement(1:N), '-ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'LineWidth', 1);
            grid(ha(2),'on'); xlabel(ha(2),'Coupled iteration'); ylabel(ha(2),'Tip Displacement (m)');
        end
    end
    
    if PLOT_ITER == 1
        figure(400);hold on;
        plot(beam_model.Node.Coord(WngNodes,2),beam_model.Node.Coord(WngNodes,3),'b--');hold on
        plot(node_pos(WngNodes,2),node_pos(WngNodes,3),'r-');hold on
    end
    
    if (DS_RESOVERALL < 5e-4)
        %logfcn(sprintf(['\n\t!! Nonlinear trim solution found in ' num2str(N) ' steps!!\n']));
        break;
    elseif (DS_RESOVERALL > 5e-4) && N == NSTEP+1
        %logfcn(sprintf('\n\t!! WARNING - Nonlinear trim solution NOT FOUND. !!\n'));
    end
    
end
% 
% results = TrefftzLineDrag(beam_model.Res.Aero.results,beam_model.Aero.state,beam_model.Aero.lattice,beam_model.Aero.ref);
% beam_model.Res.Aero.results.CD = results.Trefftz_drag_Coeff;
% beam_model.Res.Aero.results.treffts = results.treffts;

if ~isempty(beam_model.Res.Bar.CForces)
    IntLoads = storeinternalloads(beam_model.Res.Bar.CForces);
    BarForces.Fx(:,1) = IntLoads.Fx;BarForces.Fy(:,1) = IntLoads.Fy;BarForces.Fz(:,1) = IntLoads.Fz;
    BarForces.Mx(:,1) = IntLoads.Mx;BarForces.My(:,1) = IntLoads.My;BarForces.Mz(:,1) = IntLoads.Mz;
else
    IntLoads = storeinternalloads(beam_model.Res.Beam.CForces);
    BarForces.Fx(:,1) = IntLoads.Fx;BarForces.Fy(:,1) = IntLoads.Fy;BarForces.Fz(:,1) = IntLoads.Fz;
    BarForces.Mx(:,1) = IntLoads.Mx;BarForces.My(:,1) = IntLoads.My;BarForces.Mz(:,1) = IntLoads.Mz;
end

beam_model.Lattice.flexible{TRIM_INDEX} = beam_model.Aero.lattice;
%beam_model.Loads.flexible{TRIM_INDEX}   = storeloads(beam_model);
beam_model.Results.flexible{TRIM_INDEX} = beam_model.Res;

NR = beam_model.Res.NRd;
ngrid = beam_model.Info.ngrid;
ex_node = zeros(ngrid,3);
ey_node = zeros(ngrid,3);
ez_node = zeros(ngrid,3);

for k = 1:ngrid
    ex_node(k,1)=NR(1,1,k);
    ey_node(k,1)=NR(1,2,k);
    ez_node(k,1)=NR(1,3,k);
    ex_node(k,2)=NR(2,1,k);
    ey_node(k,2)=NR(2,2,k);
    ez_node(k,2)=NR(2,3,k);
    ex_node(k,3)=NR(3,1,k);
    ey_node(k,3)=NR(3,2,k);
    ez_node(k,3)=NR(3,3,k);
end

Displacements = zeros(ngrid,6);
Displacements(:,1:3) = beam_model.Res.NDispl(:,1:3);
Displacements(:,4) = -atan(ez_node(:,2)./ez_node(:,3));
Displacements(:,5) = atan(ez_node(:,1)./ez_node(:,3));
Displacements(:,6) = atan(ex_node(:,2)./ex_node(:,1));

beam_model.Res.Aero.results.gamma = beam_model.Res.Aero.gamma;
beam_model.Res.Aero.results.CN   = beam_model.Res.Aero.CN ;
Static.Displacements{TRIM_INDEX}        = Displacements;
Static.Positions{TRIM_INDEX}            = node_pos;
Static.NonDimAeroForces{TRIM_INDEX}     = [beam_model.Res.Aero.spancd,...
    beam_model.Res.Aero.spancs,beam_model.Res.Aero.spancl];
Static.p_mid_r{TRIM_INDEX}              = beam_model.Res.Aero.p_mid_r;
Static.BarForces{TRIM_INDEX}            = BarForces;
Static.State{TRIM_INDEX}                = beam_model.Res.state;
Static.Gamma{TRIM_INDEX}                = beam_model.Res.Aero.gamma;
%logfcn(sprintf('\n\n\t... Tip Displacement:  %f m.\n\n',   Displacements(beam_model.WingNodes(end),3)));

end