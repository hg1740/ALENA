function [NODEPOS,beam_model,BarForces] = nonlinear_static_solver(F,beam_model,NSTEP,ENABLE_DISP)

%% SETUP MATRICES FOR SOLVER
%beam_model.WingNodes = SortWingNodes(beam_model);

MaxIter = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%% COUNTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ngrid = beam_model.Info.ngrid;
nbar  = beam_model.Info.nbar;
nbeam = beam_model.Info.nbeam;
ndof  = beam_model.Info.ndof;
ndof2 = beam_model.Info.ndof2;

%%%%%%%%%%%%%%%%%%%%%% INITIAL STRAINS & CURVATURES %%%%%%%%%%%%%%%%%%%%%%%
BARPO  = set_initial_PO(nbar, beam_model.Bar, beam_model.Node);
BEAMPO = set_initial_PO(nbeam, beam_model.Beam, beam_model.Node);
BARKR  = zeros(2, 3, nbar);
BEAMKR = zeros(2, 3, nbeam);
%BARKR(:,:,14) = 2.7*[7.053602906274305e-05,7.483471347882127e-06,0.434240281433032;8.042102483402539e-05,5.429081866461154e-06,0.339704181207701];
%BARKR(:,:,36) = 2.7*[-7.053602906274305e-05,-7.483471347882127e-06,0.434240281433032;-8.042102483402539e-05,-5.429081866461154e-06,0.339704181207701];
% BEAMKR = 1/102.5*repmat([0,1,0;0,1,0],1,1,nbeam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERNAL FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beam_model.Res.Bar.CForces  = zeros(2, 6, nbar);
beam_model.Res.Beam.CForces = zeros(2, 6, nbeam);

%%%%%%%%%%%%%%%%%%%%% INTERNAL STRAINS AND CURVATURES %%%%%%%%%%%%%%%%%%%%%%
beam_model.Res.Bar.CStrains  = zeros(2, 6, nbar);
beam_model.Res.Beam.CStrains = zeros(2, 6, nbeam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLACEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beam_model.Res.NDispl        = zeros(ngrid, 6);
beam_model.Res.NRot          = zeros(ngrid, 3);

%%%%%%%%%%%%%%%%%%%%%%%%% NODE ROTATION MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
beam_model.Res.NRd    = beam_model.Node.R;

%%%%%%%%%%%%%%%%%%%%%%%%%% BAR ROTATION MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%
beam_model.Res.Bar.R  = beam_model.Bar.R;
beam_model.Res.Beam.R = beam_model.Beam.R;

NODEPOS = beam_model.Node.Coord;
NR      = beam_model.Res.NRd;
BARR    = beam_model.Res.Bar.R;
BARF    = allocate_barf(nbar);
BEAMR   = beam_model.Res.Beam.R;
BEAMF   = allocate_barf(nbeam);
BUSHR   = beam_model.CBush.R;
JOINTR  = beam_model.RJoint.R;
NR0     = NR;

%%%%%%%%%%%%%%%%%%%%%%%%%% Joints (RJOINTS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TotConst = 0;
lambda_joint = 0;
delta_lambda_joint = 0;
B_joint = 0;
if ~isempty(beam_model.RJoint.ID)
    % Holonomic Contraints multiplier
    lambda_joint = [];
    delta_lambda_joint = [];
    
    for i = 1:length(beam_model.RJoint.ID)
        sizenr = length(beam_model.RJoint.DOFC(i).list);
        lambda_joint = [lambda_joint;zeros(sizenr,1)];
        delta_lambda_joint = [delta_lambda_joint;zeros(sizenr,1)];
        TotConst = TotConst + sizenr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Rigid Elements (RBARS) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_rbar = 0;
delta_lambda_rbar = 0;
B_rbar = 0;
if ~isempty(beam_model.RBar.ID)
    % Holonomic Contraints multiplier
    lambda_rbar = [];
    delta_lambda_rbar = [];
    for i = 1:length(beam_model.RBar.ID)
        sizenr = length(beam_model.RBar.DOFC(i).list);
        lambda_rbar = [lambda_rbar;zeros(sizenr,1)];
        delta_lambda_rbar = [delta_lambda_rbar;zeros(sizenr,1)];
        TotConst = TotConst + sizenr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DSCALE   = 1/NSTEP;
STEPS    = linspace(DSCALE,1,1/DSCALE);
LOAD_SCALE = DSCALE;
NEWSTEP(1) = STEPS(1);
LOAD_COUNT = 0;
RES_TOL  = 10^(-3);
FRES_TOL = 10^(-4);
WngNodes = beam_model.WingNodes;
TipDisp  = 0;

AdaptiveStepping = 1;

SimpleStiff = 0;

Aerodef = zeros(ngrid, 3);
subcount = 5;

while LOAD_SCALE < 1
    %***********************************************************************************************************************
    % erase rotations
    beam_model.Res.NDispl(:, 4:6) = zeros(ngrid, 3);
    DR   = zeros(3,3,beam_model.Info.ngrid);
    gdef = zeros(beam_model.Info.ngrid, 6);
    
    LOAD_COUNT = LOAD_COUNT + 1;
    
    % ADAPTIVE SCALING
    if AdaptiveStepping == 1
        if LOAD_SCALE < 1
            
            % Increase load step
            if subcount < 4
                NEWSTEPS = ceil((1-LOAD_SCALE)/(2*DSCALE));
                DSCALE = DSCALE * 2;
            end
            
            % Decrease load step
            if subcount > 8
                NEWSTEPS = ceil((1-LOAD_SCALE)/(0.5*DSCALE));
                DSCALE = DSCALE * 0.5;
            end
            
            if subcount < 4 || subcount > 8
                NEWINTERVAL = (1-STEPS(1))/NEWSTEPS;
                STEPS       = linspace(STEPS(1) + NEWINTERVAL,1,NEWSTEPS);
                NEWSTEP(LOAD_COUNT) = STEPS(1);
            else
                NEWSTEP(LOAD_COUNT) = STEPS(1);
            end
            
            LOAD_SCALE = NEWSTEP(LOAD_COUNT);
            
        end
    else
        LOAD_SCALE = LOAD_COUNT * DSCALE;
    end
    %LOAD_SCALE = N * DSCALE;
    Fext = LOAD_SCALE .* F;
    DS_RES = 1;
    subcount = 0;
    NR0 = NR;
    CONVERGENCE = [];
    while DS_RES >RES_TOL
        subcount = subcount +1;
        
        F_flw = gf_flw_nodal(beam_model.Param.LOAD, beam_model.Info, ...
            beam_model.F_FLW, NR, beam_model.Node.DOF, LOAD_SCALE);
        
        if SimpleStiff == 1
            ik = [];
            jk = [];
            vk = [];
            
            [ibar, jbar, vbar] = set_DT_mat((1:beam_model.Info.nbar), beam_model.Bar, BARR, beam_model.Node.DOF,...
                beam_model.Node.Coord, NODEPOS, beam_model.Node.R, NR);
            
            ik = [ik;ibar];
            jk = [jk;jbar];
            vk = [vk;vbar];
            
            [ibeam, jbeam, vbeam] = set_DT_mat((1:beam_model.Info.nbeam), beam_model.Beam, BEAMR, beam_model.Node.DOF,...
                beam_model.Node.Coord, NODEPOS, beam_model.Node.R, NR);
            
            ik = [ik;ibeam];
            jk = [jk;jbeam];
            vk = [vk;vbeam];
            
            [ibush, jbush, vbush] = set_CBUSH_mat(beam_model.CBush, NR, beam_model.Res.NDispl(:,4:6),beam_model.Res.NRot);
            
            ik = [ik;ibush];
            jk = [jk;jbush];
            vk = [vk;vbush];
            
            Kbar = sparse(ik, jk, vk, ndof, ndof);
            
            J = Kbar ;
        else
            
            ik = [];
            jk = [];
            vk = [];
            % assembly BAR jacobian
            [ibar, jbar, vbar] = j_assembly(nbar, beam_model.Bar, BARR, ...
                beam_model.Node.Coord, beam_model.Node.R, NODEPOS, NR, ...
                beam_model.Node.DOF, beam_model.Res.Bar.CForces, BARF, ...
                beam_model.Res.NDispl,beam_model.F_FLW, LOAD_SCALE, ...
                beam_model.Param.LOAD, beam_model.Info);
            
            ik = [ik;ibar];
            jk = [jk;jbar];
            vk = [vk;vbar];
            
            [ibeam, jbeam, vbeam] = j_assembly(nbeam, beam_model.Beam, BEAMR, ...
                beam_model.Node.Coord, beam_model.Node.R, NODEPOS, NR, ...
                beam_model.Node.DOF, beam_model.Res.Beam.CForces, BEAMF, ...
                beam_model.Res.NDispl,beam_model.F_FLW, LOAD_SCALE, ...
                beam_model.Param.LOAD, beam_model.Info);
            
            ik = [ik;ibeam];
            jk = [jk;jbeam];
            vk = [vk;vbeam];
            
            Kbar = sparse(ik, jk, vk, ndof, ndof);
            
            J = Kbar ;
        end
        
        % assembly BAR jacobian
        
        if ~isempty(beam_model.RBE2)
            J = RBE2Assembly(beam_model.RBE2,J);
        end
        
        % element forces
        Fbar = get_actual_nodal_forces(ngrid, nbar, beam_model.Bar, BARR, ...
            beam_model.Node.DOF, NODEPOS, NR, beam_model.Res.Bar.CForces);
        
        Fbeam = get_actual_nodal_forces(ngrid, nbeam, beam_model.Beam, BEAMR, ...
            beam_model.Node.DOF, NODEPOS, NR, beam_model.Res.Beam.CForces);
        
        Fbush = get_bush_forces(ngrid,beam_model.CBush,beam_model.Res.NDispl, beam_model.Node, NR, beam_model.Res.NRot);
        
        FN = Fbar + Fbeam + Fbush;
        
        Fel = get_internal_forces(ndof, beam_model.Node.DOF, FN);
        
        R = Fext + F_flw - Fel;
        
        if ~isempty(beam_model.RBE2)
            Resid = RBE2Assembly2(beam_model.RBE2,R);
        else
            Resid=R;
        end
        dummyndof=size(Resid,1);
        
        %%%%%%%%%%%%%%%%%% HOLONOMIC CONSTRAINT GRADIENT MATRIX (RJOINT)%%%%%%%%%%%%%%
        if ~isempty(beam_model.RJoint.ID)
            % Describe a holonomic constraint run through the number of
            % Contraints
            
            B_joint  = zeros(TotConst,ndof2);
            count = 0;
            for j = 1:length(beam_model.RJoint.ID)
                nrdof = beam_model.RJoint.DOFC(j).list;
                n1 = beam_model.RJoint.Node(j,1);
                n2 = beam_model.RJoint.Node(j,2);
                for Bi = nrdof
                    count = count+1;
                    if Bi <= 3
                        B_joint(count,beam_model.Node.DOF2(n1,Bi)) = -1;
                        B_joint(count,beam_model.Node.DOF2(n2,Bi)) = 1;
                    elseif Bi == 5
                        % With respect to lamda 1
                        B_joint(count,beam_model.Node.DOF2(n1,4:6)) = (-(crossm(DR(:,:,n1)*JOINTR(:,1,1,j))*Gmat(gdef(n1,4:6)))'*DR(:,:,n2)*JOINTR(:,2,2,j))';
                        B_joint(count,beam_model.Node.DOF2(n2,4:6)) = (-(crossm(DR(:,:,n2)*JOINTR(:,2,2,j))*Gmat(gdef(n2,4:6)))'*DR(:,:,n1)*JOINTR(:,1,1,j))';
                    elseif Bi == 6
                        % With respect to lamda 2
                        B_joint(count,beam_model.Node.DOF2(n1,4:6)) = (-(crossm(DR(:,:,n1)*JOINTR(:,1,1,j))*Gmat(gdef(n1,4:6)))'*DR(:,:,n2)*JOINTR(:,3,2,j))';
                        B_joint(count,beam_model.Node.DOF2(n2,4:6)) = (-(crossm(DR(:,:,n2)*JOINTR(:,3,2,j))*Gmat(gdef(n2,4:6)))'*DR(:,:,n1)*JOINTR(:,1,1,j))';
                    end
                end
            end
            BlockZeros = zeros(TotConst,TotConst);
            Sk = [J,B_joint';B_joint,BlockZeros];
            
            % Calculate the difference between the 2 nodes of the hinge
            % Currently works on the displacement at that load step check
            % rotations here
            count = 0;
            for jk = 1:length(beam_model.RJoint.ID)
                nrdof = beam_model.RJoint.DOFC(jk).list;
                n1 = beam_model.RJoint.Node(jk,1);
                n2 = beam_model.RJoint.Node(jk,2);
                for Bi = nrdof
                    count = count+1;
                    if Bi <= 3
                        Rd(count,1) = -(gdef(n2,Bi) - gdef(n1,Bi));
                    elseif Bi == 5
                        %Rd(count,1) = -(DR(:,:,n1)*JOINTR(:,1,1,jk))'*(DR(:,:,n2)*JOINTR(:,2,2,jk));
                        Rd(count,1) = -JOINTR(:,1,1,jk)'*JOINTR(:,2,2,jk);
                    elseif Bi == 6
                        %Rd(count,1) = -(DR(:,:,n1)*JOINTR(:,1,1,jk))'*(DR(:,:,n2)*JOINTR(:,3,2,jk));
                        Rd(count,1) = -JOINTR(:,1,1,jk)'*JOINTR(:,3,2,jk);
                    end
                end
            end
            
            Resid = Resid - B_joint'*lambda_joint;
            
            Resid = [Resid;Rd];
            
        end
        
        %%%%%%%%%%%%%%%%%% HOLONOMIC CONSTRAINT GRADIENT MATRIX (RBAR) %%%%%%%%%%%%%%
        if ~isempty(beam_model.RBar.ID)
            % Describe a holonomic constraint run through the number of
            % Contraints
            
            B_rbar  = zeros(TotConst,ndof2);
            count = 0;
            for j = 1:length(beam_model.RBar.ID)
                nrdof = beam_model.RBar.DOFC(j).list;
                n1 = beam_model.RBar.Node(j,1);
                n2 = beam_model.RBar.Node(j,2);
                
                deltaR = DR(:,:,n1)*crossm(NR0(:,:,n1)*beam_model.RBar.rel(j).data')*Gmat(gdef(n1,4:6)');
                for Bi = nrdof
                    count = count+1;
                    if Bi <= 3 % Define
                        if beam_model.Node.DOF2(n1,Bi)~=0
                            B_rbar(count,beam_model.Node.DOF2(n1,Bi)) = -1;
                        end
                        if beam_model.Node.DOF2(n2,Bi)~=0
                            B_rbar(count,beam_model.Node.DOF2(n2,Bi)) = 1;
                        end
                        if beam_model.Node.DOF2(n1,4)~=0
                            B_rbar(count,beam_model.Node.DOF2(n1,4)) = deltaR(Bi,1);
                        end
                        if beam_model.Node.DOF2(n1,5)~=0
                            B_rbar(count,beam_model.Node.DOF2(n1,5)) = deltaR(Bi,2);
                        end
                        if beam_model.Node.DOF2(n1,6)~=0
                            B_rbar(count,beam_model.Node.DOF2(n1,6)) = deltaR(Bi,3);
                        end
                    elseif Bi >= 4 % Define the equality of rotational displacements
                        if beam_model.Node.DOF2(n1,Bi)~=0
                            B_rbar(count,beam_model.Node.DOF2(n1,Bi)) = -1;
                        end
                        if beam_model.Node.DOF2(n2,Bi)~=0
                            B_rbar(count,beam_model.Node.DOF2(n2,Bi)) = 1;
                        end
                    end
                end
            end
            BlockZeros = zeros(TotConst,TotConst);
            Sk = [J,B_rbar';B_rbar,BlockZeros];
            
            % Calculate the difference between the 2 nodes of the hinge
            % Currently works on the displacement at that load step check
            % rotations here
            count = 0;
            for jk = 1:length(beam_model.RBar.ID)
                nrdof = beam_model.RBar.DOFC(jk).list;
                n1 = beam_model.RBar.Node(jk,1);
                n2 = beam_model.RBar.Node(jk,2);
                Resid_bar = NODEPOS(n2,:)' - NODEPOS(n1,:)' - NR(:,:,n1)*beam_model.RBar.rel(jk).data';
                
                for Bi = nrdof
                    count = count+1;
                    if Bi <= 3
                        Rd(count,1) = -Resid_bar(Bi);
                    elseif Bi >= 4
                        Rd(count,1) = -(gdef(n2,Bi) - gdef(n1,Bi));
                    end
                end
            end
            
            Resid = Resid - B_rbar'*lambda_rbar;
            
            Resid = [Resid;Rd];
            
        end
        
        if ~isempty(beam_model.RJoint.ID)
            delta_q = Sk\Resid;
            SOL = zeros(dummyndof, 1);
            SOL(Jldof) = delta_q(1:dummyndof,1);
            delta_lambda_joint = delta_q(dummyndof+1:end,1);
            lambda_joint = lambda_joint + delta_lambda_joint;
            
        elseif ~isempty(beam_model.RBar.ID)
            delta_q = Sk\Resid;
            SOL = zeros(dummyndof, 1);
            SOL = delta_q(1:dummyndof,1);
            delta_lambda_rbar = delta_q(dummyndof+1:end,1);
            lambda_rbar = lambda_rbar + delta_lambda_rbar;
        else
            SOL = zeros(dummyndof, 1);
            SOL = J \ Resid;
        end
        
        if ~isempty(beam_model.RBE2)
            SOL = RBE2disp(beam_model.RBE2,SOL,ndof);
        end
        
        gdef = zeros(beam_model.Info.ngrid, 6);
        DR = zeros(3,3,beam_model.Info.ngrid);
        for n = 1:beam_model.Info.ngrid
            dof = beam_model.Node.DOF(n, 1:6);
            index = find(dof);
            if ~isempty(index)
                gdef(n, index) = SOL(dof(index));
            end
            DR(:,:,n) = Rmat(gdef(n, 4:6)');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% UPDATE ROTATION MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%
        NR0 = NR;
        NR     = update_node_rot(ngrid, DR, NR);
        BARR   = update_bar_rot(nbar, DR, beam_model.Bar.Conn, BARR, gdef);
        BEAMR  = update_bar_rot(nbeam, DR, beam_model.Beam.Conn, BEAMR, gdef);
        BUSHR  = update_bush_rot(beam_model.Info.ncbush, beam_model.CBush, DR, BUSHR);
        JOINTR = update_bush_rot(beam_model.Info.nrjoint, beam_model.RJoint, DR, JOINTR);
        
        beam_model.Res.NRd    = NR;
        beam_model.Res.Bar.R  = BARR;
        beam_model.Res.Beam.R = BEAMR;
        
        %%%%%%%%%%%%%%%%%%%%%% UPDATE AERO NODE POSITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
        %% This is needed to match SOL 106
        if (beam_model.Info.nrbe0 > 0)
            AERO_POS = update_aerobeam_node_nl(beam_model.Info.ngrid, beam_model.Node, gdef(:,1:3,1),...
                DR-repmat(eye(3,3),[1,1,beam_model.Info.ngrid]),NR0);
            % update coord database with slave nodes position
            for n=1:beam_model.Info.ngrid
                ne = length(beam_model.Node.Aero.Index(n).data);
                if ne
                    gdef(beam_model.Node.Aero.Index(n).data, 1:3, 1) = AERO_POS(n).data';
                end
            end
            clear AERO_POS;
        end
        
        beam_model.Res.NDispl = beam_model.Res.NDispl + gdef;
        beam_model.Res.NRot   = beam_model.Res.NRot + gdef(:,4:6);
        NODEPOS = beam_model.Node.Coord(:, 1:3) + beam_model.Res.NDispl(:,1:3);
        
        %%%%%%%%%%%%%%%%%%%%%% RECALCULATE INTERNAL LOADS %%%%%%%%%%%%%%%%%%%%%%%%%
        [beam_model.Res.Bar.CStrains, beam_model.Res.Bar.CForces, ...
            beam_model.Res.Bar.CStresses, beam_model.Res.Bar.CSM] = ...
            get_bar_force_strain_NL(nbar, beam_model.Bar, beam_model.PBar, ...
            beam_model.Mat, beam_model.Node, beam_model.Res.NDispl, NR, BARR, BARPO, BARKR);
        
        [beam_model.Res.Beam.CStrains, beam_model.Res.Beam.CForces, ...
            beam_model.Res.Beam.CStresses, beam_model.Res.Beam.CSM] = ...
            get_bar_force_strain_NL(nbeam, beam_model.Beam, beam_model.PBeam, ...
            beam_model.Mat, beam_model.Node, beam_model.Res.NDispl, NR, BEAMR, BEAMPO, BEAMKR);
        
        DS_RES = sqrt(sum(sum(gdef(:,1:3) .* gdef(:,1:3),2),1));
        F_RES  = sqrt(sum(Resid.*Resid));
        
        %         CONVERGENCE = [CONVERGENCE,DS_RES];
        %
        %         figure(100);
        %         semilogy(1:subcount,CONVERGENCE);
        
        if (DS_RES < RES_TOL)
            break;
        end
        
        if subcount >= MaxIter
            warning(['Iteration Count has been exceeded at Loadstep: ' num2str(N) '. Try increasing the number of loadsteps']);
            break;
        end
    end
    snapnow;
    TipDisp(LOAD_COUNT+1) = beam_model.Res.NDispl(WngNodes(end));
    
    % store curvatures
    BARKR  = beam_model.Res.Bar.CStrains(:, 4:6, :);
    BEAMKR = beam_model.Res.Beam.CStrains(:, 4:6, :);
    
    % store global element forces
    BARF  =  store_global_elem_force(nbar, BARR, beam_model.Res.Bar.CForces);
    BEAMF =  store_global_elem_force(nbeam,BEAMR,beam_model.Res.Beam.CForces);
    
    if ENABLE_DISP==1
        fprintf(['\n\t-Loadstep #' num2str(LOAD_COUNT) '...done in ' num2str(subcount) ' iterations']);
    end
    beam_model.Res.Bar.Colloc = bar_defo_colloc(nbar, beam_model.Bar, beam_model.Node.DOF, NODEPOS, beam_model.Res.NRd);
    beam_model.Res.Beam.Colloc = bar_defo_colloc(nbeam, beam_model.Beam, beam_model.Node.DOF, NODEPOS, beam_model.Res.NRd);
    
    %     bary = [];for i = 1:beam_model.WingNodes(end)-1;bary = [bary;beam_model.Res.Bar.Colloc(:,2,i)];end
    %     barstrain = [];for i = 1:beam_model.WingNodes(end)-1;barstrain = [barstrain;beam_model.Res.Bar.CStrains(:,6,i)];end
    %     figure(700);hold on;plot(bary,barstrain,'ko-');
end

%%%%%%%%%%%%%%%%%%% RECALCULATE INTERNAL LOADS ALONG THE BEAM AXIS  %%%%%%%%%%%%%%%%%%%%%%
% Recalculate the beam rotation matrix along the beam axis:

%[~,beam_model.Res.Beam.CForces] = BeamRotationMatrix(beam_model);


% %This is the correct way of updating the RBE2s (SOL400)
% if (beam_model.Info.nrbe0 > 0)
%     AERO_POS = update_aerobeam_node(beam_model.Info.ngrid, beam_model.Node, beam_model.Res.NDispl(:,1:3,1),...
%         NR-repmat(eye(3,3),[1,1,beam_model.Info.ngrid]));
%     %            update coord database with slave nodes position
%     for n=1:beam_model.Info.ngrid
%         ne = length(beam_model.Node.Aero.Index(n).data);
%         if ne
%             Aerodef(beam_model.Node.Aero.Index(n).data, 1:3, 1) = AERO_POS(n).data';
%         end
%     end
%     clear AERO_POS;
% end
%
% beam_model.Res.NDispl(:,1:3) = beam_model.Res.NDispl(:,1:3) + Aerodef(:,1:3,1);
%
% % Total Lagrangian node coords
% NODEPOS = beam_model.Node.Coord(:, 1:3) + beam_model.Res.NDispl(:,1:3);

% Storing the internal loads at the midpoint of each of the beams
Fx = [];
Fy = [];
Fz = [];
Mx = [];
My = [];
Mz = [];
if beam_model.Info.nbar>0
    for j = 1:size(beam_model.Bar.Colloc,3)
        Fx    = [Fx;mean([beam_model.Res.Bar.CForces(1,1,j),beam_model.Res.Bar.CForces(2,1,j)])];
        Fy    = [Fy;mean([beam_model.Res.Bar.CForces(1,2,j),beam_model.Res.Bar.CForces(2,2,j)])];
        Fz    = [Fz;mean([beam_model.Res.Bar.CForces(1,3,j),beam_model.Res.Bar.CForces(2,3,j)])];
        Mx    = [Mx;mean([beam_model.Res.Bar.CForces(1,4,j),beam_model.Res.Bar.CForces(2,4,j)])];
        My    = [My;mean([beam_model.Res.Bar.CForces(1,5,j),beam_model.Res.Bar.CForces(2,5,j)])];
        Mz    = [Mz;mean([beam_model.Res.Bar.CForces(1,6,j),beam_model.Res.Bar.CForces(2,6,j)])];
    end
end
if beam_model.Info.nbeam>0
    for j = 1:size(beam_model.Beam.Colloc,3)
        Fx    = [Fx;mean([beam_model.Res.Beam.CForces(1,1,j),beam_model.Res.Beam.CForces(2,1,j)])];
        Fy    = [Fy;mean([beam_model.Res.Beam.CForces(1,2,j),beam_model.Res.Beam.CForces(2,2,j)])];
        Fz    = [Fz;mean([beam_model.Res.Beam.CForces(1,3,j),beam_model.Res.Beam.CForces(2,3,j)])];
        Mx    = [Mx;mean([beam_model.Res.Beam.CForces(1,4,j),beam_model.Res.Beam.CForces(2,4,j)])];
        My    = [My;mean([beam_model.Res.Beam.CForces(1,5,j),beam_model.Res.Beam.CForces(2,5,j)])];
        Mz    = [Mz;mean([beam_model.Res.Beam.CForces(1,6,j),beam_model.Res.Beam.CForces(2,6,j)])];
    end
end

BarForces.Fx = Fx;
BarForces.Fy = Fy;
BarForces.Fz = Fz;
BarForces.Mx = Mx;
BarForces.My = My;
BarForces.Mz = Mz;

%beam_model.Loads.nlinstatic = storebarloads(beam_model);
% figure(503);
% hold off;
% plot(beam_model.Node.Coord(WngNodes,2) + beam_model.Res.NDispl(WngNodes,2),beam_model.Res.NDispl(WngNodes,1),'b^--');
% hold on;
% plot(beam_model.Node.Coord(WngNodes,2) + beam_model.Res.NDispl(WngNodes,2),beam_model.Res.NDispl(WngNodes,2),'bo--');
% hold on;
% plot(beam_model.Node.Coord(WngNodes,2) + beam_model.Res.NDispl(WngNodes,2),beam_model.Res.NDispl(WngNodes,3),'b+--');
% legend('x disp','y disp','z disp','Location','Best');

beam_model.Res.SOL = SOL;
beam_model.Res.K = J;

beam_model.Res.Bar.F0  = BARF;
beam_model.Res.Beam.F0 = BEAMF;
beam_model.Res.Bar.K0  = BARKR;
beam_model.Res.Beam.K0 = BEAMKR;

end
