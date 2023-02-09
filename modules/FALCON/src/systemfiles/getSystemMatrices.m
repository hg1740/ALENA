function SYS = getSystemMatrices(X, Matrices, Aero, Sim, AnalysisParam)

%Grab variables for shorter indexing
numElem = Matrices.n_elem;    %How many elements in each component?
numNode = Matrices.n_node;    %How many structural nodes in each component?
nPart   = numel(numElem);     %How many components?

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %Get the coordinates and orientations at each node
        [coords, CaB] = strains2coords_all(X.x_f_p1, Matrices);
            
        %Preallocate
        %   - System global stiffness matrix (e.g. [Atot] * {X} - {f} = 0)
        Atot      = arrayfun(@(i) zeros(numNode(i) * 6, numElem(i) * 6), 1 : nPart, 'Unif', false);
        %   - Derivative of system w.r.t local loads (states)
        %   e.g. dQstiftot = d[Atot]/d{X}
        dQstiftot = arrayfun(@(i) zeros(numNode(i) * 6, numElem(i) * 6), 1 : nPart, 'Unif', false);
        %   - Loads
        f           = arrayfun(@(nE) zeros(nE * 6, 1), numElem, 'Unif', false);
        f_grav      = arrayfun(@(nE) zeros(nE * 6, 1), numElem, 'Unif', false);
        f_grav_glob = arrayfun(@(nE) zeros(nE * 6, 1), numElem, 'Unif', false);
        AoA_qs_t    = arrayfun(@(nE) zeros(nE, 1)    , numElem, 'Unif', false); 
        Btot_r      = arrayfun(@(nE) zeros(6, nE * 6), numElem, 'Unif', false); 
        f_glob      = arrayfun(@(nN) zeros(nN * 6, 1), numNode, 'Unif', false);
        f_grav_pt   = arrayfun(@(nN) zeros(nN * 6, 1), numNode, 'Unif', false);
        F_r         = arrayfun(@(nN) zeros(6, nN * 6), numNode, 'Unif', false);         
        
        %Loop through parts and return the system matrices
        for jj = 1 : nPart
            
            %Grab part data
            nElem = numElem(jj);
            nNode = numNode(jj);
            
            for ii = 1 : nElem %For each element...
                
                %Index numbers for this element
                ind2 = (1 : 12) + ((ii - 1) * 6);
                ind3 = (1 : 6)  + ((ii - 1) * 6);
                
                %Get the element coordinate system
                if Matrices.Parent{jj}(1)
                    CGB_e = Matrices.CGa*CaB{Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii));
                else
                    CGB_e = Matrices.CGa*CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii));
                end
                
                %Grab the element loads
                %   - X.x_f_pl is the vector of nodal/elemental forces
                elemLoads = X.x_f_p1{jj}(ind3);
                
                %% Structural                
                
                %Grab data
                complianceMatrix = Matrices.Aa{jj}(:, :, ii);    %Compliance matrix without shape function terms
                e1               = Matrices.e1tot{jj}(ind3);
                Ab               = Matrices.Ab{jj}(:,:,ii);      %Ab = (A_1)' * dS + (A_2)' * (dS^2 / 2)
                
                %'f2ftilde' returns the matrix where F & M are the beam
                %forces and moments in the element ref. frame
                %
                % [0      , F_tilde ;
                %  F_tilde, M_tilde]
                xftilde     = f2ftilde(elemLoads);
                
                %Calculates the derivative of xftilde in order to pass to the
                %Newton-Raphson solver. This is:
                %
                % [C] * [F + eF ]
                %       [  M    ]
                %
                xftilde_var = v2vtilde(complianceMatrix * (elemLoads + e1));
                
                % A = [C] * [0 F_tilde ; F_tilde, M]
                A      =  Ab * xftilde  * complianceMatrix;
                % dQ = [C] * (d/dM * [0 F_tilde ; F_tilde, M] + d/dF * [0
                % F_tilde ; F_tilde, M]) (Partial derivative)
                dQstif =  Ab * xftilde_var;
                
                %Create matrix for component
                Atot{jj}(ind2,ind3)      = Atot{jj}(ind2,ind3)      + A;
                dQstiftot{jj}(ind2,ind3) = dQstiftot{jj}(ind2,ind3) + dQstif;
                
                %% Aero
                if Sim.aero_flag && Matrices.isLS{jj}
                    %                 if Aero.inflow_switch == 1%iter_counter == 0 &&
                    %                     gust_vec{jj}(ind4,:) = getInflow;
                    %                 elseif Aero.inflow_switch == 0
                    %                     gust_vec{jj}(ind4,:) = [0;0;0];
                    %                 end
                    %
                    %                 CA0B_tot_blkdiag = [Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)); zeros(3) eye(3)];
                    %                 CBA0_tot_blkdiag = CA0B_tot_blkdiag';
                    %
                    %                 v_A0     = CA0B_tot_blkdiag*- [CGB_e' zeros(3); zeros(3) zeros(3)]*[gust_vec{jj}(ind4);0;0;0];
                    %
                    %                 ydot =  v_A0(2);
                    %                 zdot =  v_A0(3);
                    
                    %                 x_v = [CaB(:,:,Matrices.s_mp_ind{jj}(ii))',zeros(3);zeros(3),CaB(:,:,Matrices.s_mp_ind{jj}(ii))']*X.x_va_p1;
                    %
                    %                 CA0B_tot_blkdiag = [Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)); zeros(3) eye(3)];
                    %                 CBA0_tot_blkdiag = CA0B_tot_blkdiag';
                    %
                    %                 v_A0     = CA0B_tot_blkdiag*x_v;
                    %
                    %                 ydot =  v_A0(2);
                    %                 zdot =  v_A0(3);
                    
                    CA0B_tot_blkdiag = [Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)); zeros(3) eye(3)];
                    CBA0_tot_blkdiag = CA0B_tot_blkdiag';
                    
                    v_A0     = CA0B_tot_blkdiag*[CGB_e' zeros(3); zeros(3) zeros(3)]*[Matrices.CGa,zeros(3);zeros(3),Matrices.CGa]*X.x_va_p1;
                    
                    ydot =  v_A0(2);
                    zdot =  v_A0(3);
                    
                    tiplossFactor = (Aero.lift_dist{jj}(ii)+Aero.lift_dist{jj}(ii+1))/2;
                    chord         = Aero.chord{jj}(ii);
                    
                    AoA_qs_t{jj}(ii,1) = - zdot/ydot;
                    
                    if 0%AoA_qs_t{jj}(ii,1) > Aero.stall_angle/180*pi
                        Lift               = pi*Aero.rho*ydot^2*Aero.stall_angle/180*pi*tiplossFactor*chord;
                        AoA_qs_t{jj}(ii,1) = Aero.stall_angle/180*pi;
                    else
                        Lift               =  - pi*Aero.rho*ydot*zdot*tiplossFactor*chord - 0.5*Aero.rho*ydot*ydot*Aero.dCLdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor*chord;
                    end
                    Drag           = Lift*AoA_qs_t{jj}(ii,1);
                    Mom            = 0.5*Aero.rho*ydot*ydot*Aero.dCMdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor*chord*chord;
                    aeroForce      = CBA0_tot_blkdiag*[0;Drag;Lift;Mom;0;0];
                    
                    f{jj}(ind3,1)  = aeroForce;
                else
                    f{jj}(ind3,1)      = zeros(6,1);
                    AoA_qs_t{jj}(ii,1) = 0;
                end
                %%
                %Calculate the position and orientation of any point where two
                %beams are connected (parent/child)
                if Matrices.Parent{jj}(1)
                    Ra    =  coords{Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2))) + ...
                        CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*(coords{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) + Matrices.p0{jj});
                    CaB_d = [CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) zeros(3); ...
                        zeros(3) CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii))];
                else
                    Ra    = coords{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) + Matrices.p0{jj};
                    CaB_d = [CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) zeros(3); zeros(3) CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii))];
                end
                Ra_t1 = [eye(3)  ,zeros(3);skew(Ra),eye(3)  ];
                %             Ra    = coords(:,:,Matrices.s_mp_ind{jj}(ii)) + Matrices.p0{jj} + CaB(:,:,Matrices.s_mp_ind{jj}(ii))*Matrices.cp{jj}(:,ii);
                %             Ra_t2 = [eye(3)  ,zeros(3);skew(Ra),eye(3)  ];
                
                %Calculates the global forces about the origin of the aircraft
                %frame - assumes the force is constant over the element
                Btot_r{jj}(:,ind3)  = Ra_t1 * CaB_d * (Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii));
                %             Btot_r_aero{jj}(:,ind3)  = Btot_r_grav{jj}(:,ind3);%Ra_t2 * CaB_d * (Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii));
                %% Gravity
                if Sim.grav_fact
                    %Aa2 is the mass matrix
                    gravForce               = Matrices.Aa2{jj}(:,:,ii)*[CGB_e' zeros(3); zeros(3) CGB_e']*Sim.grav_vec;
                    f_grav{jj}(ind3,1)      = gravForce;
                    f_grav_glob{jj}(ind3,1) = [CGB_e' zeros(3); zeros(3) CGB_e']*gravForce;
                else
                    f_grav{jj}(ind3,1)      = zeros(6,1);
                    f_grav_glob{jj}(ind3,1) = zeros(6,1);
                end
                %%
            end
            
            %%
            for ii = 1:nNode
                ind3  = [1:6]  + (ii-1)*6;
                if Matrices.Parent{jj}(1)
                    CGB_e = Matrices.CGa*CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii));
                    Ra    =  coords{Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2))) + ...
                        CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*(coords{jj}(:,:,Matrices.s_node_ind{jj}(ii)) + Matrices.p0{jj});
                    CaB_d = [CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii)) zeros(3); ...
                        zeros(3) CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}'*CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii))];
                else
                    CGB_e = Matrices.CGa*CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii));
                    %Nodal positions
                    Ra    = coords{jj}(:,:,Matrices.s_node_ind{jj}(ii)) + Matrices.p0{jj};
                    %Nodal orientations
                    CaB_d = [CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii)) zeros(3); zeros(3) CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii))];
                end
                %f_glob are vector of loads in global coordinate system which
                %have now been transformed into element system
                %f_a are loads in aircraft frame
                f_glob{jj}(ind3,1) = [CGB_e',zeros(3);zeros(3),CGB_e']*Matrices.f_a{jj}(ind3);
                
                Ra_t    = [eye(3)  ,zeros(3);skew(Ra),eye(3)  ];
                F_r{jj}(:,ind3)  = Ra_t * CaB_d ; % loads at the nodes
                %             F_rg{jj}(:,ind3) = Ra_t;
                
                %Loads from point masses
                gravForce_pt               = Matrices.M_pt{jj}(:,:,ii)*[CGB_e' zeros(3); zeros(3) CGB_e']*Sim.grav_vec;
                f_grav_pt{jj}(ind3,1)      = gravForce_pt;
                %             f_grav_pt_glob{jj}(ind3,1) = [CGB_e' zeros(3); zeros(3) CGB_e']*gravForce_pt;
            end
            
        end
        
        %Set up sparse matrices
        Atot      = cellfun(@sparse, Atot     , 'Unif', false);
        dQstiftot = cellfun(@sparse, dQstiftot, 'Unif', false);
        
        %Stash
        SYS.AoA_qs_t  = AoA_qs_t;
        SYS.Btot_r    = Btot_r  ;
        SYS.F_r       = F_r;
        SYS.Atot      = Atot;
        SYS.dQstiftot = dQstiftot;
        SYS.f         = f;
        SYS.f_glob    = f_glob;
        SYS.f_grav    = f_grav;
        SYS.f_grav_pt = f_grav_pt;
        
    otherwise
        
        if Sim.rb_flag
            pa    = X.x_va_p1(1:3);
            omega = X.x_va_p1(4:6);
            quat  = X.x_qa_p1;
            
            OmegaQuat_a  = [ 0          omega(1)  omega(2)  omega(3);
                -omega(1)   0        -omega(3)  omega(2);
                -omega(2)   omega(3)  0        -omega(1);
                -omega(3)  -omega(2)  omega(1)  0       ];
            
            CGa_e   = Quat2Rot(X.x_qa_p1);
            CGa_tot  = [CGa_e,zeros(3)];
            
            if ~Sim.speed_up
                dOmegaQuat_e = [ quat(2)  quat(3)  quat(4);
                    -quat(1)  quat(4) -quat(3);
                    -quat(4) -quat(1)  quat(2);
                    quat(3) -quat(2) -quat(1)];
                
                dOmegaQuat_a = [zeros(4,3),dOmegaQuat_e];
                
                dCGa_e1 = 2*[ quat(1) -quat(4)  quat(3);
                    quat(4)  quat(1) -quat(2);
                    -quat(3)  quat(2)  quat(1)];
                
                dCGa_e2 = 2*[ quat(2)  quat(3)  quat(4);
                    quat(3) -quat(2) -quat(1);
                    quat(4)  quat(1) -quat(2)];
                
                dCGa_e3 = 2*[-quat(3)  quat(2)  quat(1);
                    quat(2)  quat(3)  quat(4);
                    -quat(1)  quat(4) -quat(3)];
                
                dCGa_e4 = 2*[-quat(4) -quat(1)  quat(2);
                    quat(1) -quat(4)  quat(3);
                    quat(2)  quat(3)  quat(4)];
                
                dCaG_e1_blkdiag = [dCGa_e1' zeros(3);zeros(3) dCGa_e1'];
                dCaG_e2_blkdiag = [dCGa_e2' zeros(3);zeros(3) dCGa_e2'];
                dCaG_e3_blkdiag = [dCGa_e3' zeros(3);zeros(3) dCGa_e3'];
                dCaG_e4_blkdiag = [dCGa_e4' zeros(3);zeros(3) dCGa_e4'];
                
                dCGa_tot = [dCGa_e1*pa dCGa_e2*pa dCGa_e3*pa dCGa_e4*pa];
            else
                dOmegaQuat_a = zeros(4,6);
                dCGa_tot     = zeros(3,4);
            end
            
            omega_tilde   = [skew(omega) zeros(3);skew(pa) skew(omega)];
            PH            = Matrices.M_rb*X.x_va_p1;
            domega_tilde  = [zeros(3) -skew(PH(1:3));-skew(PH(1:3)) -skew(PH(4:6))];
            
            if Sim.grav_fact
                gravForce     = Matrices.M_rb*[CGa_e' zeros(3); zeros(3) CGa_e']*Sim.grav_vec;
                f_grav_r      = gravForce;
                f_grav_glob_r = [CGa_e zeros(3); zeros(3) CGa_e]*gravForce;
                
                if ~Sim.speed_up
                    dF_grav_tot_r = [Matrices.M_rb*dCaG_e1_blkdiag*Sim.grav_vec...
                        Matrices.M_rb*dCaG_e2_blkdiag*Sim.grav_vec...
                        Matrices.M_rb*dCaG_e3_blkdiag*Sim.grav_vec...
                        Matrices.M_rb*dCaG_e4_blkdiag*Sim.grav_vec];
                else
                    dF_grav_tot_r = zeros(6,4);
                end
            else
                f_grav_r      = zeros(6,1);
                dF_grav_tot_r = zeros(6,4);
            end
        end
        
        if ~Sim.speed_up
            dQgyr_r_dpa = zeros(6,3);
            dQgyr_r_dqa = zeros(6,4);
            dF_grav_dpa = zeros(6,3);
            dF_grav_dqa = zeros(6,4);
            dF_tot_dpa  = zeros(6,3);
            dF_tot_dqa  = zeros(6,4);
        end
        
        for jj = 1:nPart
            Atot{jj}        = sparse(zeros(nNode*6,nElem*6));
            Ctot{jj}        = sparse(zeros(nNode*6,nElem*6));
            Etot{jj}        = sparse(zeros(nNode*6,nElem*6));
            
            Qgyr_r{jj}      = sparse(zeros(                        6,nElem*6));
            M_r{jj}         = sparse(zeros(                        6,nElem*6));
            
            F_r{jj}         = sparse(zeros(                        6,nElem*6+6));
            F_rg{jj}        = sparse(zeros(                        6,nElem*6+6));
            
            dQgyr_r{jj}     = sparse(zeros(                        6,nElem*6));
            dQgyr_r_dv{jj}  = sparse(zeros(                        6,                    6));
            
            dQstiftot{jj}   = sparse(zeros(nNode*6,nElem*6));
            dQgyrtot{jj}    = sparse(zeros(nNode*6,nElem*6));
            dTgamtot{jj}    = sparse(zeros(nNode*6,nElem*6));
            
            OmegaQuat{jj}   = sparse(zeros( nElem   *4));
            dOmegaQuat{jj}  = sparse(zeros( nElem   *4,nElem*6));
            
            CGB_tot{jj}     = sparse(zeros( nElem   *3,nElem*6));
            dCGB_tot{jj}    = sparse(zeros( nElem   *3,nElem*4));
            
            if Sim.aero_flag
                A_x{jj}         = zeros(nElem*2,nElem*2);
                A_x_1{jj}       = zeros(nElem*2,1);
                
                dA_x{jj}        = zeros(nElem*2,nElem*2);
                dA_x_dwg{jj}    = zeros(nElem*2,nElem*3);
                dA_x_1{jj}      = zeros(nElem*2,nElem*6);
                
                dF_x{jj}        = zeros(nElem*6,nElem*2);
                
                dA_x_1_dwg{jj}  = zeros(nElem*2,nElem*3);
            else
                A_x{jj}         = [];
                A_x_1{jj}       = zeros(0,1);
                
                dA_x{jj}        = [];
                dA_x_dwg{jj}    = [];
                dA_x_1{jj}      = [];
                
                dF_x{jj}        = zeros(nElem*6,0);
                
                dF_wg_tot       = [];
                dF_u_tot        = [];
                dA_x_1_dwg      = [];
            end
            
            dF_grav_tot{jj} = zeros(nElem*6,nElem*4);
            dF_tot{jj}      = zeros(nElem*6,nElem*6);
            
            Btot_r{jj}      = zeros(6                      ,nElem*6);
            dV_dva{jj}      = zeros(nElem*6,                      6);
            
            AoA_qs_t{jj}    = zeros(nElem,1);
            
            f{jj}           = zeros(nElem,1);
            
            CGa = Quat2Rot(X.x_qa_p1);
            
            for ii = 1:nElem
                %% Structural
                ind1 = [1:4]  + (ii-1)*4;
                ind2 = [1:12] + (ii-1)*6;
                ind3 = [1:6]  + (ii-1)*6;
                ind4 = [1:3]  + (ii-1)*3;
                ind5 = [1:2]  + (ii-1)*2;
                
                if 0
                    xftilde = f2ftilde(X.x_f_p1{jj}(ind3),0);
                    xvtilde = v2vtilde(X.x_v_p1{jj}(ind3),0);
                    
                    xftilde_var = v2vtilde(Matrices.Aa{jj}(:,:,ii) *(X.x_f_p1{jj}(ind3)+Matrices.e1tot{jj}(ind3)),0);
                    xvtilde_var = f2ftilde(Matrices.Aa2{jj}(:,:,ii)* X.x_v_p1{jj}(ind3),0);
                    
                    A  =  Matrices.Ab{jj}(:,:,ii) * xftilde  * Matrices.Aa{ jj}(:,:,ii);
                    C  =  Matrices.Ab{jj}(:,:,ii) * xvtilde  * Matrices.Aa2{jj}(:,:,ii);
                    E  = -Matrices.Ab{jj}(:,:,ii) * xvtilde' * Matrices.Aa{ jj}(:,:,ii);
                    
                    dQstif =  Matrices.Ab{jj}(:,:,ii) * xftilde_var;
                    dQgyr  =  Matrices.Ab{jj}(:,:,ii) * xvtilde_var;
                    dTgam  = -Matrices.Ab{jj}(:,:,ii) * xftilde_var';
                else
                    IkronVO     =  kron(eye(6),X.x_v_p1{jj}(ind3));
                    IkronMF     =  kron(eye(6),X.x_f_p1{jj}(ind3));
                    
                    A           =  Matrices.A_pre1{jj}(:,:,ii)*IkronMF;
                    C           =  Matrices.C_pre1{jj}(:,:,ii)*IkronVO;
                    E           =  Matrices.E_pre1{jj}(:,:,ii)*IkronVO;
                    
                    C_pt        =  Matrices.C_pt_pre1{jj}(:,:,ii)*IkronVO;
                    
                    dQstif      = -Matrices.A_pre2{jj}(:,:,ii)*IkronMF - Matrices.A_pre3{jj}(:,:,ii);
                    dQgyr       = -Matrices.C_pre2{jj}(:,:,ii)*IkronVO;
                    dTgam       = -Matrices.E_pre2{jj}(:,:,ii)*IkronMF - Matrices.E_pre3{jj}(:,:,ii);
                    
                    dQgyr_pt    = -Matrices.C_pt_pre2{jj}(:,:,ii)*IkronVO;
                    
                    xvtilde     =  Matrices.E1tilde*IkronVO;
                    xvtilde_var =  Matrices.E2tildeIkronM{jj}(:,:,ii)*IkronVO;
                end
                
                Atot{jj}(ind2,ind3)       = Atot{jj}(ind2,ind3)       + A;
                Ctot{jj}(ind2,ind3)       = Ctot{jj}(ind2,ind3)       + C;
                Etot{jj}(ind2,ind3)       = Etot{jj}(ind2,ind3)       + E;
                
                Ctot{jj}(ind3+6,ind3)     = Ctot{jj}(ind3+6,ind3)     + C_pt;
                
                dQstiftot{jj}(ind2,ind3)  = dQstiftot{jj}(ind2,ind3)  + dQstif;
                dQgyrtot{jj}(ind2,ind3)   = dQgyrtot{jj}(ind2,ind3)   + dQgyr;
                dTgamtot{jj}(ind2,ind3)   = dTgamtot{jj}(ind2,ind3)   + dTgam;
                
                dQgyrtot{jj}(ind3+6,ind3) = dQgyrtot{jj}(ind3+6,ind3) + dQgyr_pt;
                
                pa    = X.x_v_p1{jj}(ind3(1:3));%PosOrientNodal(ind3(1:3));
                omega = X.x_v_p1{jj}(ind3(4:6));%PosOrientNodal(ind3(4:6));
                
                quat  = X.x_q_p1{jj}(ind1);
                
                OmegaQuat_e  = [ 0          omega(1)  omega(2)  omega(3);
                    -omega(1)   0        -omega(3)  omega(2);
                    -omega(2)   omega(3)  0        -omega(1);
                    -omega(3)  -omega(2)  omega(1)  0       ];
                
                OmegaQuat{jj}(ind1,ind1)  = OmegaQuat_e;
                
                CGB_e   = Quat2Rot(quat);
                
                CGB_tot{jj}(ind4,ind3)  = [CGB_e,zeros(3)];
                
                if ~Sim.speed_up
                    dOmegaQuat_e = [ quat(2)  quat(3)  quat(4);
                        -quat(1)  quat(4) -quat(3);
                        -quat(4) -quat(1)  quat(2);
                        quat(3) -quat(2) -quat(1)];
                    
                    dOmegaQuat{jj}(ind1,ind3) = [zeros(4,3),dOmegaQuat_e];
                    
                    dCGB_e1 = 2*[ quat(1) -quat(4)  quat(3);
                        quat(4)  quat(1) -quat(2);
                        -quat(3)  quat(2)  quat(1)];
                    
                    dCGB_e2 = 2*[ quat(2)  quat(3)  quat(4);
                        quat(3) -quat(2) -quat(1);
                        quat(4)  quat(1) -quat(2)];
                    
                    dCGB_e3 = 2*[-quat(3)  quat(2)  quat(1);
                        quat(2)  quat(3)  quat(4);
                        -quat(1)  quat(4) -quat(3)];
                    
                    dCGB_e4 = 2*[-quat(4) -quat(1)  quat(2);
                        quat(1) -quat(4)  quat(3);
                        quat(2)  quat(3)  quat(4)];
                    
                    dCBG_e1_blkdiag = [dCGB_e1' zeros(3);zeros(3) dCGB_e1'];
                    dCBG_e2_blkdiag = [dCGB_e2' zeros(3);zeros(3) dCGB_e2'];
                    dCBG_e3_blkdiag = [dCGB_e3' zeros(3);zeros(3) dCGB_e3'];
                    dCBG_e4_blkdiag = [dCGB_e4' zeros(3);zeros(3) dCGB_e4'];
                    
                    dCGB_tot{jj}(ind4,ind1) = [dCGB_e1*pa dCGB_e2*pa dCGB_e3*pa dCGB_e4*pa];
                else
                    dOmegaQuat{jj}(ind1,ind3) = zeros(4,6);
                    dCGB_tot{jj}(ind4,ind1)   = zeros(3,4);
                end
                
                if Sim.rb_flag || Sim.rbloads_flag
                    CaB = CGa'*CGB_e;
                    
                    Ra        = CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1);
                    Ra_node   = CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*[1;0;0]*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))/2);
                    CaB_d     = [CaB     ,zeros(3);zeros(3),     CaB     ];
                    Ra_t      = [eye(3)  ,zeros(3);skew(Ra),     eye(3)  ];
                    Ra_t_node = [eye(3)  ,zeros(3);skew(Ra_node),eye(3)  ];
                    
                    M_r{jj}(:,ind3)     = (Ra_t*      CaB_d*        Matrices.Aa2{ jj}(:,:,ii  ))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)) + ...
                        (Ra_t_node* CaB_d*        Matrices.M_pt{jj}(:,:,ii+1))* Matrices.T_pt{jj}(:,:,ii);
                    
                    Qgyr_r{jj}(:,ind3)  = (Ra_t*      CaB_d*xvtilde*Matrices.Aa2{      jj}(:,:,ii))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)) + ...
                        (Ra_t_node* CaB_d*        Matrices.C_pt_pre1{jj}(:,:,ii))*IkronVO;
                    
                    dQgyr_r{jj}(:,ind3) = -Ra_t*      CaB_d*xvtilde_var                            *(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)) + ...
                        +Ra_t_node* CaB_d*        Matrices.C_pt_pre2{jj}(:,:,ii) *IkronVO;
                    
                    dV_dva{jj}(ind3,:)  = CaB_d'*Ra_t';
                    
                    Btot_r{jj}(:,ind3)  = Ra_t * CaB_d * (Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii));
                    
                    if ~Sim.speed_up
                        dRa1          = [1;0;0];
                        dRa_t1        = [zeros(3),zeros(3);skew(dRa1),zeros(3)];
                        dRa2          = [0;1;0];
                        dRa_t2        = [zeros(3),zeros(3);skew(dRa2),zeros(3)];
                        dRa3          = [0;0;1];
                        dRa_t3        = [zeros(3),zeros(3);skew(dRa3),zeros(3)];
                        CaG_d         = [CGa'    ,zeros(3);zeros(3),CGa'      ];
                        dRa_dqa       = [dCGa_e1'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1)...
                            dCGa_e2'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1)...
                            dCGa_e3'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1)...
                            dCGa_e4'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1)];
                        dQgyr_r_dRa             = [(dRa_t1*CaB_d*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                            (dRa_t2*CaB_d*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                            (dRa_t3*CaB_d*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                        dQgyr_r_dp{jj}(:,ind4)  =               dQgyr_r_dRa*CGa';
                        dQgyr_r_dq{jj}(:,ind1)  = [(Ra_t*CaG_d*dCBG_e1_blkdiag'*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*CaG_d*dCBG_e2_blkdiag'*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*CaG_d*dCBG_e3_blkdiag'*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*CaG_d*dCBG_e4_blkdiag'*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                        dQgyr_r_dpa             = dQgyr_r_dpa - dQgyr_r_dRa*CGa';
                        dQgyr_r_dqa             = dQgyr_r_dqa + dQgyr_r_dRa*dRa_dqa + ...
                            [(Ra_t*dCaG_e1_blkdiag*blkdiag(CGB_e,CGB_e)*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*dCaG_e2_blkdiag*blkdiag(CGB_e,CGB_e)*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*dCaG_e3_blkdiag*blkdiag(CGB_e,CGB_e)*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                            (Ra_t*dCaG_e4_blkdiag*blkdiag(CGB_e,CGB_e)*xvtilde*Matrices.Aa2{jj}(:,:,ii))*X.x_v_p1{jj}(ind3)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                    end
                end
                %% Aero
                if Sim.aero_flag && Matrices.isLS{jj}
                    %                 v_elem     = Matrices.CBA0{jj}(:,:,ii)*pa;
                    %                 omega_elem = Matrices.CBA0{jj}(:,:,ii)*omega;
                    if Aero.gust_switch == 1 %iter_counter == 0 &&
                        %                     gust_vec{jj}(ind4,:) = [interp3(Aero.gust.x_coord,Aero.gust.y_coord,Aero.gust.z_coord,Aero.gust.x,x_p_p1{jj}(ind4(1)),x_p_p1{jj}(ind4(2)),x_p_p1{jj}(ind4(3)),'linear',0);
                        %                                             interp3(Aero.gust.x_coord,Aero.gust.y_coord,Aero.gust.z_coord,Aero.gust.y,x_p_p1{jj}(ind4(1)),x_p_p1{jj}(ind4(2)),x_p_p1{jj}(ind4(3)),'linear',0);
                        %                                             interp3(Aero.gust.x_coord,Aero.gust.y_coord,Aero.gust.z_coord,Aero.gust.z,x_p_p1{jj}(ind4(1)),x_p_p1{jj}(ind4(2)),x_p_p1{jj}(ind4(3)),'linear',0)];
                        gust_vec{jj}(ind4,:) = getGust(X.x_p_p1{jj}(ind4(1)),X.x_p_p1{jj}(ind4(2)),X.x_p_p1{jj}(ind4(3)),Aero.gust);
                    elseif Aero.gust_switch == 0
                        gust_vec{jj}(ind4,:) = [0;0;0];
                    end
                    if isfield(Aero,'inflow_switch') && Aero.inflow_switch == 1
                        gust_vec{jj}(ind4,:) = gust_vec{jj}(ind4,:) + Aero.inflow_vel;
                    end
                    
                    if 1
                        CA0B_tot_blkdiag = [Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)); zeros(3) eye(3)];
                        CBA0_tot_blkdiag = CA0B_tot_blkdiag';
                    else
                        CA0B_tot_blkdiag =  [Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)+[0;Matrices.chord{jj}(ii)/2;0]); zeros(3) eye(3)];
                        CBA0_tot_blkdiag = ([Matrices.CBA0{jj}(:,:,ii)' zeros(3);zeros(3) Matrices.CBA0{jj}(:,:,ii)']*[eye(3) -skew(Matrices.cp{jj}(:,ii)                               ); zeros(3) eye(3)])';
                    end
                    
                    v_A0     = CA0B_tot_blkdiag*(X.x_v_p1{jj}(ind3) - [CGB_e' zeros(3); zeros(3) zeros(3)]*[gust_vec{jj}(ind4);0;0;0]);
                    %                 v_elem     = v_A0(1:3);
                    %                 omega_elem = v_A0(4:6);
                    %
                    %                 ydot = -v_elem(2);
                    %                 zdot = -v_elem(3);
                    %                 adot =  omega_elem(1);
                    ydot =  v_A0(2);
                    zdot =  v_A0(3);
                    adot =  v_A0(4);
                    
                    %                 rotmat = diag([1 -1 -1 1 1 1]);
                    
                    tiplossFactor = (Aero.lift_dist{jj}(ii)+Aero.lift_dist{jj}(ii+1))/2;
                    
                    chord         = Matrices.chord{jj}(ii);
                    %                 xs            = Matrices.xs{jj}(ii);
                    
                    %                 AoA_qs_t{jj}(ii,1) = 0.5*( - zdot - (chord/2-xs)*adot + 4*ydot^2/chord*[a1*b1 a2*b2]*x_x_p1{jj}(ind5))/ydot;
                    AoA_qs_t{jj}(ii,1) = 0.5*( - zdot + 4*ydot^2/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*X.x_x_p1{jj}(ind5))/ydot;
                    
                    if 0%AoA_qs_t{jj}(ii,1) > Aero.stall_angle/180*pi
                        Lift               = pi*Aero.rho*ydot^2*Aero.stall_angle/180*pi*tiplossFactor*chord;
                        AoA_qs_t{jj}(ii,1) = Aero.stall_angle/180*pi;
                    else
                        %                     Lift           =  + 0.5*pi*Aero.rho*ydot*( - zdot - (chord/2-xs)*adot + 4*ydot^2/chord*[a1*b1 a2*b2]*x_x_p1{jj}(ind5))*tiplossFactor*chord;
                        Lift               =  + 0.5*pi*Aero.rho*ydot*( - zdot + 4*ydot^2/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*X.x_x_p1{jj}(ind5))*tiplossFactor*chord - 0.5*Aero.rho*ydot*ydot*Aero.dCLdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor*chord;
                    end
                    Drag           =  Lift*AoA_qs_t{jj}(ii,1);
                    %                 Mom            = -pi*Aero.rho*chord^3/8/2*ydot*adot*tiplossFactor - Lift*(xs-chord/4);
                    Mom            = -pi*Aero.rho*chord^3/8/2*ydot*adot*tiplossFactor + 0.5*Aero.rho*ydot*ydot*Aero.dCMdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor*chord*chord;
                    aeroForce      = CBA0_tot_blkdiag*[0;Drag;Lift;Mom;0;0];
                    %                 aeroForce_glob = blkdiag(C_twist',C_twist')*blkdiag(CaB_e,CaB_e)*aeroForce;
                    
                    f{jj}(ind3,1)  = aeroForce;
                    %                 f_glob(ind3) = aeroForce_glob;
                    
                    if 1%~Sim.speed_up
                        dF_e      = zeros(6);
                        %                     dF_e(2,2) = +0.5*pi*Aero.rho*(- zdot + 12*ydot^2/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*X.x_x_p1{jj}(ind5))*             tiplossFactor*chord;
                        %                     dF_e(3,3) = -0.5*pi*Aero.rho*ydot                                                                                                               *             tiplossFactor*chord;
                        %                     dF_e(4,2) = - pi*Aero.rho*chord^3/8/2*adot*tiplossFactor;%
                        %                     dF_e(4,4) = - pi*Aero.rho*chord^3/8/2*ydot*tiplossFactor;%
                        dF_e(3,2) = -0.5*pi*Aero.rho*zdot*chord*tiplossFactor + 6*pi*Aero.rho*ydot^2*([Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*X.x_x_p1{jj}(ind5))*tiplossFactor - Aero.rho*ydot*chord*Aero.dCLdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor;
                        dF_e(3,3) = -0.5*pi*Aero.rho*ydot*chord*tiplossFactor;
                        dF_e(2,2) =  AoA_qs_t{jj}(ii,1)*dF_e(3,2) + (0.5*zdot/ydot^2 + 2/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*X.x_x_p1{jj}(ind5))*Lift;
                        dF_e(2,3) =  AoA_qs_t{jj}(ii,1)*dF_e(3,3) -  0.5/ydot                                                                                                           *Lift;
                        dF_e(4,2) = -pi*Aero.rho*chord^3/8/2*adot*tiplossFactor + Aero.rho*ydot*Aero.dCMdd{jj}(:,ii)*Aero.delta_flap{jj}(ii)*tiplossFactor*chord*chord;
                        dF_e(4,4) = -pi*Aero.rho*chord^3/8/2*ydot*tiplossFactor;
                        dF_tot{jj}(ind3,ind3) = CBA0_tot_blkdiag * dF_e * CA0B_tot_blkdiag;
                    else
                        dF_tot{jj}(ind3,ind3) = zeros(6);
                    end
                    
                    if ~Sim.speed_up
                        dF_wg_tot{jj}(ind3,ind4) = - CBA0_tot_blkdiag * dF_e * CA0B_tot_blkdiag * [CGB_e';zeros(3)];
                        dF_u_e                   =   zeros(6,1);
                        dF_u_e(3)                = - 0.5*Aero.rho*ydot*ydot*Aero.dCLdd{jj}(:,ii)*tiplossFactor*chord;
                        dF_u_e(4)                =   0.5*Aero.rho*ydot*ydot*Aero.dCMdd{jj}(:,ii)*tiplossFactor*chord*chord;
                        dF_u_tot{jj}(ind3,ii)    =   CBA0_tot_blkdiag * dF_u_e;
                    else
                        dF_wg_tot{jj}(ind3,ind4) = zeros(6,3);
                        dF_u_tot{jj}(ind3,ii)    = zeros(6,1);
                    end
                    A_x{jj}(ind5,ind5)    = -2*ydot/chord*diag([Aero.Leishman.b1 Aero.Leishman.b2]);
                    if 0%AoA_qs_t{jj}(ii,1) > Aero.stall_angle/180*pi
                        A_x_1{jj}(ind5,1) = Aero.stall_angle/180*pi*[1;1];
                    else
                        %                     A_x_1{jj}(ind5,1) = ( - zdot - (chord/4-xs)*adot)/ydot*[1;1];
                        A_x_1{jj}(ind5,1) = - zdot/ydot*[1;1];
                    end
                    
                    %                 dA_x_1_e              = zeros(2,6);
                    %                 dA_x_1_e(:,2)         = -(- zdot - (chord/4-xs)*adot)/ydot^2*[1;1];
                    %                 dA_x_1_e(:,3)         =                             1/ydot  *[1;1];
                    %                 dA_x_1_e(:,4)         = -(chord/4-xs)                /ydot  *[1;1];
                    %                 dA_x_1{jj}(ind5,ind3) = dA_x_1_e;%*blkdiag(C_twist,C_twist);
                    dA_x_1_e              = zeros(2,6);
                    dA_x_1_e(:,2)         = zdot/ydot^2*[1;1];
                    dA_x_1_e(:,3)         =  - 1/ydot  *[1;1];
                    dA_x_1{jj}(ind5,ind3) = dA_x_1_e*CA0B_tot_blkdiag;
                    
                    if ~Sim.speed_up
                        dA_x_e = zeros(2,6);
                        dA_x_e(:,2) = -2/chord*diag([Aero.Leishman.b1 Aero.Leishman.b2])*X.x_x_p1{jj}(ind5);
                        dA_x{jj}(ind5,ind3)       =   dA_x_e  *CA0B_tot_blkdiag;
                        dA_x_dwg{jj}(ind5,ind4)   = - dA_x_e  *CA0B_tot_blkdiag * [CGB_e';zeros(3)];
                        dA_x_1_dwg{jj}(ind5,ind4) = - dA_x_1_e*CA0B_tot_blkdiag * [CGB_e';zeros(3)];
                        %                 else
                        %                     dA_x_1_dwg{jj}(ind5,ind4) = zeros(2,3);
                    end
                    
                    dF_x_e              =  zeros(6,2);
                    dF_x_e(3,:)         = +0.5*pi*Aero.rho*ydot*(4*ydot^2/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2])*tiplossFactor*chord;
                    dF_x_e(2,:)         =  AoA_qs_t{jj}(ii,1)*dF_x_e(3,:) + 2*ydot/chord*[Aero.Leishman.a1*Aero.Leishman.b1 Aero.Leishman.a2*Aero.Leishman.b2]*Lift;
                    dF_x{jj}(ind3,ind5) =  CBA0_tot_blkdiag * dF_x_e;%blkdiag(C_twist',C_twist')*dF_x_e;
                else
                    f{jj}(ind3,1)  = zeros(6,1);
                end
                %% Gravity
                if Sim.grav_fact
                    gravForce               = Matrices.Aa2{jj}(:,:,ii)*[CGB_e' zeros(3); zeros(3) CGB_e']*Sim.grav_vec;
                    f_grav{jj}(ind3,1)      = gravForce;
                    f_grav_glob{jj}(ind3,1) = [CGB_e' zeros(3); zeros(3) CGB_e']*gravForce;
                    
                    if ~Sim.speed_up
                        dF_grav_tot{jj}(ind3,ind1) = [Matrices.Aa2{jj}(:,:,ii)*dCBG_e1_blkdiag*Sim.grav_vec...
                            Matrices.Aa2{jj}(:,:,ii)*dCBG_e2_blkdiag*Sim.grav_vec...
                            Matrices.Aa2{jj}(:,:,ii)*dCBG_e3_blkdiag*Sim.grav_vec...
                            Matrices.Aa2{jj}(:,:,ii)*dCBG_e4_blkdiag*Sim.grav_vec];
                        if Sim.rb_flag
                            dF_grav_tot2{jj}(:,ind1) = [(Ra_t*CaG_d*dCBG_e1_blkdiag'*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*CaG_d*dCBG_e2_blkdiag'*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*CaG_d*dCBG_e3_blkdiag'*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*CaG_d*dCBG_e4_blkdiag'*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                            dF_grav_dRa              = [(dRa_t1*CaB_d*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                                (dRa_t2*CaB_d*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                                (dRa_t3*CaB_d*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                            dF_grav_dp{jj}(:,ind4)   =               dF_grav_dRa*CGa';
                            dF_grav_dpa              = dF_grav_dpa - dF_grav_dRa*CGa';
                            dF_grav_dqa              = dF_grav_dqa + dF_grav_dRa*dRa_dqa + ...
                                [(Ra_t*dCaG_e1_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*dCaG_e2_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*dCaG_e3_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                                (Ra_t*dCaG_e4_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce)*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                        end
                    else
                        dF_grav_tot{jj}(ind3,ind1) = zeros(6,4);
                    end
                else
                    f_grav{jj}(ind3,1)      = zeros(6,1);
                    if ~Sim.speed_up
                        dF_grav_tot{jj}(ind3,ind1)   = zeros(6,4);
                        if Sim.rb_flag
                            dF_grav_tot2{jj}(:,ind1) = zeros(6,4);
                            dF_grav_dp{jj}(:,ind4)   = zeros(6,3);
                        end
                    end
                end
                %%
                if Sim.aero_flag && Matrices.isLS{jj} && ~Sim.speed_up && Sim.rb_flag
                    dF_tot_dq{jj}(:,ind1) = [(Ra_t*CaG_d*dCBG_e1_blkdiag'*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*CaG_d*dCBG_e2_blkdiag'*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*CaG_d*dCBG_e3_blkdiag'*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*CaG_d*dCBG_e4_blkdiag'*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                    dF_tot_dRa            = [(dRa_t1*CaB_d*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                        (dRa_t2*CaB_d*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)), ...
                        (dRa_t3*CaB_d*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                    dF_tot_dp{jj}(:,ind4) =              dF_tot_dRa*CGa';
                    dF_tot_dpa            = dF_tot_dpa - dF_tot_dRa*CGa';
                    dF_tot_dqa            = dF_tot_dqa + dF_tot_dRa*dRa_dqa + ...
                        [(Ra_t*dCaG_e1_blkdiag*blkdiag(CGB_e,CGB_e)*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*dCaG_e2_blkdiag*blkdiag(CGB_e,CGB_e)*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*dCaG_e3_blkdiag*blkdiag(CGB_e,CGB_e)*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii)),...
                        (Ra_t*dCaG_e4_blkdiag*blkdiag(CGB_e,CGB_e)*f{jj}(ind3,1))*(Matrices.s{jj}(ii+1)-Matrices.s{jj}(ii))];
                else
                    dF_tot_dq{jj}(:,ind1) = zeros(6,4);
                    dF_tot_dp{jj}(:,ind4) = zeros(6,3);
                end
            end
            %%
            
            for ii = 1:nNode
                ind3 = [1:6] + (ii-1)*6;
                if ii == 1
                    Ra    = Matrices.p0{jj};
                    CaB   = eye(3);
                else
                    x_halfnode = [(Matrices.s{jj}(ii)-Matrices.s{jj}(ii-1))/2;0;0];
                    ind1 = [1:4] + (ii-2)*4;
                    ind4 = [1:3] + (ii-2)*3;
                    CGB_e   = Quat2Rot(X.x_q_p1{jj}(ind1));
                    CaB     = CGa'*CGB_e;
                    Ra      = CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*x_halfnode);
                    if ~Sim.speed_up
                        quat = Quat2Rot(X.x_q_p1{jj}(ind1));
                        
                        dCGB_e1 = 2*[ quat(1) -quat(4)  quat(3);
                            quat(4)  quat(1) -quat(2);
                            -quat(3)  quat(2)  quat(1)];
                        
                        dCGB_e2 = 2*[ quat(2)  quat(3)  quat(4);
                            quat(3) -quat(2) -quat(1);
                            quat(4)  quat(1) -quat(2)];
                        
                        dCGB_e3 = 2*[-quat(3)  quat(2)  quat(1);
                            quat(2)  quat(3)  quat(4);
                            -quat(1)  quat(4) -quat(3)];
                        
                        dCGB_e4 = 2*[-quat(4) -quat(1)  quat(2);
                            quat(1) -quat(4)  quat(3);
                            quat(2)  quat(3)  quat(4)];
                        
                        dCBG_e1_blkdiag = [dCGB_e1' zeros(3);zeros(3) dCGB_e1'];
                        dCBG_e2_blkdiag = [dCGB_e2' zeros(3);zeros(3) dCGB_e2'];
                        dCBG_e3_blkdiag = [dCGB_e3' zeros(3);zeros(3) dCGB_e3'];
                        dCBG_e4_blkdiag = [dCGB_e4' zeros(3);zeros(3) dCGB_e4'];
                    end
                end
                CaB_d   = [CaB     ,zeros(3);zeros(3),CaB     ];
                Ra_t    = [eye(3)  ,zeros(3);skew(Ra),eye(3)  ];
                if Sim.rb_flag
                    F_r{jj}(:,ind3)  = Ra_t * CaB_d ;
                    F_rg{jj}(:,ind3) = Ra_t;
                else
                    F_r  = [];
                    F_rg = [];
                end
                if Sim.grav_fact
                    gravForce_pt               = Matrices.M_pt{jj}(:,:,ii)*[CGB_e' zeros(3); zeros(3) CGB_e']*Sim.grav_vec;
                    f_grav_pt{jj}(ind3,1)      = gravForce_pt;
                    f_grav_pt_glob{jj}(ind3,1) = [CGB_e' zeros(3); zeros(3) CGB_e']*gravForce_pt;
                    
                    if ~Sim.speed_up && ii > 1 && Sim.rb_flag
                        dRa1          = [1;0;0];
                        dRa_t1        = [zeros(3),zeros(3);skew(dRa1),zeros(3)];
                        dRa2          = [0;1;0];
                        dRa_t2        = [zeros(3),zeros(3);skew(dRa2),zeros(3)];
                        dRa3          = [0;0;1];
                        dRa_t3        = [zeros(3),zeros(3);skew(dRa3),zeros(3)];
                        CaG_d         = [CGa'    ,zeros(3);zeros(3),CGa'      ];
                        dRa_dqa       = [dCGa_e1'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*x_halfnode)...
                            dCGa_e2'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*x_halfnode)...
                            dCGa_e3'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*x_halfnode)...
                            dCGa_e4'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+CGB_e*x_halfnode)];
                        dRa_dq        = [CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+dCGB_e1*x_halfnode)...
                            CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+dCGB_e2*x_halfnode)...
                            CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+dCGB_e3*x_halfnode)...
                            CGa'*(X.x_p_p1{jj}(ind4)-X.x_pa_p1+dCGB_e4*x_halfnode)];
                        dF_grav_pt_tot{jj}(ind3,ind1) = [Matrices.M_pt{jj}(:,:,ii)*dCBG_e1_blkdiag*Sim.grav_vec...
                            Matrices.M_pt{jj}(:,:,ii)*dCBG_e2_blkdiag*Sim.grav_vec...
                            Matrices.M_pt{jj}(:,:,ii)*dCBG_e3_blkdiag*Sim.grav_vec...
                            Matrices.M_pt{jj}(:,:,ii)*dCBG_e4_blkdiag*Sim.grav_vec];
                        if Sim.rb_flag
                            dF_grav_dRa                 = [(dRa_t1*CaB_d*gravForce_pt), ...
                                (dRa_t2*CaB_d*gravForce_pt), ...
                                (dRa_t3*CaB_d*gravForce_pt)];
                            dF_grav_pt_tot2{jj}(:,ind1) = dF_grav_dRa*dRa_dq + ...
                                [(Ra_t*CaG_d*dCBG_e1_blkdiag'*gravForce_pt),...
                                (Ra_t*CaG_d*dCBG_e2_blkdiag'*gravForce_pt),...
                                (Ra_t*CaG_d*dCBG_e3_blkdiag'*gravForce_pt),...
                                (Ra_t*CaG_d*dCBG_e4_blkdiag'*gravForce_pt)];
                            dF_grav_pt_dp{jj}(:,ind4)   =               dF_grav_dRa*CGa';
                            dF_grav_dpa                 = dF_grav_dpa - dF_grav_dRa*CGa';
                            dF_grav_dqa                 = dF_grav_dqa + dF_grav_dRa*dRa_dqa + ...
                                [(Ra_t*dCaG_e1_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce_pt),...
                                (Ra_t*dCaG_e2_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce_pt),...
                                (Ra_t*dCaG_e3_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce_pt),...
                                (Ra_t*dCaG_e4_blkdiag*blkdiag(CGB_e,CGB_e)*gravForce_pt)];
                        end
                    else
                        dF_grav_pt_tot{jj}(ind3,ind1) = zeros(6,4);
                    end
                else
                    f_grav_pt{jj}(ind3,1)      = zeros(6,1);
                    if ~Sim.speed_up
                        dF_grav_pt_tot{jj}(ind3,ind1)   = zeros(6,4);
                        if Sim.rb_flag
                            dF_grav_pt_tot2{jj}(:,ind1) = zeros(6,4);
                            dF_grav_pt_dp{jj}(:,ind4)   = zeros(6,3);
                        end
                    end
                end
                f_glob{jj}(ind3,1)         = [CGa',zeros(3);zeros(3),CGa']*Matrices.f_a{jj}(ind3);
            end
        end
        %% Output the matrices to the SYS structure
        
        SYS.AoA_qs_t      = AoA_qs_t;
        SYS.Btot_r        = Btot_r  ;
        SYS.Atot          = Atot;
        SYS.dQstiftot     = dQstiftot;
        SYS.dQgyrtot      = dQgyrtot;
        SYS.Ctot          = Ctot;
        SYS.Etot          = Etot;
        SYS.dTgamtot      = dTgamtot;
        
        if Sim.rb_flag
            SYS.M_r           = M_r;
            SYS.Qgyr_r        = Qgyr_r;
            SYS.dQgyr_r       = dQgyr_r;
            SYS.F_r           = F_r;
            SYS.F_rg          = F_rg;
            SYS.dV_dva        = dV_dva;
        end
        
        SYS.OmegaQuat     = OmegaQuat;
        SYS.CGB_tot       = CGB_tot;
        SYS.dOmegaQuat    = dOmegaQuat;
        SYS.dCGB_tot      = dCGB_tot;
        
        if Sim.rb_flag
            SYS.omega_tilde   = omega_tilde;
            SYS.domega_tilde  = domega_tilde;
            SYS.OmegaQuat_a   = OmegaQuat_a;
            SYS.dOmegaQuat_a  = dOmegaQuat_a;
            SYS.CGa_tot       = CGa_tot;
            SYS.dCGa_tot      = dCGa_tot;
        end
        
        SYS.A_x           = A_x;
        SYS.dA_x          = dA_x;
        SYS.dA_x_wg       = dA_x_dwg;
        SYS.A_x_1         = A_x_1;
        SYS.dA_x_1        = dA_x_1;
        SYS.dF_x          = dF_x;
        
        SYS.f              = f;
        SYS.f_glob         = f_glob;
        SYS.f_grav         = f_grav;
        SYS.f_grav_pt      = f_grav_pt;
        SYS.dF_tot         = dF_tot;
        SYS.dF_grav_tot    = dF_grav_tot;
        SYS.dF_grav_pt_tot = dF_grav_pt_tot;
        if Sim.rb_flag
            SYS.f_grav_r      = f_grav_r;
            SYS.dF_grav_tot_r = dF_grav_tot_r;
            SYS.CGa_rot       = [CGa zeros(3); zeros(3) CGa];
        end
        
        SYS.dF_wg     = dF_wg_tot;
        SYS.dF_u      = dF_u_tot;
        SYS.dA_x_1_wg = dA_x_1_dwg;
        
        if ~Sim.speed_up && Sim.rb_flag
            SYS.dQgyr_r_dp      = dQgyr_r_dp;
            SYS.dQgyr_r_dq      = dQgyr_r_dq;
            SYS.dQgyr_r_dpa     = dQgyr_r_dpa;
            SYS.dQgyr_r_dqa     = dQgyr_r_dqa;
            SYS.dF_grav_tot2    = dF_grav_tot2;
            SYS.dF_grav_dp      = dF_grav_dp;
            SYS.dF_grav_pt_tot2 = dF_grav_pt_tot2;
            SYS.dF_grav_pt_dp   = dF_grav_pt_dp;
            SYS.dF_grav_dpa     = dF_grav_dpa;
            SYS.dF_grav_dqa     = dF_grav_dqa;
            SYS.dF_tot_dq       = dF_tot_dq;
            SYS.dF_tot_dp       = dF_tot_dp;
            SYS.dF_tot_dpa      = dF_tot_dpa;
            SYS.dF_tot_dqa      = dF_tot_dqa;
        end
        % elseif Sim.Soln == 2.2
        %     dB_dx1    = kron(Matrices.I_quad ,X.x_p1         );
        % %     dB_dx2    = kron(Matrices.e1_quad,Matrices.I_quad);
        %     SYS.dB_dx = Matrices.B_x_P_quad_p_I*dB_dx1;% + Matrices.B_quad*dB_dx2;
        
end

end