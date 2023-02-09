function [r_full,Q,S_sparse,D_scale] = getResidual(X,Matrices,SYS,Sim, AnalysisParam, bUpdateTangent, tt)

%Default is that the system tangent matrix is not calculated
D_scale  = [];
S_sparse = [];
            
%Number of structural nodes in each part
nNode = Matrices.n_node;
%Number of beam elements in each part
nElem = Matrices.n_elem;
%Number of parts
nPart = numel(nElem);
%Total number of states
nStates = Sim.Ind.x_f_ind{end}(end);

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        r_full          = [];
        Q.q_p1_full     = [];
        
        %Preallocate
        %   - LHS of system
        Qstif = arrayfun(@(nN) zeros(6 * nN, 1), nNode, 'Unif', false);
        %   - Residual
        r_f = arrayfun(@(nE) zeros(6 * nE, 1), nElem, 'Unif', false);  
        %   - System tangent matrix f' = df(q)/dq
        if bUpdateTangent
            S_full = zeros(nStates);
        end
        
        for jj = 1 : nPart
            
            %x_f_p1 is the current value of the states (F & M)
            q = X.x_f_p1{jj};
            
            %Qstif is the LHS of the static equation of motion 
            % [K]{x} = {f} -> res = [K]{x} - {f}
            % {F' ; M'} - [0, F_tilde ;F_tilde, M_tilde] * [C]{F+eF ; M}
            %
            % N.B. 'Dtot' is the spatial derivative vector of the states
            %       Dtot * {F ; M} = {F' ; M'}
            Qstif{jj} = SYS.Atot{jj}*(q + Matrices.e1tot{jj}) - Matrices.Dtot{jj} * q;
            
            %R1_bc denotes that the root of the component is fixed
            %Btot_mp is an integration matrix
            %Anything that is state dependent is generated in 'SYS' i.e. state
            %dependent aero loads
            % * * r_f is the residual * *
            % - Matrices.f is a follower loads
            % - SYS.f_global are point loads from global system (non-follower)
            % - SYS.f are state dependent distributed loads (aero)
            % - SYS.f_grav are state dependent gravity/inertial
            % - SYS.f_grav_pt are state dependent point mass loads (currently not used)
            r_f{jj} = r_f{jj} + Matrices.R1_bc{jj} * ( ...
                Qstif{jj} - Matrices.f{jj} - SYS.f_glob{jj} - ...
                Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});
            
            %Update the tangent matrix
            if bUpdateTangent
                S_full(Sim.Ind.x_f_ind{jj},Sim.Ind.x_f_ind{jj})                         = Matrices.R1_bc{jj}*(SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj});
                if Matrices.Parent{jj}(1)
                    S_full(Sim.Ind.x_f_ind{Matrices.Parent{jj}(1)},Sim.Ind.x_f_ind{jj}) = Matrices.R3_bc{jj}*(SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj});
                end
            end
            
            %Account for connection between members?
            if Matrices.Parent{jj}(1)
                r_f{Matrices.Parent{jj}(1)} = r_f{Matrices.Parent{jj}(1)} + Matrices.R3_bc{jj}*(Qstif{jj} - Matrices.f{jj} - SYS.f_glob{jj} - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});
            end
        end
        
        %Full system state vector (stacked for each member)
        Q.q_p1_full = vertcat(X.x_f_p1{:});
        
        %Full system residual
        r_full = vertcat(r_f{:});
        
        %Scale the tangent matrix and make matrices sparse
        if bUpdateTangent
            %Numerical conditioning
            D_scale  = sparse(diag(1./max(abs(S_full),[],2)));
            S_sparse = sparse(S_full);
        end
        
    otherwise
        
        if Sim.Soln == 2
            r_full          = [];
            r_va            = zeros(6,1);
            
            Q.q_p1_full     = [];
            Q.q_dot_p1_full = [];
            
            r_f = cell(nPart,1);
            r_v = cell(nPart,1);
            for jj = 1:nPart
                r_f{jj} = zeros(6*Matrices.n_elem(jj),1);
                r_v{jj} = zeros(6*Matrices.n_elem(jj),1);
            end
            
            if Sim.tangent_flag
                S_full    = zeros(Sim.Ind.flex_ind{end}(end));
                if Sim.rb_flag
                    S_fr_full = zeros(Sim.Ind.flex_ind{end}(end),6+4+3);
                    S_rf_full = zeros(6+4+3,Sim.Ind.flex_ind{end}(end));
                end
            end
            for jj = 1:nPart
                if 1%Sim.time_vec(tt)<9999
                    mult = 1;
                    %             mult = (200+1e-1*(1-cos(4*pi*Sim.time_vec(tt+1)))/2)/200;
                    %             mult = (1-cos(4*pi*Sim.time_vec(tt+1)))/2;
                else
                    mult = 0;
                end
                
                Qstif{jj} = SYS.Atot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj}) - Matrices.Dtot{jj}*X.x_f_p1{jj};
                Qgyr{jj}  = SYS.Ctot{jj}* X.x_v_p1{jj};
                
                %         Qdamp{jj} = (SYS.Atot{jj} - Matrices.Dtot{jj})*X.x_f_dot_p1{jj};
                Qdamp{jj} = (- SYS.dQstiftot{jj} - Matrices.Dtot{jj});
                
                Tgam{jj}  = - Matrices.Dtot{jj}*X.x_v_p1{jj} + SYS.Etot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj});%
                
                r_f{jj} = r_f{jj} + Matrices.R1_bc{jj}*(Matrices.M{jj} * X.x_v_dot_p1{jj} + Qgyr{jj} + Qstif{jj} + Sim.struct_damp*Qdamp{jj}* X.x_f_dot_p1{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});% - M{jj}*x0.x_v{jj}
                r_v{jj} = r_v{jj} + Matrices.R2_bc{jj}*(Matrices.T1{jj}* X.x_f_dot_p1{jj} + Tgam{jj});
                r_q{jj} =           X.x_q_dot_p1{jj} + 0.5*SYS.OmegaQuat{jj}*X.x_q_p1{jj};
                r_p{jj} =           X.x_p_dot_p1{jj} - SYS.CGB_tot{jj}*X.x_v_p1{jj};
                r_x{jj} =           X.x_x_dot_p1{jj} - SYS.A_x{jj}*X.x_x_p1{jj} - SYS.A_x_1{jj};
                
                if Matrices.Parent{jj}(1)
                    parent_ind = Matrices.Parent{jj}(1);
                    r_f{parent_ind} = r_f{parent_ind} + Matrices.R3_bc{jj}*(Matrices.M{jj} * X.x_v_dot_p1{jj} + Qgyr{jj} + Qstif{jj} + Sim.struct_damp*Qdamp{jj}* X.x_f_dot_p1{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});
                    r_v{jj}         = r_v{jj}         + Matrices.R4_bc{jj}*(Matrices.T1{parent_ind}* X.x_f_dot_p1{parent_ind} - Matrices.Dtot{parent_ind}*X.x_v_p1{parent_ind} + SYS.Etot{parent_ind}*(X.x_f_p1{parent_ind}+Matrices.e1tot{parent_ind}));
                else
                    r_v{jj} = r_v{jj} + Matrices.R2_bc{jj}*(Matrices.vel_in{jj}*X.x_va_p1);
                end
                
                %         r_full          = [r_full;r_f{jj};r_v{jj};r_q{jj};r_p{jj};r_x{jj}];
                
                Q.q_p1_full     = [Q.q_p1_full    ;X.x_v_p1{jj}    ;X.x_f_p1{jj}    ;X.x_q_p1{jj}    ;X.x_p_p1{jj}    ;X.x_x_p1{jj}    ];
                Q.q_dot_p1_full = [Q.q_dot_p1_full;X.x_v_dot_p1{jj};X.x_f_dot_p1{jj};X.x_q_dot_p1{jj};X.x_p_dot_p1{jj};X.x_x_dot_p1{jj}];
                
                if Sim.tangent_flag
                    n_elem = Matrices.n_elem(jj);
                    if Sim.aero_flag
                        n_aero = 2;
                        S_full(Sim.Ind.flex_ind{jj},Sim.Ind.flex_ind{jj}) = ...
                            [Matrices.R1_bc{jj}*[ Matrices.M{jj}/Sim.gam1/Sim.time.dt + SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} + Sim.struct_damp*Qdamp{jj}/Sim.gam1/Sim.time.dt, -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}                  , zeros(n_elem*6+6,n_elem*3)         , -Matrices.Btot_mp{jj}*SYS.dF_x{jj}];
                            Matrices.R2_bc{jj}*[-Matrices.Dtot{jj}      - SYS.dTgamtot{jj}                                                                   , Matrices.T1{jj}/Sim.gam1/Sim.time.dt  + SYS.Etot{jj}                                                 , zeros(n_elem*6+6,n_elem*4)                                 , zeros(n_elem*6+6,n_elem*3)         , zeros(n_elem*6+6,n_elem*n_aero)]  ;
                            0.5*SYS.dOmegaQuat{jj}                                                                                                           , zeros(n_elem*4,n_elem*6)                                                                             , eye(n_elem*4)/Sim.gam1/Sim.time.dt + 0.5*SYS.OmegaQuat{jj} , zeros(n_elem*4,n_elem*3)           , zeros(n_elem*4,n_elem*n_aero)     ;
                            -SYS.CGB_tot{jj}                                                                                                                  , zeros(n_elem*3,n_elem*6)                                                                             , -SYS.dCGB_tot{jj}                                          , eye(n_elem*3)/Sim.gam1/Sim.time.dt , zeros(n_elem*3,n_elem*n_aero)     ;
                            -SYS.dA_x_1{jj}                                                                                                                   , zeros(n_elem*n_aero,n_elem*6)                                                                        , zeros(n_elem*n_aero,(n_elem  )*4)                          , zeros(n_elem*n_aero,n_elem*3)      , eye(  n_elem*n_aero)/Sim.gam1/Sim.time.dt - SYS.A_x{jj}];
                        if Matrices.Parent{jj}(1)
                            n_elem_p = Matrices.n_elem(parent_ind);
                            S_full(Sim.Ind.x_v_ind{parent_ind},Sim.Ind.flex_ind{jj}) = Matrices.R3_bc{jj}*[ Matrices.M{jj}/Sim.gam1/Sim.time.dt + SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} + Sim.struct_damp*Qdamp{jj}/Sim.gam1/Sim.time.dt, -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}                  , zeros(n_elem*6+6,n_elem*3)         , Matrices.Btot_mp{jj}*SYS.dF_x{jj}];
                            S_full(Sim.Ind.x_f_ind{jj},Sim.Ind.flex_ind{parent_ind}) = Matrices.R4_bc{jj}*[-Matrices.Dtot{parent_ind}      - SYS.dTgamtot{parent_ind}                                                   , Matrices.T1{parent_ind}/Sim.gam1/Sim.time.dt  + SYS.Etot{parent_ind}                                 , zeros(n_elem_p*6+6,n_elem_p*4)                             , zeros(n_elem_p*6+6,n_elem_p*3)     , zeros(n_elem_p*6+6,n_elem_p*n_aero)];
                        end
                    else
                        S_full(Sim.Ind.flex_ind{jj},Sim.Ind.flex_ind{jj}) = ...
                            [Matrices.R1_bc{jj}*[ Matrices.M{jj}/Sim.gam1/Sim.time.dt + SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} + Sim.struct_damp*Qdamp{jj}/Sim.gam1/Sim.time.dt, -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}                  , zeros(n_elem*6+6,n_elem*3)         ];
                            Matrices.R2_bc{jj}*[-Matrices.Dtot{jj}      - SYS.dTgamtot{jj}                                                                   , Matrices.T1{jj}/Sim.gam1/Sim.time.dt  + SYS.Etot{jj}                                                 , zeros(n_elem*6+6,n_elem*4)                                 , zeros(n_elem*6+6,n_elem*3)         ];
                            0.5*SYS.dOmegaQuat{jj}                                                                                                           , zeros(n_elem*4,n_elem*6)                                                                             , eye(n_elem*4)/Sim.gam1/Sim.time.dt + 0.5*SYS.OmegaQuat{jj} , zeros(n_elem*4,n_elem*3)           ;
                            -SYS.CGB_tot{jj}                                                                                                                  , zeros(n_elem*3,n_elem*6)                                                                             , -SYS.dCGB_tot{jj}                                          , eye(n_elem*3)/Sim.gam1/Sim.time.dt ];
                        if Matrices.Parent{jj}(1)
                            n_elem_p = Matrices.n_elem(parent_ind);
                            S_full(Sim.Ind.x_v_ind{parent_ind},Sim.Ind.flex_ind{jj}) = Matrices.R3_bc{jj}*[ Matrices.M{jj}/Sim.gam1/Sim.time.dt + SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} + Sim.struct_damp*Qdamp{jj}/Sim.gam1/Sim.time.dt, -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}                  , zeros(n_elem*6+6,n_elem*3)];
                            S_full(Sim.Ind.x_f_ind{jj},Sim.Ind.flex_ind{parent_ind}) = Matrices.R4_bc{jj}*[-Matrices.Dtot{parent_ind}      - SYS.dTgamtot{parent_ind}                                                   , Matrices.T1{parent_ind}/Sim.gam1/Sim.time.dt  + SYS.Etot{parent_ind}                                 , zeros(n_elem_p*6+6,n_elem_p*4)                             , zeros(n_elem_p*6+6,n_elem_p*3)];
                        end
                    end
                    
                    if Sim.rb_flag
                        if Sim.aero_flag
                            if ~Matrices.Parent{jj}(1)
                                S_fr_full(Sim.Ind.flex_ind{jj},:) = ...
                                    [zeros(n_elem*6,6),                      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                                    Matrices.R2_bc{jj}*Matrices.vel_in{jj}, zeros(n_elem*6,4),      zeros(n_elem*6,3);
                                    zeros(n_elem*4,6),                      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                                    zeros(n_elem*3,6),                      zeros(n_elem*3,4),      zeros(n_elem*3,3);
                                    zeros(n_elem*n_aero,6),                 zeros(n_elem*n_aero,4), zeros(n_elem*n_aero,3)];
                            end
                            S_rf_full(:,Sim.Ind.flex_ind{jj}) = ...
                                [SYS.M_r{jj}/Sim.gam1/Sim.time.dt + SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj}  - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3), zeros(6,n_elem*n_aero);
                                zeros(4,n_elem*6)                                                                                   , zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3), zeros(4,n_elem*n_aero);
                                zeros(3,n_elem*6)                                                                                   , zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3), zeros(3,n_elem*n_aero)];
                        else
                            if ~Matrices.Parent{jj}(1)
                                S_fr_full(Sim.Ind.flex_ind{jj},:) = ...
                                    [zeros(n_elem*6,6),                      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                                    Matrices.R2_bc{jj}*Matrices.vel_in{jj}, zeros(n_elem*6,4),      zeros(n_elem*6,3);
                                    zeros(n_elem*4,6),                      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                                    zeros(n_elem*3,6),                      zeros(n_elem*3,4),      zeros(n_elem*3,3)];
                            end
                            S_rf_full(:,Sim.Ind.flex_ind{jj}) = ...
                                [SYS.M_r{jj}/Sim.gam1/Sim.time.dt + SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj}  - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3);
                                zeros(4,n_elem*6)                                                                                   , zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3);
                                zeros(3,n_elem*6)                                                                                   , zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3)];
                        end
                        %                 S_fr_full = [S_fr_full;S_fr];
                        %                 S_rf_full = [S_rf_full,S_rf];
                    end
                end
                if Sim.rb_flag
                    r_va =   r_va + SYS.M_r{jj}*X.x_v_dot_p1{jj} + SYS.Qgyr_r{jj}*X.x_v_p1{jj} - mult*(SYS.F_r{jj}*Matrices.f{jj} + SYS.F_rg{jj}*SYS.f_glob{jj}) - SYS.Btot_r{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.F_r{jj}*SYS.f_grav_pt{jj};
                end
            end
            for jj = 1:nPart
                r_full        = [r_full;r_f{jj};r_v{jj};r_q{jj};r_p{jj};r_x{jj}];
            end
            
            if Sim.rb_flag
                r_va =   r_va + Matrices.M_rb*X.x_va_dot_p1 + SYS.omega_tilde*Matrices.M_rb*X.x_va_p1 - SYS.f_grav_r;
                r_qa =                        X.x_qa_dot_p1 + 0.5*SYS.OmegaQuat_a*X.x_qa_p1;
                r_pa =                        X.x_pa_dot_p1 -     SYS.CGa_tot    *X.x_va_p1;
                
                r_full = [r_full;r_va;r_qa;r_pa];
                
                Q.q_p1_full     = [Q.q_p1_full    ;X.x_va_p1    ;X.x_qa_p1    ;X.x_pa_p1    ];
                Q.q_dot_p1_full = [Q.q_dot_p1_full;X.x_va_dot_p1;X.x_qa_dot_p1;X.x_pa_dot_p1];
            end
            
            if Sim.tangent_flag
                if Sim.rb_flag
                    S_rr_full = [Matrices.M_rb/Sim.gam1/Sim.time.dt + SYS.omega_tilde*Matrices.M_rb + SYS.domega_tilde,   SYS.dF_grav_tot_r,                                 zeros(6,3)     ;
                        SYS.dOmegaQuat_a,                                                                        eye(4)/Sim.gam1/Sim.time.dt + 0.5*SYS.OmegaQuat_a, zeros(4,3)     ;
                        -SYS.CGa_tot,                                                                            -SYS.dCGa_tot,                                      eye(3)/Sim.gam1/Sim.time.dt];
                    
                    S_full  = [S_full   , S_fr_full;
                        S_rf_full, S_rr_full];
                end
                
                D_scale  = sparse(diag(1./max(abs(S_full),[],2)));
                S_sparse = sparse(S_full);
            else
                D_scale  = [];
                S_sparse = [];
            end
            
        elseif Sim.Soln == 2.1
            r_full{1}          = [];
            r_full{2}          = [];
            r_va               = zeros(6,1);
            
            Q.q_p1_full{1}     = [];
            Q.q_dot_p1_full{1} = [];
            Q.q_p1_full{2}     = [];
            Q.q_dot_p1_full{2} = [];
            
            if Sim.tangent_flag
                S_full{1}    = [];
                S_fr_full{1} = [];
                S_rf_full{1} = [];
                S_full{2}    = [];
                S_fr_full{2} = [];
                S_rf_full{2} = [];
            end
            for jj = 1:nPart
                if Sim.time_vec(tt)<9999
                    mult = 1;
                    %                 mult = sin(2.26*2*pi*time(tt));
                    %                 mult = 0.5*(sin(pi*time(tt)-pi/2)+1);
                    %                 mult = cos(2*pi*time(tt)+270*pi/180);
                else
                    mult = 0;
                end
                
                Qstif{jj} = SYS.Atot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj}) - Matrices.Dtot{jj}*X.x_f_p1{jj};
                Qgyr{jj}  = SYS.Ctot{jj}* X.x_v_p1{jj};
                
                Qdamp{jj} = (SYS.Atot{jj} - Matrices.Dtot{jj})*X.x_f_dot_p1{jj};
                
                Tgam{jj}  = - Matrices.Dtot{jj}*X.x_v_p1{jj} + SYS.Etot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj});%
                
                r_f{jj} = Matrices.R1_bc{jj}*(Matrices.M{jj} * X.x_v_dot_p1{jj} + Qgyr{jj} + Qstif{jj} + Sim.struct_damp*Qdamp{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}));% - M{jj}*x0.x_v{jj}
                r_v{jj} = Matrices.R2_bc{jj}*(Matrices.T1{jj}* X.x_f_dot_p1{jj} + Tgam{jj} + Matrices.vel_in{jj}*X.x_va_p1);
                r_q{jj} = X.x_q_dot_p1{jj} + 0.5*SYS.OmegaQuat{jj}*X.x_q_p1{jj};
                r_p{jj} = X.x_p_dot_p1{jj} - SYS.CGB_tot{jj}*X.x_v_p1{jj};
                r_x{jj} = X.x_x_dot_p1{jj} - SYS.A_x{jj}*X.x_x_p1{jj} - SYS.A_x_1{jj};
                
                r_full{1}       = [r_full{1};r_f{jj};r_v{jj};r_x{jj}];
                r_full{2}       = [r_full{2};r_q{jj};r_p{jj}        ];
                
                Q.q_p1_full{1}     = [Q.q_p1_full{1}    ;X.x_v_p1{jj}    ;X.x_f_p1{jj}    ; X.x_x_p1{jj}    ];
                Q.q_dot_p1_full{1} = [Q.q_dot_p1_full{1};X.x_v_dot_p1{jj};X.x_f_dot_p1{jj}; X.x_x_dot_p1{jj}];
                Q.q_p1_full{2}     = [Q.q_p1_full{2}    ;X.x_q_p1{jj}    ;X.x_p_p1{jj}                      ];
                Q.q_dot_p1_full{2} = [Q.q_dot_p1_full{2};X.x_q_dot_p1{jj};X.x_p_dot_p1{jj}                  ];
                
                if Sim.tangent_flag
                    n_elem = Matrices.n_elem(jj);
                    if Sim.aero_flag
                        n_aero = 2;
                    else
                        n_aero = 0;
                    end
                    
                    S1 = [Matrices.R1_bc{jj}*[ Matrices.M{jj}/Sim.gam1/Sim.time.dt + SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} , Matrices.Btot_mp{jj}*SYS.dF_x{jj}];
                        Matrices.R2_bc{jj}*[-Matrices.Dtot{jj}      - SYS.dTgamtot{jj}                                                                   , Matrices.T1{jj}/Sim.gam1/Sim.time.dt  + SYS.Etot{jj} , zeros(n_elem*6+6,n_elem*n_aero)]  ;
                        -SYS.dA_x_1{jj}                                                                                                                  , zeros(n_elem*n_aero,n_elem*6)                        , eye(  n_elem*n_aero)/Sim.gam1/Sim.time.dt - SYS.A_x{jj}];
                    
                    S_fr1 = [zeros(n_elem*6,6),                    ;
                        Matrices.R2_bc{jj}*Matrices.vel_in{jj};
                        zeros(n_elem*n_aero,6),               ];
                    
                    S_rf1 = [SYS.M_r{jj}/Sim.gam1/Sim.time.dt + SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj}  - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), zeros(6,n_elem*n_aero)];
                    
                    S_full{1}    = blkdiag(S_full{1},S1);
                    S_fr_full{1} = [S_fr_full{1};S_fr1];
                    S_rf_full{1} = [S_rf_full{1},S_rf1];
                    
                    S2 = [eye(n_elem*4)/Sim.gam1/Sim.time.dt + 0.5*SYS.OmegaQuat{jj} , zeros(n_elem*4,n_elem*3);
                        -SYS.dCGB_tot{jj}                                          , eye(n_elem*3)/Sim.gam1/Sim.time.dt];
                    
                    S_fr2 = [zeros(n_elem*4,4),      zeros(n_elem*4,3);
                        zeros(n_elem*3,4),      zeros(n_elem*3,3)];
                    
                    S_rf2 = [zeros(4,n_elem*4), zeros(4,n_elem*3);
                        zeros(3,n_elem*4), zeros(3,n_elem*3)];
                    
                    S_full{2}    = blkdiag(S_full{2},S2);
                    S_fr_full{2} = [S_fr_full{2};S_fr2];
                    S_rf_full{2} = [S_rf_full{2},S_rf2];
                end
                r_va =   r_va + SYS.M_r{jj}*X.x_v_dot_p1{jj} + SYS.Qgyr_r{jj}*X.x_v_p1{jj} - mult*(SYS.F_r{jj}*Matrices.f{jj} + SYS.F_rg{jj}*SYS.f_glob{jj}) - SYS.Btot_r{jj}*(SYS.f{jj}+SYS.f_grav{jj});
            end
            
            r_va =   r_va + Matrices.M_rb*X.x_va_dot_p1 + SYS.omega_tilde*Matrices.M_rb*X.x_va_p1 - SYS.f_grav_r;
            r_qa =                        X.x_qa_dot_p1 + 0.5*SYS.OmegaQuat_a*X.x_qa_p1;
            r_pa =                        X.x_pa_dot_p1 -     SYS.CGa_tot    *X.x_va_p1;
            
            r_full{1} = [r_full{1};r_va];
            r_full{2} = [r_full{2};r_qa;r_pa];
            
            if Sim.tangent_flag
                S_rr_full1 = [Matrices.M_rb/Sim.gam1/Sim.time.dt + SYS.omega_tilde + SYS.domega_tilde];
                
                S_full{1}  = [S_full{1}   , S_fr_full{1};
                    S_rf_full{1}, S_rr_full1];
                
                D_scale{1}  = sparse(diag(1./max(abs(S_full{1}),[],2)));
                S_sparse{1} = sparse(S_full{1});
                
                S_rr_full2 = [eye(4)/Sim.gam1/Sim.time.dt + 0.5*SYS.OmegaQuat_a, zeros(4,3)                  ;
                    -SYS.dCGa_tot,                                      eye(3)/Sim.gam1/Sim.time.dt];
                
                S_full{2} = [S_full{2}   , S_fr_full{2};
                    S_rf_full{2}, S_rr_full2];
                
                D_scale{2}  = sparse(diag(1./max(abs(S_full{2}),[],2)));
                S_sparse{2} = sparse(S_full{2});
            else
                D_scale{1}  = [];
                S_sparse{1} = [];
                D_scale{2}  = [];
                S_sparse{2} = [];
            end
            
            Q.q_p1_full{1}     = [Q.q_p1_full{1}    ;X.x_va_p1];
            Q.q_dot_p1_full{1} = [Q.q_dot_p1_full{1};X.x_va_dot_p1];
            Q.q_p1_full{2}     = [Q.q_p1_full{2}    ;X.x_qa_p1    ;X.x_pa_p1    ];
            Q.q_dot_p1_full{2} = [Q.q_dot_p1_full{2};X.x_qa_dot_p1;X.x_pa_dot_p1];
        elseif Sim.Soln == 2.2
            xkronx     = kron(X.x_p1,X.x_p1);
            r_full     = Matrices.C_quad*X.x_dot_p1 + Matrices.A_quad*X.x_p1 + Matrices.B_quad*xkronx;
            Q.q_p1     = X.x_p1;
            Q.q_dot_p1 = X.x_dot_p1;
            S_sparse   = Matrices.C_quad/Sim.gam1/Sim.time.dt + Matrices.A_quad + SYS.dB_dx;
            D_scale    = [];%eye(size(S_sparse));%sparse(diag(1./max(abs(S_sparse),[],2)));
        end
        
end

end