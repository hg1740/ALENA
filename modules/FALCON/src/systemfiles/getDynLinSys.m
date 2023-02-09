function [A,B1,B2,C,D1,D2] = getDynLinSys(X,Matrices,SYS,Sim)

if Sim.Soln == 1
    r_full          = [];
    Q.q_p1_full     = [];
    
    if Sim.tangent_flag
        S_full      = [];
    end
    
    for jj = 1:length(Matrices.n_elem)     
        Qstif{jj} = SYS.Atot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj}) - Matrices.Dtot{jj}*X.x_f_p1{jj};
        
        r_f{jj} = Matrices.R1_bc{jj}*(Qstif{jj} - Matrices.f{jj} - SYS.f_glob{jj} - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}));
        
        r_full        = [r_full;r_f{jj}];
        Q.q_p1_full   = [Q.q_p1_full;X.x_f_p1{jj}];
        
        if Sim.tangent_flag
            S = Matrices.R1_bc{jj}*(SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj});
            S_full    = blkdiag(S_full,S);
        end
    end
    
    % Scale the tangent matrix and make matrices sparse
    if Sim.tangent_flag
        D_scale  = sparse(diag(1./max(abs(S_full),[],2)));
        S_sparse = sparse(S_full);
    else
        D_scale  = [];
        S_sparse = [];
    end
elseif Sim.Soln == 2
    if Sim.method == 1
        A_ = [];
        B_ = [];
        C_ = [];
        D_ = [];
        E_ = [];
        F_ = [];
        G_ = [];
        H_ = [];
        I_ = [];
        J_ = [];
        K_ = [];
        L_ = [];
        M_ = [];
        N_ = zeros(6);
        O_ = [];
        P_ = zeros(6);
        Q_ = [];
        R_ = [];
        S_ = [];
        T_ = [];  
        U_ = [];
        V_ = []; 
        W_ = [];
    else
        Q1_full    = [];
        Q1_fr_full = [];
        Q1_rf_full = [];
        
        Q2_full    = [];
        Q2_fr_full = [];
        Q2_rf_full = [];
        
        Q3_full    = [];
        Q3_rf_full = [];
        
        Q4_full    = [];
        Q4_rf_full = [];
    end

    for jj = 1:length(Matrices.n_elem)        
        n_elem = Matrices.n_elem(jj);
        
        if Sim.method == 1
            A_ = blkdiag(A_, Matrices.R1_bc{jj}*  - Matrices.M{jj});
            B_ = blkdiag(B_, Matrices.R1_bc{jj}*( SYS.Ctot{jj} - SYS.dQgyrtot{jj}  - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} ));
            C_ = blkdiag(C_, Matrices.R1_bc{jj}*( SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj}                   ));
            D_ = blkdiag(D_, Matrices.R1_bc{jj}*(              - Matrices.Btot_mp{jj}*SYS.dF_x{jj}                       ));
            E_ = blkdiag(E_, Matrices.R2_bc{jj}*  - Matrices.T1{jj}                      );
            F_ = blkdiag(F_, Matrices.R2_bc{jj}*( - Matrices.Dtot{jj} - SYS.dTgamtot{jj}));
            G_ = blkdiag(G_, Matrices.R2_bc{jj}*(   SYS.Etot{jj}                        ));
            H_ = [H_;Matrices.R2_bc{jj}*Matrices.vel_in{jj}];
            I_ = blkdiag(I_, - eye(  n_elem*2));
            J_ = blkdiag(J_, - SYS.dA_x_1{jj});
            K_ = blkdiag(K_, - SYS.A_x{jj}   );
            if Sim.rb_flag
                L_ = [L_, - SYS.M_r{jj}];
                M_ = [M_,   SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj} - SYS.Btot_r{jj}*SYS.dF_tot{jj}];
                %         P_ =  P_  - SYS.Btot_r{jj}*SYS.dF_tot{jj}*SYS.dV_dva{jj};
                O_ = [O_,                                    - SYS.Btot_r{jj}*SYS.dF_x{jj}  ];
%                 Q_ = [Q_; - Matrices.R1_bc{jj}*Matrices.Btot_mp{jj}*SYS.dF_tot{jj}*SYS.dV_dva{jj}];
%                 R_ = [R_; - SYS.dA_x_1{jj}*SYS.dV_dva{jj}];
                %         N_ = N_ - SYS.M_r{jj}*SYS.dV_dva{jj};
            end
            S_ = blkdiag(S_, - Matrices.R1_bc{jj}*Matrices.Btot_mp{jj}*SYS.dF_wg{jj});
            T_ = blkdiag(T_, - Matrices.R1_bc{jj}*Matrices.Btot_mp{jj}*SYS.dF_u{jj} );
            
            if Sim.rb_flag
                U_ = [U_, - SYS.Btot_r{jj}*SYS.dF_wg{jj}];
                W_ = [W_, - SYS.Btot_r{jj}*SYS.dF_u{jj}];
            end
            
            V_ = blkdiag(V_, -SYS.dA_x_1_wg{jj});
        else
%             Q1 = [Matrices.R1_bc{jj}*[ Matrices.M{jj}              , zeros( n_elem*6+6, n_elem*6) , zeros(n_elem*6+6,(n_elem  )*2)];
%                   Matrices.R2_bc{jj}*[zeros( n_elem*6+6, n_elem*6) , Matrices.T1{jj}              , zeros(n_elem*6+6,(n_elem  )*2)];
%                   zeros(n_elem*2,n_elem*6)                         , zeros((n_elem  )*2,n_elem*6) , eye(  n_elem*2)               ];
%             
%             Q1_full    = blkdiag(Q1_full,Q1);
%             
%             Q2 = [Matrices.R1_bc{jj}*[ SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} , - Matrices.Btot_mp{jj}*SYS.dF_x{jj}] ;
%                   Matrices.R2_bc{jj}*[ - Matrices.Dtot{jj} - SYS.dTgamtot{jj}                                , SYS.Etot{jj}                                         ,   zeros(n_elem*6+6,n_elem*2)]        ;
%                   -SYS.dA_x_1{jj}                                                                            , zeros(n_elem*2,n_elem*6)                             , - SYS.A_x{jj}                        ];
%             
%             Q2_full    = blkdiag(Q2_full,Q2);

            Qdamp{jj} = (- SYS.dQstiftot{jj} - Matrices.Dtot{jj});
            
            if Sim.aero_flag
                n_aero = 2;
                Q1 = ...
                    [Matrices.R1_bc{jj}*[ Matrices.M{jj}              , Sim.struct_damp*Qdamp{jj}     , zeros(n_elem*6+6,n_elem*4)    , zeros(n_elem*6+6,n_elem*3)    , zeros(n_elem*6+6,n_elem*n_aero)];
                     Matrices.R2_bc{jj}*[ zeros( n_elem*6+6, n_elem*6), Matrices.T1{jj}               , zeros(n_elem*6+6,n_elem*4)    , zeros(n_elem*6+6,n_elem*3)    , zeros(n_elem*6+6,n_elem*n_aero)];
                     zeros(n_elem*4,n_elem*6)                         , zeros(n_elem*4,n_elem*6)      , eye(n_elem*4)                 , zeros(n_elem*4,n_elem*3)      , zeros(n_elem*4,n_elem*n_aero)   ;
                     zeros(n_elem*3,n_elem*6)                         , zeros(n_elem*3,n_elem*6)      , zeros(n_elem*3,n_elem*4)      , eye(n_elem*3)                 , zeros(n_elem*3,n_elem*n_aero)   ;
                     zeros(n_elem*n_aero,n_elem*6)                    , zeros(n_elem*n_aero,n_elem*6) , zeros(n_elem*n_aero,n_elem*4) , zeros(n_elem*n_aero,n_elem*3) , eye(n_elem*n_aero)             ];
                Q1_full    = blkdiag(Q1_full,Q1);
                Q2 = ...
                    [Matrices.R1_bc{jj}*[ SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} , -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}-SYS.dF_grav_pt_tot{jj} , zeros(n_elem*6+6,n_elem*3)    , -Matrices.Btot_mp{jj}*SYS.dF_x{jj}];
                     Matrices.R2_bc{jj}*[-Matrices.Dtot{jj} - SYS.dTgamtot{jj}                                  , SYS.Etot{jj}                                         , zeros(n_elem*6+6,n_elem*4)                                       , zeros(n_elem*6+6,n_elem*3)    , zeros(n_elem*6+6,n_elem*n_aero)   ];
                     0.5*SYS.dOmegaQuat{jj}                                                                     , zeros(n_elem*4,n_elem*6)                             , 0.5*SYS.OmegaQuat{jj}                                            , zeros(n_elem*4,n_elem*3)      , zeros(n_elem*4,n_elem*n_aero)      ;
                    -SYS.CGB_tot{jj}                                                                            , zeros(n_elem*3,n_elem*6)                             , -SYS.dCGB_tot{jj}                                                , zeros(n_elem*3)               , zeros(n_elem*3,n_elem*n_aero)      ;
                    -SYS.dA_x_1{jj} - SYS.dA_x{jj}                                                              , zeros(n_elem*n_aero,n_elem*6)                        , zeros(n_elem*n_aero,n_elem*4)                                    , zeros(n_elem*n_aero,n_elem*3) , - SYS.A_x{jj}                     ];
                Q2_full    = blkdiag(Q2_full,Q2);
                Q3_full    = blkdiag(Q3_full,...
                              [-Matrices.R1_bc{jj}*Matrices.Btot_mp{jj}*SYS.dF_u{jj};
                                zeros(n_elem*6,n_elem);
                                zeros(n_elem*4,n_elem);
                                zeros(n_elem*3,n_elem);
                                zeros(n_elem*2,n_elem)]);
                Q4_full    = blkdiag(Q4_full,...
                              [-Matrices.R1_bc{jj}*Matrices.Btot_mp{jj}*SYS.dF_wg{jj};
                                zeros(n_elem*6,n_elem*3);
                                zeros(n_elem*4,n_elem*3);
                                zeros(n_elem*3,n_elem*3);
                               -SYS.dA_x_1_wg{jj}-SYS.dA_x_wg{jj}]);
                Q3_rf_full = [Q3_rf_full,...
                              [- SYS.Btot_r{jj}*SYS.dF_u{jj}; zeros(4,n_elem);  zeros(3,n_elem)  ]];            
                Q4_rf_full = [Q4_rf_full,...
                              [- SYS.Btot_r{jj}*SYS.dF_wg{jj};zeros(4,n_elem*3);zeros(3,n_elem*3)]];              
            else 
                Q1 = ...
                    [Matrices.R1_bc{jj}*[ Matrices.M{jj}              , Sim.struct_damp*Qdamp{jj}     , zeros(n_elem*6+6,n_elem*4)    , zeros(n_elem*6+6,n_elem*3)];
                     Matrices.R2_bc{jj}*[ zeros( n_elem*6+6, n_elem*6), Matrices.T1{jj}               , zeros(n_elem*6+6,n_elem*4)    , zeros(n_elem*6+6,n_elem*3)];
                     zeros(n_elem*4,n_elem*6)                         , zeros(n_elem*4,n_elem*6)      , eye(n_elem*4)                 , zeros(n_elem*4,n_elem*3)   ;
                     zeros(n_elem*3,n_elem*6)                         , zeros(n_elem*3,n_elem*6)      , zeros(n_elem*3,n_elem*4)      , eye(n_elem*3)             ];
                Q1_full    = blkdiag(Q1_full,Q1);
                Q2 = ...
                    [Matrices.R1_bc{jj}*[ SYS.Ctot{jj} - SYS.dQgyrtot{jj} - Matrices.Btot_mp{jj}*SYS.dF_tot{jj} , SYS.Atot{jj} - Matrices.Dtot{jj} - SYS.dQstiftot{jj} , -Matrices.Btot_mp{jj}*SYS.dF_grav_tot{jj}-SYS.dF_grav_pt_tot{jj} , zeros(n_elem*6+6,n_elem*3)];
                     Matrices.R2_bc{jj}*[-Matrices.Dtot{jj} - SYS.dTgamtot{jj}                                  , SYS.Etot{jj}                                         , zeros(n_elem*6+6,n_elem*4)                                       , zeros(n_elem*6+6,n_elem*3)];
                     0.5*SYS.dOmegaQuat{jj}                                                                     , zeros(n_elem*4,n_elem*6)                             , 0.5*SYS.OmegaQuat{jj}                                            , zeros(n_elem*4,n_elem*3)   ;
                    -SYS.CGB_tot{jj}                                                                            , zeros(n_elem*3,n_elem*6)                             , -SYS.dCGB_tot{jj}                                                , zeros(n_elem*3)           ];
                Q2_full    = blkdiag(Q2_full,Q2);
                Q3_full    = blkdiag(Q3_full,...
                              [-eye( n_elem*6,n_elem*6);
                               zeros(n_elem*6,n_elem*6);
                               zeros(n_elem*4,n_elem*6);
                               zeros(n_elem*3,n_elem*6)]);
                if Sim.rb_flag           
                    Q3_rf_full = [Q3_rf_full,...
                                 [- SYS.F_r{jj}*Matrices.R1_bc{jj}'; zeros(4,n_elem*6);  zeros(3,n_elem*6)  ]];
                end
            end
            
            if Sim.rb_flag
%                 Q1_fr = [zeros(n_elem*6,6);
%                          zeros(n_elem*6,6);
%                          zeros(n_elem*2,6)];
%                 
%                 Q1_rf = [SYS.M_r{jj}, zeros(6,n_elem*6), zeros(6,n_elem*2)];
%                 
%                 Q1_fr_full = [Q1_fr_full;Q1_fr];
%                 Q1_rf_full = [Q1_rf_full,Q1_rf];
%                 
%                 Q2_fr = [zeros(n_elem*6,6),                    ;
%                          Matrices.R2_bc{jj}*Matrices.vel_in{jj};
%                          zeros(n_elem*2,6),                    ];
%                 
%                 Q2_rf = [SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj} - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), - SYS.Btot_r{jj}*SYS.dF_x{jj}];
%                 
%                 Q2_fr_full = [Q2_fr_full;Q2_fr];
%                 Q2_rf_full = [Q2_rf_full,Q2_rf];
                if Sim.aero_flag   
                    Q1_fr = [zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*4,6),      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                             zeros(n_elem*3,6),      zeros(n_elem*3,4),      zeros(n_elem*3,3);
                             zeros(n_elem*n_aero,6), zeros(n_elem*n_aero,4), zeros(n_elem*n_aero,3)];
                    
                    Q1_rf = [SYS.M_r{jj}      , zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3), zeros(6,n_elem*n_aero);
                            zeros(4,n_elem*6) , zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3), zeros(4,n_elem*n_aero);
                            zeros(3,n_elem*6) , zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3), zeros(3,n_elem*n_aero)];
                    
                    Q1_fr_full = [Q1_fr_full;Q1_fr];
                    Q1_rf_full = [Q1_rf_full,Q1_rf];
                    
                    Q2_fr = [zeros(n_elem*6,6),                      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             Matrices.R2_bc{jj}*Matrices.vel_in{jj}, zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*4,6),                      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                             zeros(n_elem*3,6),                      zeros(n_elem*3,4),      zeros(n_elem*3,3);
                             zeros(n_elem*n_aero,6),                 zeros(n_elem*n_aero,4), zeros(n_elem*n_aero,3)];
                    
                    Q2_rf = [SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj}  - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), -SYS.Btot_r{jj}*SYS.dF_grav_tot{jj}-SYS.F_r{jj}*SYS.dF_grav_pt_tot{jj}+SYS.dQgyr_r_dq{jj}-SYS.dF_grav_tot2{jj}-SYS.dF_grav_pt_tot2{jj}-SYS.dF_tot_dq{jj}, SYS.dQgyr_r_dp{jj}-SYS.dF_grav_dp{jj}-SYS.dF_grav_pt_dp{jj}-SYS.dF_tot_dp{jj}, -SYS.Btot_r{jj}*SYS.dF_x{jj};
                             zeros(4,n_elem*6)                                                , zeros(4,n_elem*6), zeros(4,n_elem*4),                                                                                                                                         zeros(4,n_elem*3),                                                              zeros(4,n_elem*n_aero);
                             zeros(3,n_elem*6)                                                , zeros(3,n_elem*6), zeros(3,n_elem*4),                                                                                                                                         zeros(3,n_elem*3),                                                              zeros(3,n_elem*n_aero)];
                    
                    Q2_fr_full = [Q2_fr_full;Q2_fr];
                    Q2_rf_full = [Q2_rf_full,Q2_rf];
                else
                    Q1_fr = [zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*4,6),      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                             zeros(n_elem*3,6),      zeros(n_elem*3,4),      zeros(n_elem*3,3)];
                    
                    Q1_rf = [SYS.M_r{jj}       , zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3);
                             zeros(4,n_elem*6) , zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3);
                             zeros(3,n_elem*6) , zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3)];
                    
                    Q1_fr_full = [Q1_fr_full;Q1_fr];
                    Q1_rf_full = [Q1_rf_full,Q1_rf];
                    
                    Q2_fr = [zeros(n_elem*6,6),                      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             Matrices.R2_bc{jj}*Matrices.vel_in{jj}, zeros(n_elem*6,4),      zeros(n_elem*6,3);
                             zeros(n_elem*4,6),                      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                             zeros(n_elem*3,6),                      zeros(n_elem*3,4),      zeros(n_elem*3,3)];
                    
                    Q2_rf = [SYS.Qgyr_r{jj} + SYS.dQgyr_r{jj}  - SYS.Btot_r{jj}*SYS.dF_tot{jj}, zeros(6,n_elem*6), -SYS.Btot_r{jj}*SYS.dF_grav_tot{jj}-SYS.dF_grav_tot2{jj}, zeros(6,n_elem*3);
                             zeros(4,n_elem*6)                                                , zeros(4,n_elem*6), zeros(4,n_elem*4),                                        zeros(4,n_elem*3);
                             zeros(3,n_elem*6)                                                , zeros(3,n_elem*6), zeros(3,n_elem*4),                                        zeros(3,n_elem*3)];
                    
                    Q2_fr_full = [Q2_fr_full;Q2_fr];
                    Q2_rf_full = [Q2_rf_full,Q2_rf];
                end
            end
        end
    end
    
    if Sim.method == 1
        if Sim.rb_flag
            N_ = N_ - (Matrices.M_rb + 1e-6*eye(6));%Q1_rr_full;

            P_ = P_ + (SYS.omega_tilde*(Matrices.M_rb + 1e-6*eye(6)) + SYS.domega_tilde);%;
        end
    else
        if Sim.rb_flag
%             Q1_rr_full = Matrices.M_rb + 1e-6*eye(6);
%             
%             Q1_full    = [Q1_full   , Q1_fr_full;
%                           Q1_rf_full, Q1_rr_full];
%             
%             Q2_rr_full =  (SYS.omega_tilde*(Matrices.M_rb + 1e-6*eye(6)) + SYS.domega_tilde);%;
%             
%             Q2_full    = [Q2_full   , Q2_fr_full;
%                           Q2_rf_full, Q2_rr_full];
                      
            Q1_rr_full = [Matrices.M_rb + 1e-6*eye(6), zeros(6,4), zeros(6,3);
                          zeros(4,6),                  eye(4),     zeros(4,3);
                          zeros(3,6),                  zeros(3,4), eye(3)   ];
            
            Q1_full    = [Q1_full   , Q1_fr_full;
                          Q1_rf_full, Q1_rr_full];
            
            Q2_rr_full =  [SYS.omega_tilde*(Matrices.M_rb + 1e-6*eye(6)) + SYS.domega_tilde,-SYS.dF_grav_tot_r-SYS.dF_grav_dqa-SYS.dF_tot_dqa,   SYS.dQgyr_r_dpa-SYS.dF_grav_dpa-SYS.dF_tot_dpa;
                           0.5*SYS.dOmegaQuat_a,                                             0.5*SYS.OmegaQuat_a,                                zeros(4,3)                                    ;
                          -SYS.CGa_tot,                                                     -SYS.dCGa_tot,                                       zeros(3)                                     ];
            
            Q2_full    = [Q2_full   , Q2_fr_full;
                          Q2_rf_full, Q2_rr_full];
                    
            Q3_full    = [Q3_full;Q3_rf_full];
            Q4_full    = [Q4_full;Q4_rf_full];
        end
    end
    
%     D_scale  = sparse(diag(1./max(abs(Q1_full),[],2)));        
            
%     A = -Q1_full\Q2_full;
%     
%     test = eig(full(A));
%     scatter(real(test),imag(test)/2/pi,'k.')
%     grid on
%     set(gca,'YLim',[0 150])
%     hold on
%     
%     B = [];
%     C = [];
%     D = [];
    
    %% RB modes
%     n_v = length(A_);
%     
%     invC_ = inv(C_);
%     
%     Q1_rb = [E_*invC_*A_ , -(E_*invC_*B_ + G_*invC_*A_), zeros(n_v,6);
%              zeros(n_v)  ,                     eye(n_v), zeros(n_v,6);
%              zeros(6,n_v),                           L_,           N_];
%          
%     Q2_rb = [zeros(n_v)  ,   F_ - G_*invC_*B_          ,           H_;
%              eye(n_v)    ,    zeros(n_v)               , zeros(n_v,6);
%              zeros(6,n_v),                           M_,           P_];
%          
%     A_rb = Q1_rb\Q2_rb;
%     [vecs_rb,vals_rb] = eig(full(A_rb));
%     vals_rb        = diag(vals_rb);
%     [vals_rb,ind]  = sort(vals_rb);
%     vecs_rb        = vecs_rb(:,ind);
%     scatter(real(vals_rb),imag(vals_rb)/2/pi,'k.')
    
    %% Aeroelastic modes
    if Sim.method == 1
        n_v = length(A_);
        n_x = length(K_);
        
        n_w   = size(S_,2);
        n_u   = size(T_,2);
        
        invC_ = inv(C_);
        
        if ~Sim.rb_flag
            Q1_ae = [E_*invC_*A_   , -E_*invC_*B_ - G_*invC_*A_,   -E_*invC_*D_;
                     zeros(n_v)    ,                   eye(n_v), zeros(n_v,n_x);
                     zeros(n_x,n_v),             zeros(n_x,n_v),             I_];
            
            Q2_ae = [zeros(n_v)    ,           F_ - G_*invC_*B_,   -G_*invC_*D_;
                     eye(n_v)      ,                 zeros(n_v), zeros(n_v,n_x);
                     zeros(n_x,n_v),                         J_,             K_];
            
            Q3_ae = [-G_*invC_*T_, E_*invC_*T_;zeros(n_v,2*n_u);   zeros(n_x,2*n_u)];
            Q4_ae = [-G_*invC_*S_, E_*invC_*S_;zeros(n_v,2*n_w);V_,zeros(n_x,  n_w)];
            
            A_ae  = Q1_ae\Q2_ae;
            
            B1_ae = Q1_ae\Q3_ae;
            B2_ae = Q1_ae\Q4_ae;
            
            C_ae  = [ invC_*A_, -invC_*B_, -invC_*D_];
            
            D1_ae = [-invC_*T_,zeros(n_v,n_u)];
            D2_ae = [-invC_*S_,zeros(n_v,n_w)];
        else
            Q1_ae = [E_*invC_*A_   , -E_*invC_*B_ - G_*invC_*A_,   -E_*invC_*D_, zeros(n_v,6); %-E_*invC_*Q_
                     zeros(n_v)    ,                   eye(n_v), zeros(n_v,n_x), zeros(n_v,6);
                     zeros(n_x,n_v),             zeros(n_x,n_v),             I_, zeros(n_x,6);
                     zeros(6,n_v)  ,                         L_,   zeros(6,n_x),           N_];
            
            Q2_ae = [zeros(n_v)    ,           F_ - G_*invC_*B_,   -G_*invC_*D_,           H_; %-G_*invC_*Q_
                     eye(n_v)      ,                 zeros(n_v), zeros(n_v,n_x), zeros(n_v,6);
                     zeros(n_x,n_v),                         J_,             K_, zeros(n_x,6); %R_
                     zeros(6,n_v)  ,                         M_,             O_,           P_];%zeros(6,6)
            
            Q3_ae = [-G_*invC_*T_, E_*invC_*T_;zeros(n_v,2*n_u);   zeros(n_x,2*n_u);W_,zeros(6,n_u)];
            Q4_ae = [-G_*invC_*S_, E_*invC_*S_;zeros(n_v,2*n_w);V_,zeros(n_x,  n_w);U_,zeros(6,n_w)];
            
            A_ae  = Q1_ae\Q2_ae;
            
            B1_ae = Q1_ae\Q3_ae;
            B2_ae = Q1_ae\Q4_ae;
            
            C_ae  = [ invC_*A_, -invC_*B_, -invC_*D_, zeros(n_v,6)];
            
            D1_ae = [-invC_*T_,zeros(n_v,n_u)];
            D2_ae = [-invC_*S_,zeros(n_v,n_w)];
        end
        
    else
        A_ae  = -Q1_full\Q2_full;
        B1_ae = -Q1_full\Q3_full;
        if Sim.aero_flag
            B2_ae = -Q1_full\Q4_full;
        else
            B2_ae =  [];
        end
            
        C_ae  = eye(size(A_ae));
        
        D1_ae = zeros(size(Q3_full));
        D2_ae = zeros(size(Q4_full));
    end
    A  = A_ae;
    B1 = B1_ae;
    B2 = B2_ae;
    C  = C_ae;
    D1 = D1_ae;
    D2 = D2_ae;
end