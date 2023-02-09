function X_out = getInitialRates(X,X_out,Matrices,SYS,Sim)

if Sim.Soln == 1

elseif Sim.Soln == 2
%     r_full          = [];
%     r_va            = zeros(6,1);
    
    Q.q_p1_full     = [];
    Q.q_dot_p1_full = [];
    
    Q1_full    = [];
    Q1_fr_full = [];
    Q1_rf_full = [];
    
    Q2_full    = [];
    Q2_fr_full = [];
    Q2_rf_full = [];
    
    Q3_full    = [];
    Q3_rr      = zeros(6,1);
        
    for jj = 1:length(Matrices.n_elem)
        mult = 1;

%         Qstif{jj} = SYS.Atot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj}) - Matrices.Dtot{jj}*X.x_f_p1{jj};
%         Qgyr{jj}  = SYS.Ctot{jj}* X.x_v_p1{jj};
%         
%         Qdamp{jj} = (- SYS.dQstiftot{jj} - Matrices.Dtot{jj});
%         
%         Tgam{jj}  = - Matrices.Dtot{jj}*X.x_v_p1{jj} + SYS.Etot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj});%
%         
%         r_f{jj} = Matrices.R1_bc{jj}*(Matrices.M{jj} * X.x_v_dot_p1{jj} + Qgyr{jj} + Qstif{jj} + Sim.struct_damp*Qdamp{jj}* X.x_f_dot_p1{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}));% - M{jj}*x0.x_v{jj}
%         r_v{jj} = Matrices.R2_bc{jj}*(Matrices.T1{jj}* X.x_f_dot_p1{jj} + Tgam{jj} + Matrices.vel_in{jj}*X.x_va_p1);
%         r_q{jj} = X.x_q_dot_p1{jj} + 0.5*SYS.OmegaQuat{jj}*X.x_q_p1{jj};
%         r_p{jj} = X.x_p_dot_p1{jj} - SYS.CGB_tot{jj}*X.x_v_p1{jj};
%         r_x{jj} = X.x_x_dot_p1{jj} - SYS.A_x{jj}*X.x_x_p1{jj} - SYS.A_x_1{jj};
%         
%         r_full          = [r_full;r_f{jj};r_v{jj};r_q{jj};r_p{jj};r_x{jj}];
        
        Q.q_p1_full     = [Q.q_p1_full    ;X.x_v_p1{jj}    ;X.x_f_p1{jj}    ;X.x_q_p1{jj}    ;X.x_p_p1{jj}    ;X.x_x_p1{jj}    ];
        Q.q_dot_p1_full = [Q.q_dot_p1_full;X.x_v_dot_p1{jj};X.x_f_dot_p1{jj};X.x_q_dot_p1{jj};X.x_p_dot_p1{jj};X.x_x_dot_p1{jj}];
        
        
        n_elem = Matrices.n_elem(jj);
        
        if Sim.aero_flag
            Q1 = [R1_bc{jj}*[M{jj}                      , zeros(n_elem*6+6,n_elem*6) , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3)  , zeros(n_elem*6+6,n_elem*2)] ;
                  R2_bc{jj}*[zeros(n_elem*6+6,n_elem*6) , T1{jj}                     , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3)  , zeros(n_elem*6+6,n_elem*2)] ;
                  zeros(n_elem*4,n_elem*6)              , zeros(n_elem*4,n_elem*6)   , eye(n_elem*4)              , zeros(n_elem*4,n_elem*3)    , zeros(n_elem*4,  n_elem*2)  ;
                  zeros(n_elem*3,n_elem*6)              , zeros(n_elem*3,n_elem*6)   , zeros(n_elem*3,n_elem*4)   , eye(n_elem*3)               , zeros(n_elem*3,  n_elem*2)  ;
                  zeros(n_elem*2,n_elem*6)              , zeros(n_elem*2,n_elem*6)   , zeros(n_elem*2,n_elem*4)   , zeros(n_elem*2,n_elem*3)    , eye(  n_elem*2)             ];
            
            Q2 = [R1_bc{jj}*[   Ctot{jj}                , Atot{jj} - Dtot{jj}        , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3)  , zeros(n_elem*6+6,n_elem*2)]  ;
                  R2_bc{jj}*[- Dtot2{jj}                , Etot{jj}                   , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3)  , zeros(n_elem*6+6,n_elem*2)]  ;
                  zeros(n_elem*4,n_elem*6)              , zeros(n_elem*4,n_elem*6)   , 0.5*OmegaQuat{jj}          , zeros(n_elem*4,n_elem*3)    , zeros(n_elem*4  ,n_elem*2)   ;
                  -CGB_tot{jj}                          , zeros(n_elem*3,n_elem*6)   , zeros(n_elem*3,n_elem*4)   , zeros(n_elem*3,n_elem*3)    , zeros(n_elem*3  ,n_elem*2)   ;
                  zeros(n_elem*2,n_elem*6)              , zeros(n_elem*2,n_elem*6)   , zeros(n_elem*2,n_elem*4)   , zeros(n_elem*2,n_elem*3)    ,  - A_x{jj}                  ];
            
            Q3_full    = [Q3_full; + R1_bc{jj}*(- mult*(Matrices.f{jj}+f_glob{jj}) - Btot{jj}*(f{jj}+f_grav{jj}) + Atot{jj}*e1tot{jj});
                          R2_bc{jj}*Etot{jj}*e1tot{jj};
                          zeros(n_elem*4,1);
                          zeros(n_elem*3,1);
                          - A_x_1{jj}];
        else
            Q1 = [Matrices.R1_bc{jj}*[Matrices.M{jj}            , zeros(n_elem*6+6,n_elem*6) , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3) ] ;
                  Matrices.R2_bc{jj}*[zeros(n_elem*6+6,n_elem*6), Matrices.T1{jj}            , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3) ] ;
                  zeros(n_elem*4,n_elem*6)                      , zeros(n_elem*4,n_elem*6)   , eye(n_elem*4)              , zeros(n_elem*4,n_elem*3)     ;
                  zeros(n_elem*3,n_elem*6)                      , zeros(n_elem*3,n_elem*6)   , zeros(n_elem*3,n_elem*4)   , eye(n_elem*3)               ];
            
            Q2 = [Matrices.R1_bc{jj}*[  SYS.Ctot{jj}       , SYS.Atot{jj} - Matrices.Dtot{jj} , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3) ] ;
                  Matrices.R2_bc{jj}*[- Matrices.Dtot{jj}  , SYS.Etot{jj}                     , zeros(n_elem*6+6,n_elem*4) , zeros(n_elem*6+6,n_elem*3) ] ;
                  zeros(n_elem*4,n_elem*6)                 , zeros(n_elem*4,n_elem*6)         , 0.5*SYS.OmegaQuat{jj}      , zeros(n_elem*4,n_elem*3)     ;
                  -SYS.CGB_tot{jj}                         , zeros(n_elem*3,n_elem*6)         , zeros(n_elem*3,n_elem*4)   , zeros(n_elem*3,n_elem*3)    ];
            
            Q3_full    = [Q3_full; + Matrices.R1_bc{jj}*(- mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) + SYS.Atot{jj}*Matrices.e1tot{jj});
                                     Matrices.R2_bc{jj}*SYS.Etot{jj}*Matrices.e1tot{jj};
                                     zeros(n_elem*4,1);
                                     zeros(n_elem*3,1)];
        end
        Q1_full    = blkdiag(Q1_full,Q1);
        Q2_full    = blkdiag(Q2_full,Q2);
        
        if Sim.rb_flag
            Q1_fr = [zeros(n_elem*6,6) , zeros(n_elem*6,4), zeros(n_elem*6,3);
                     zeros(n_elem*6,6) , zeros(n_elem*6,4), zeros(n_elem*6,3);
                     zeros(n_elem*4,6) , zeros(n_elem*4,4), zeros(n_elem*4,3);
                     zeros(n_elem*3,6) , zeros(n_elem*3,4), zeros(n_elem*3,3);
                     zeros(n_elem*2,6) , zeros(n_elem*2,4), zeros(n_elem*2,3)];
            
            Q1_rf = [M_r{jj}          , zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3), zeros(6,n_elem*2);
                     zeros(4,n_elem*6), zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3), zeros(4,n_elem*2);
                     zeros(3,n_elem*6), zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3), zeros(3,n_elem*2)];
            
            Q2_fr = [zeros(n_elem*6,6)    , zeros(n_elem*6,4), zeros(n_elem*6,3);
                     R2_bc{jj}*vel_in{jj} , zeros(n_elem*6,4), zeros(n_elem*6,3);
                     zeros(n_elem*4,6)    , zeros(n_elem*4,4), zeros(n_elem*4,3);
                     zeros(n_elem*3,6)    , zeros(n_elem*3,4), zeros(n_elem*3,3);
                     zeros(n_elem*2,6)    , zeros(n_elem*2,4), zeros(n_elem*2,3)];
            
            Q2_rf = [Qgyr_r{jj}       , zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3), zeros(6,n_elem*2);
                     zeros(4,n_elem*6), zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3), zeros(4,n_elem*2);
                     zeros(3,n_elem*6), zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3), zeros(3,n_elem*2)];
                
            Q1_fr_full = [Q1_fr_full;Q1_fr];
            Q1_rf_full = [Q1_rf_full,Q1_rf];
            Q2_fr_full = [Q2_fr_full;Q2_fr];
            Q2_rf_full = [Q2_rf_full,Q2_rf];
        end
        if Sim.rb_flag
            Q3_rr      = Q3_rr     - mult*(F_r{jj}*Matrices.f{jj} + F_rg{jj}*f_glob{jj}) - Btot_r{jj}*(f{jj}+f_grav{jj}); 
        end
    end
       
    if Sim.rb_flag
        Q.q_p1_full     = [Q.q_p1_full    ;X.x_va_p1    ;X.x_qa_p1    ;X.x_pa_p1    ];
        Q.q_dot_p1_full = [Q.q_dot_p1_full;X.x_va_dot_p1;X.x_qa_dot_p1;X.x_pa_dot_p1];
        
        Q1_rr_full = [Matrices.M_rb + 0e-6*eye(6), zeros(6,4),  zeros(6,3) ;
                      zeros(4,6)                 , eye(4)    ,  zeros(4,3) ;
                      zeros(3,6)                 , zeros(3,4),  eye(3)    ];
        
        Q1_full  = [Q1_full   , Q1_fr_full;
                    Q1_rf_full, Q1_rr_full];
        
        Q2_rr_full = [omega_tilde , zeros(6,4)      , zeros(6,3) ;
                      zeros(4,6)  , 0.5*OmegaQuat_a , zeros(4,3) ;
                      -CGa_tot    , zeros(3,4)      , zeros(3)   ];
        
        Q2_full  = [Q2_full   , Q2_fr_full;
                    Q2_rf_full, Q2_rr_full];
        
        Q3_rr      = Q3_rr  - f_grav_r;
        Q3_full    = [Q3_full;Q3_rr;zeros(4,1);zeros(3,1)];
    end
        
    Q.q_dot_p1_full = - Q1_full\(Q2_full*Q.q_p1_full + Q3_full);

    for jj = 1:length(Matrices.n_elem)       
        X_out.x_f_dot{jj}(:,1) = Q.q_dot_p1_full(Sim.Ind.x_f_ind{jj},1);
        X_out.x_v_dot{jj}(:,1) = Q.q_dot_p1_full(Sim.Ind.x_v_ind{jj},1);
        X_out.x_q_dot{jj}(:,1) = Q.q_dot_p1_full(Sim.Ind.x_q_ind{jj},1);
        X_out.x_p_dot{jj}(:,1) = Q.q_dot_p1_full(Sim.Ind.x_p_ind{jj},1);
        X_out.x_x_dot{jj}(:,1) = Q.q_dot_p1_full(Sim.Ind.x_x_ind{jj},1);
    end
    if Sim.rb_flag
        X_out.x_va_dot(:,1) = Q.q_dot_p1_full(Sim.Ind.x_va_ind,1);
        X_out.x_qa_dot(:,1) = Q.q_dot_p1_full(Sim.Ind.x_qa_ind,1);
        X_out.x_pa_dot(:,1) = Q.q_dot_p1_full(Sim.Ind.x_pa_ind,1);
    end
end