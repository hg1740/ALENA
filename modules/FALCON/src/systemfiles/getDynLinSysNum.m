function [A,B1,B2,C,D1,D2] = getDynLinSysNum(X_p1,Matrices,Aero,Sim)

if Sim.Soln == 1
    error('Not supported yet!')
elseif Sim.Soln == 2
    Sim.eps = 1e-6;
    tic
    SYS  = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
    vec0 = getSysVec(X_p1,Matrices,SYS,Sim);
    toc
    X_p1_store = X_p1;
    A = zeros(length(vec0));
    for jj = 1:length(Matrices.n_elem)
        for kk = 1:length(Sim.Ind.x_v_ind{jj})
            tic
            X_p1.x_v_p1{jj}(kk)          = X_p1.x_v_p1{jj}(kk) + Sim.eps;
            SYS                          = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                         = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_v_ind{jj}(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                         = X_p1_store;
            toc
        end
        for kk = 1:length(Sim.Ind.x_f_ind{jj})
            tic
            X_p1.x_f_p1{jj}(kk)          = X_p1.x_f_p1{jj}(kk) + Sim.eps;
            SYS                          = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                         = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_f_ind{jj}(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                         = X_p1_store;
            toc
        end
        for kk = 1:length(Sim.Ind.x_q_ind{jj})
            tic
            X_p1.x_q_p1{jj}(kk)          = X_p1.x_q_p1{jj}(kk) + Sim.eps;
            SYS                          = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                         = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_q_ind{jj}(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                         = X_p1_store;
            toc
        end
        for kk = 1:length(Sim.Ind.x_p_ind{jj})
            tic
            X_p1.x_p_p1{jj}(kk)          = X_p1.x_p_p1{jj}(kk) + Sim.eps;
            SYS                          = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                         = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_p_ind{jj}(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                         = X_p1_store;
            toc
        end
        if Sim.aero_flag
            for kk = 1:length(Sim.Ind.x_x_ind{jj})
                tic
                X_p1.x_x_p1{jj}(kk)          = X_p1.x_x_p1{jj}(kk) + Sim.eps;
                SYS                          = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
                vec1                         = getSysVec(X_p1,Matrices,SYS,Sim);
                A(:,Sim.Ind.x_x_ind{jj}(kk)) = (vec1-vec0)/Sim.eps;
                X_p1                         = X_p1_store;
                toc
            end
        end
    end
    if Sim.rb_flag
        for kk = 1:length(Sim.Ind.x_va_ind)
            tic
            X_p1.x_va_p1(kk)          = X_p1.x_va_p1(kk) + Sim.eps;
            SYS                       = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                      = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_va_ind(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                      = X_p1_store;
            toc
        end
        for kk = 1:length(Sim.Ind.x_qa_ind)
            tic
            X_p1.x_qa_p1(kk)          = X_p1.x_qa_p1(kk) + Sim.eps;
            SYS                       = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                      = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_qa_ind(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                      = X_p1_store;
            toc
        end
        for kk = 1:length(Sim.Ind.x_pa_ind)
            tic
            X_p1.x_pa_p1(kk)          = X_p1.x_pa_p1(kk) + Sim.eps;
            SYS                       = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
            vec1                      = getSysVec(X_p1,Matrices,SYS,Sim);
            A(:,Sim.Ind.x_pa_ind(kk)) = (vec1-vec0)/Sim.eps;
            X_p1                      = X_p1_store;
            toc
        end
    end
    [B1,B2,C,D1,D2] = deal([]);
end

function vec_out = getSysVec(X,Matrices,SYS,Sim)
r_full          = [];
r_va            = zeros(6,1);

mult = 1;

r_f = cell(length(Matrices.n_elem),1);
r_v = cell(length(Matrices.n_elem),1);
for jj = 1:length(Matrices.n_elem)
    r_f{jj} = zeros(6*Matrices.n_elem(jj),1);
    r_v{jj} = zeros(6*Matrices.n_elem(jj),1);
end

S_full    = zeros(Sim.Ind.flex_ind{end}(end));
if Sim.rb_flag
    S_fr_full = zeros(Sim.Ind.flex_ind{end}(end),6+4+3);
    S_rf_full = zeros(6+4+3,Sim.Ind.flex_ind{end}(end));
end

for jj = 1:length(Matrices.n_elem)
    Qstif{jj} = SYS.Atot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj}) - Matrices.Dtot{jj}*X.x_f_p1{jj};
    Qgyr{jj}  = SYS.Ctot{jj}* X.x_v_p1{jj};
    
    Qdamp{jj} = (- SYS.dQstiftot{jj} - Matrices.Dtot{jj});
    
    Tgam{jj}  = - Matrices.Dtot{jj}*X.x_v_p1{jj} + SYS.Etot{jj}*(X.x_f_p1{jj}+Matrices.e1tot{jj});%
    
    r_f{jj} = r_f{jj} + Matrices.R1_bc{jj}*(Qgyr{jj} + Qstif{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});% - M{jj}*x0.x_v{jj}
    r_v{jj} = r_v{jj} + Matrices.R2_bc{jj}*(Tgam{jj});
    r_q{jj} =           0.5*SYS.OmegaQuat{jj}*X.x_q_p1{jj};
    r_p{jj} =           - SYS.CGB_tot{jj}*X.x_v_p1{jj};
    r_x{jj} =           - SYS.A_x{jj}*X.x_x_p1{jj} - SYS.A_x_1{jj};
    
    if Matrices.Parent{jj}(1)
        parent_ind = Matrices.Parent{jj}(1);
        r_f{parent_ind} = r_f{parent_ind} + Matrices.R3_bc{jj}*(Qgyr{jj} + Qstif{jj} - mult*(Matrices.f{jj}+SYS.f_glob{jj}) - Matrices.Btot_mp{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.f_grav_pt{jj});
        r_v{jj}         = r_v{jj}         + Matrices.R4_bc{jj}*(- Matrices.Dtot{parent_ind}*X.x_v_p1{parent_ind} + SYS.Etot{parent_ind}*(X.x_f_p1{parent_ind}+Matrices.e1tot{parent_ind}));
    else
        r_v{jj} = r_v{jj} + Matrices.R2_bc{jj}*(Matrices.vel_in{jj}*X.x_va_p1);
    end
    
    n_elem = Matrices.n_elem(jj);
    if Sim.aero_flag
        n_aero = 2;
        S_full(Sim.Ind.flex_ind{jj},Sim.Ind.flex_ind{jj}) = ...
            [Matrices.R1_bc{jj}*[ Matrices.M{jj},                  Sim.struct_damp*Qdamp{jj},     zeros(n_elem*6+6,   n_elem*4), zeros(n_elem*6+6   ,n_elem*3), zeros(n_elem*6+6,n_elem*n_aero)];
             Matrices.R2_bc{jj}*[ zeros(n_elem*6+6,     n_elem*6), Matrices.T1{jj},               zeros(n_elem*6+6,   n_elem*4), zeros(n_elem*6+6   ,n_elem*3), zeros(n_elem*6+6,n_elem*n_aero)];
                                  zeros(n_elem*4,       n_elem*6), zeros(n_elem*4,     n_elem*6), eye(n_elem*4),                 zeros(n_elem*4     ,n_elem*3), zeros(n_elem*4,n_elem*n_aero)   ;
                                  zeros(n_elem*3,       n_elem*6), zeros(n_elem*3,     n_elem*6), zeros(n_elem*3,     n_elem*4), eye(n_elem*3),                 zeros(n_elem*3,n_elem*n_aero)   ;
                                  zeros(n_elem*n_aero,  n_elem*6), zeros(n_elem*n_aero,n_elem*6), zeros(n_elem*n_aero,n_elem*4), zeros(n_elem*n_aero,n_elem*3), eye(  n_elem*n_aero)           ];
        if Matrices.Parent{jj}(1)
            n_elem_p = Matrices.n_elem(parent_ind);
            S_full(Sim.Ind.x_v_ind{parent_ind},Sim.Ind.flex_ind{jj}) = Matrices.R3_bc{jj}*[ Matrices.M{jj},                   Sim.struct_damp*Qdamp{jj}, zeros(n_elem*6+6,     n_elem*4),   zeros(n_elem*6+6     ,n_elem*3),   zeros(n_elem*6+6,  n_elem*n_aero  )];
            S_full(Sim.Ind.x_f_ind{jj},Sim.Ind.flex_ind{parent_ind}) = Matrices.R4_bc{jj}*[ zeros(n_elem_p*6+6, n_elem_p*6) , Matrices.T1{parent_ind},   zeros(n_elem_p*6+6,   n_elem_p*4), zeros(n_elem_p*6+6   ,n_elem_p*3), zeros(n_elem_p*6+6,n_elem_p*n_aero)];
        end
    else
        S_full(Sim.Ind.flex_ind{jj},Sim.Ind.flex_ind{jj}) = ...
            [Matrices.R1_bc{jj}*[ Matrices.M{jj},                  Sim.struct_damp*Qdamp{jj},     zeros(n_elem*6+6,   n_elem*4), zeros(n_elem*6+6   ,n_elem*3), ];
             Matrices.R2_bc{jj}*[ zeros(n_elem*6+6,     n_elem*6), Matrices.T1{jj},               zeros(n_elem*6+6,   n_elem*4), zeros(n_elem*6+6   ,n_elem*3), ];
                                  zeros(n_elem*4,       n_elem*6), zeros(n_elem*4,     n_elem*6), eye(n_elem*4),                 zeros(n_elem*4     ,n_elem*3),  ;
                                  zeros(n_elem*3,       n_elem*6), zeros(n_elem*3,     n_elem*6), zeros(n_elem*3,     n_elem*4), eye(n_elem*3)                  ];
        if Matrices.Parent{jj}(1)
            n_elem_p = Matrices.n_elem(parent_ind);
            S_full(Sim.Ind.x_v_ind{parent_ind},Sim.Ind.flex_ind{jj}) = Matrices.R3_bc{jj}*[ Matrices.M{jj},                   Sim.struct_damp*Qdamp{jj}, zeros(n_elem*6+6,     n_elem*4),   zeros(n_elem*6+6     ,n_elem*3)  ];
            S_full(Sim.Ind.x_f_ind{jj},Sim.Ind.flex_ind{parent_ind}) = Matrices.R4_bc{jj}*[ zeros(n_elem_p*6+6, n_elem_p*6) , Matrices.T1{parent_ind},   zeros(n_elem_p*6+6,   n_elem_p*4), zeros(n_elem_p*6+6   ,n_elem_p*3)];
        end
    end
    
    if Sim.rb_flag
        if Sim.aero_flag
            if ~Matrices.Parent{jj}(1)
                S_fr_full(Sim.Ind.flex_ind{jj},:) = ...
                    [zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                     zeros(n_elem*6,6),      zeros(n_elem*6,4),      zeros(n_elem*6,3);
                     zeros(n_elem*4,6),      zeros(n_elem*4,4),      zeros(n_elem*4,3);
                     zeros(n_elem*3,6),      zeros(n_elem*3,4),      zeros(n_elem*3,3);
                     zeros(n_elem*n_aero,6), zeros(n_elem*n_aero,4), zeros(n_elem*n_aero,3)];
            end
            S_rf_full(:,Sim.Ind.flex_ind{jj}) = ...
                [SYS.M_r{jj}       , zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3), zeros(6,n_elem*n_aero);
                 zeros(4,n_elem*6) , zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3), zeros(4,n_elem*n_aero);
                 zeros(3,n_elem*6) , zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3), zeros(3,n_elem*n_aero)];
        else
            if ~Matrices.Parent{jj}(1)
                S_fr_full(Sim.Ind.flex_ind{jj},:) = ...
                    [zeros(n_elem*6,6), zeros(n_elem*6,4), zeros(n_elem*6,3);
                     zeros(n_elem*6,6), zeros(n_elem*6,4), zeros(n_elem*6,3);
                     zeros(n_elem*4,6), zeros(n_elem*4,4), zeros(n_elem*4,3);
                     zeros(n_elem*3,6), zeros(n_elem*3,4), zeros(n_elem*3,3)];
            end
            S_rf_full(:,Sim.Ind.flex_ind{jj}) = ...
                [SYS.M_r{jj},       zeros(6,n_elem*6), zeros(6,n_elem*4), zeros(6,n_elem*3);
                 zeros(4,n_elem*6), zeros(4,n_elem*6), zeros(4,n_elem*4), zeros(4,n_elem*3);
                 zeros(3,n_elem*6), zeros(3,n_elem*6), zeros(3,n_elem*4), zeros(3,n_elem*3)];
        end
    end

    if Sim.rb_flag
        r_va =   r_va + SYS.Qgyr_r{jj}*X.x_v_p1{jj} - mult*(SYS.F_r{jj}*Matrices.f{jj} + SYS.F_rg{jj}*SYS.f_glob{jj}) - SYS.Btot_r{jj}*(SYS.f{jj}+SYS.f_grav{jj}) - SYS.F_r{jj}*SYS.f_grav_pt{jj};
    end
end
for jj = 1:length(Matrices.n_elem)
    r_full        = [r_full;r_f{jj};r_v{jj};r_q{jj};r_p{jj};r_x{jj}];
end

if Sim.rb_flag
    r_va =   r_va + SYS.omega_tilde*Matrices.M_rb*X.x_va_p1 - SYS.f_grav_r;
    r_qa =          0.5*SYS.OmegaQuat_a*X.x_qa_p1;
    r_pa =             -SYS.CGa_tot    *X.x_va_p1;
    
    r_full = [r_full;r_va;r_qa;r_pa];

    S_rr_full = [Matrices.M_rb+eye(6)*1e-6, zeros(6,4), zeros(6,3);
                 zeros(4,6),                eye(4),     zeros(4,3);
                 zeros(3,6),                zeros(3,4), eye(3)   ];
    
    S_full  = [S_full   , S_fr_full;
               S_rf_full, S_rr_full];
end

vec_out = -S_full\r_full;