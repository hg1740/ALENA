function X_p1 = updateState(X_p1,Q,Matrices,S_sparse,r_full,D_scale,Sim)

if Sim.Soln == 1

elseif Sim.Soln == 2
    dq = - (D_scale*S_sparse)\(D_scale*r_full);
    
    Q.q_p1_full     = Q.q_p1_full     + 1.0*dq;
    Q.q_dot_p1_full = Q.q_dot_p1_full + 1.0*dq/Sim.gam1/Sim.time.dt;
    
    for jj = 1:length(Matrices.n_elem)
        X_p1.x_f_p1{jj}     = Q.q_p1_full(Sim.Ind.x_f_ind{jj},1);
        X_p1.x_v_p1{jj}     = Q.q_p1_full(Sim.Ind.x_v_ind{jj},1);
        X_p1.x_q_p1{jj}     = Q.q_p1_full(Sim.Ind.x_q_ind{jj},1);
        X_p1.x_p_p1{jj}     = Q.q_p1_full(Sim.Ind.x_p_ind{jj},1);
        X_p1.x_x_p1{jj}     = Q.q_p1_full(Sim.Ind.x_x_ind{jj},1);
        
        for ii = 1:Matrices.n_elem(jj)
            ind1 = [1:4] + (ii-1)*4;
            X_p1.x_q_p1{jj}(ind1) = X_p1.x_q_p1{jj}(ind1)/norm(X_p1.x_q_p1{jj}(ind1));
        end
        
        X_p1.x_f_dot_p1{jj} = Q.q_dot_p1_full(Sim.Ind.x_f_ind{jj},1);
        X_p1.x_v_dot_p1{jj} = Q.q_dot_p1_full(Sim.Ind.x_v_ind{jj},1);
        X_p1.x_q_dot_p1{jj} = Q.q_dot_p1_full(Sim.Ind.x_q_ind{jj},1);
        X_p1.x_p_dot_p1{jj} = Q.q_dot_p1_full(Sim.Ind.x_p_ind{jj},1);
        X_p1.x_x_dot_p1{jj} = Q.q_dot_p1_full(Sim.Ind.x_x_ind{jj},1);
    end
    
    if Sim.rb_flag
        X_p1.x_va_p1     = Q.q_p1_full(    Sim.Ind.x_va_ind,1);
        X_p1.x_qa_p1     = Q.q_p1_full(    Sim.Ind.x_qa_ind,1);
        X_p1.x_pa_p1     = Q.q_p1_full(    Sim.Ind.x_pa_ind,1);
        X_p1.x_va_dot_p1 = Q.q_dot_p1_full(Sim.Ind.x_va_ind,1);
        X_p1.x_qa_dot_p1 = Q.q_dot_p1_full(Sim.Ind.x_qa_ind,1);
        X_p1.x_pa_dot_p1 = Q.q_dot_p1_full(Sim.Ind.x_pa_ind,1);
        
        X_p1.x_qa_p1     = X_p1.x_qa_p1/norm(X_p1.x_qa_p1);
    end
elseif Sim.Soln == 2.1
    dq1 = - (D_scale{1}*S_sparse{1})\(D_scale{1}*r_full{1});
    
    Q.q_p1_full{1}     = Q.q_p1_full{1}     + 1.0*dq1;
    Q.q_dot_p1_full{1} = Q.q_dot_p1_full{1} + 1.0*dq1/Sim.gam1/Sim.time.dt;
    
    dq2 = - (D_scale{2}*S_sparse{2})\(D_scale{2}*r_full{2});
    
    Q.q_p1_full{2}     = Q.q_p1_full{2}     + 1.0*dq2;
    Q.q_dot_p1_full{2} = Q.q_dot_p1_full{2} + 1.0*dq2/Sim.gam1/Sim.time.dt;
    
    for jj = 1:length(Matrices.n_elem)
        X_p1.x_f_p1{jj}     = Q.q_p1_full{1}(Sim.Ind.x_f_ind{jj},1);
        X_p1.x_v_p1{jj}     = Q.q_p1_full{1}(Sim.Ind.x_v_ind{jj},1);
        X_p1.x_q_p1{jj}     = Q.q_p1_full{2}(Sim.Ind.x_q_ind{jj},1);
        X_p1.x_p_p1{jj}     = Q.q_p1_full{2}(Sim.Ind.x_p_ind{jj},1);
        X_p1.x_x_p1{jj}     = Q.q_p1_full{1}(Sim.Ind.x_x_ind{jj},1);
        
        for ii = 1:Matrices.n_elem(jj)
            ind1 = [1:4] + (ii-1)*4;
            X_p1.x_q_p1{jj}(ind1) = X_p1.x_q_p1{jj}(ind1)/norm(X_p1.x_q_p1{jj}(ind1));
        end
        
        X_p1.x_f_dot_p1{jj} = Q.q_dot_p1_full{1}(Sim.Ind.x_f_ind{jj},1);
        X_p1.x_v_dot_p1{jj} = Q.q_dot_p1_full{1}(Sim.Ind.x_v_ind{jj},1);
        X_p1.x_q_dot_p1{jj} = Q.q_dot_p1_full{2}(Sim.Ind.x_q_ind{jj},1);
        X_p1.x_p_dot_p1{jj} = Q.q_dot_p1_full{2}(Sim.Ind.x_p_ind{jj},1);
        X_p1.x_x_dot_p1{jj} = Q.q_dot_p1_full{1}(Sim.Ind.x_x_ind{jj},1);
    end
    
    X_p1.x_va_p1     = Q.q_p1_full{1}(    Sim.Ind.x_va_ind,1);
    X_p1.x_qa_p1     = Q.q_p1_full{2}(    Sim.Ind.x_qa_ind,1);
    X_p1.x_pa_p1     = Q.q_p1_full{2}(    Sim.Ind.x_pa_ind,1);
    X_p1.x_va_dot_p1 = Q.q_dot_p1_full{1}(Sim.Ind.x_va_ind,1);
    X_p1.x_qa_dot_p1 = Q.q_dot_p1_full{2}(Sim.Ind.x_qa_ind,1);
    X_p1.x_pa_dot_p1 = Q.q_dot_p1_full{2}(Sim.Ind.x_pa_ind,1);
    
    X_p1.x_qa_p1     = X_p1.x_qa_p1/norm(X_p1.x_qa_p1);
elseif Sim.Soln == 2.2
    dq = - S_sparse\r_full;
    
    X_p1.x_p1       = Q.q_p1     + 1.0*dq;
    X_p1.x_dot_p1   = Q.q_dot_p1 + 1.0*dq/Sim.gam1/Sim.time.dt;
end