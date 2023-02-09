function X = stepState(x0,Matrices,Aero,Sim, AnalysisParam, tt)

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %x_f.p1 is the current estimate of the states
        X.x_f_p1  = x0.x_f;
        X.x_va_p1 = x0.x_va;
        
    otherwise
        
        if Sim.Soln == 2 || Sim.Soln == 2.1 %Dynamic
            X.x_f_p1     = cell(length(Matrices.n_elem),1);
            X.x_v_p1     = cell(length(Matrices.n_elem),1);
            X.x_q_p1     = cell(length(Matrices.n_elem),1);
            X.x_p_p1     = cell(length(Matrices.n_elem),1);
            X.x_x_p1     = cell(length(Matrices.n_elem),1);
            
            X.x_f_dot_p1 = cell(length(Matrices.n_elem),1);
            X.x_v_dot_p1 = cell(length(Matrices.n_elem),1);
            X.x_q_dot_p1 = cell(length(Matrices.n_elem),1);
            X.x_p_dot_p1 = cell(length(Matrices.n_elem),1);
            X.x_x_dot_p1 = cell(length(Matrices.n_elem),1);
            
            for jj = 1:length(Matrices.n_elem)
                X.x_f_p1{jj} = x0.x_f{jj}(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_f_dot{jj}(:,tt);
                X.x_v_p1{jj} = x0.x_v{jj}(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_v_dot{jj}(:,tt);
                X.x_q_p1{jj} = x0.x_q{jj}(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_q_dot{jj}(:,tt);
                X.x_p_p1{jj} = x0.x_p{jj}(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_p_dot{jj}(:,tt);
                X.x_x_p1{jj} = x0.x_x{jj}(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_x_dot{jj}(:,tt);
                
                X.x_f_dot_p1{jj} = zeros(size(X.x_f_p1{jj}));
                X.x_v_dot_p1{jj} = zeros(size(X.x_v_p1{jj}));
                X.x_q_dot_p1{jj} = zeros(size(X.x_q_p1{jj}));
                X.x_p_dot_p1{jj} = zeros(size(X.x_p_p1{jj}));
                X.x_x_dot_p1{jj} = zeros(size(X.x_x_p1{jj}));
            end
            if Sim.rb_flag
                X.x_va_p1 = x0.x_va(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_va_dot(:,tt);
                X.x_qa_p1 = x0.x_qa(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_qa_dot(:,tt);
                X.x_pa_p1 = x0.x_pa(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_pa_dot(:,tt);
                
                X.x_va_dot_p1 = zeros(size(X.x_va_p1));
                X.x_qa_dot_p1 = zeros(size(X.x_qa_p1));
                X.x_pa_dot_p1 = zeros(size(X.x_pa_p1));
            else
                X.x_va_p1 = x0.x_va(:,tt+1);
                X.x_qa_p1 = x0.x_qa(:,tt+1);
                X.x_pa_p1 = x0.x_pa(:,tt+1);
            end
        elseif Sim.Soln == 2.2
            X.x_p1     = x0.x(:,tt) + Sim.time.dt*(1-Sim.gam1) * x0.x_dot(:,tt);
            X.x_dot_p1 = zeros(size(X.x_p1));
        end
        
end

end