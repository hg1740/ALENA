function X = initState(x0,Matrices,Aero,Sim, AnalysisParam)
%initState Generates the initial state vector.

if strcmpi(AnalysisParam.AnalysisType, 'static') %Static
    
    %Velocities and forces map straight across
    X.x_f  = x0.x_f;
    X.x_va = x0.x_va;
    
    return
        
end

if Sim.Soln == 2 || Sim.Soln == 2.1 
    X.x_f     = cell(length(Matrices.n_elem),1);
    X.x_v     = cell(length(Matrices.n_elem),1);
    X.x_q     = cell(length(Matrices.n_elem),1);
    X.x_p     = cell(length(Matrices.n_elem),1);
    X.x_x     = cell(length(Matrices.n_elem),1);
    
    X.x_f_dot = cell(length(Matrices.n_elem),1);
    X.x_v_dot = cell(length(Matrices.n_elem),1);
    X.x_q_dot = cell(length(Matrices.n_elem),1);
    X.x_p_dot = cell(length(Matrices.n_elem),1);
    X.x_x_dot = cell(length(Matrices.n_elem),1);
    
    if Sim.aero_flag
        b1      = Aero.Leishman.b1;
        b2      = Aero.Leishman.b2;
        
        rho         = Aero.rho;
        stall_angle = Aero.stall_angle;
    end
    
    for jj = 1:length(Matrices.n_elem)
        X.x_f{jj}      = zeros(length(x0.x_f{jj}),length(Sim.time_vec));
        X.x_v{jj}      = zeros(length(x0.x_v{jj}),length(Sim.time_vec));
        X.x_q{jj}      = zeros(length(x0.x_q{jj}),length(Sim.time_vec));
        X.x_p{jj}      = zeros(length(x0.x_p{jj}),length(Sim.time_vec));
        
        X.x_f{jj}(:,1) = x0.x_f{jj};
        X.x_v{jj}(:,1) = x0.x_v{jj};
        X.x_q{jj}(:,1) = x0.x_q{jj};
        X.x_p{jj}(:,1) = x0.x_p{jj};
        
        X.x_f_dot{jj}  = zeros(length(x0.x_f{jj}),length(Sim.time_vec));
        X.x_v_dot{jj}  = zeros(length(x0.x_v{jj}),length(Sim.time_vec));
        X.x_q_dot{jj}  = zeros(length(x0.x_q{jj}),length(Sim.time_vec));
        X.x_p_dot{jj}  = zeros(length(x0.x_p{jj}),length(Sim.time_vec));
        
        R1_bc{jj} = sparse(zeros(length(x0.x_f{jj}),length(x0.x_f{jj})+6));
        R2_bc{jj} = sparse(zeros(length(x0.x_f{jj}),length(x0.x_f{jj})+6));
        R1_bc{jj}(:,7:end)   = eye(length(x0.x_f{jj}));
        R2_bc{jj}(:,1:end-6) = eye(length(x0.x_f{jj}));
        
        for ii = 1:Matrices.n_elem(jj)
            ind3 = [1:3] + (ii-1)*3;
            ind4 = [1:4] + (ii-1)*4;
            ind6 = [1:6] + (ii-1)*6;
            Ct0  = Quat2Rot(x0.x_q{jj}(ind4));
            X.x_p_dot{jj}(ind3,1) = Ct0*x0.x_v{jj}(ind6(1:3));
        end
        
        if Sim.aero_flag
            if ~isfield(x0,'x_x')
                A_x{jj}     = zeros((Matrices.n_elem(jj))*2,(Matrices.n_elem(jj))*2);
                A_x_1{jj}   = zeros((Matrices.n_elem(jj))*2,1);
                
                for ii = 1:Matrices.n_elem(jj)
                    ind3 = [1:6] + (ii-1)*6;
                    ind5 = [1:2] + (ii-1)*2;
                    
                    v_A0       = blkdiag(Matrices.CBA0{jj}(:,:,ii)',Matrices.CBA0{jj}(:,:,ii)')*[eye(3) -skew(Matrices.cp{jj}(:,ii)); zeros(3) eye(3)]*X.x_v{jj}(ind3,1);
                    
                    ydot =  v_A0(2);
                    zdot =  v_A0(3);
                    adot =  v_A0(4);
                    
                    if ydot == 0
                        ydot = 1e-10;
                    end
                    
                    chord         = Matrices.chord{jj}(ii);
                    
                    A_x{jj}(ind5,ind5) = -2*ydot/chord*diag([b1 b2]);
                    A_x_1{jj}(ind5,1)  = - zdot/ydot*[1;1];
                end
                x_x0{jj}   = -A_x{jj}\A_x_1{jj};
            end
            
            X.x_x{jj}      = zeros(length(x_x0{jj}),length(Sim.time_vec));
            X.x_x{jj}(:,1) = x_x0{jj};
            X.x_x_dot{jj}  = zeros(length(x_x0{jj}),length(Sim.time_vec));
        else
            X.x_x{jj}      = zeros(0,length(Sim.time_vec));
            X.x_x_dot{jj}  = zeros(0,length(Sim.time_vec));
        end
        
%         Q_BC_r{jj}     = sparse(zeros(6,1));
    end
    
    if Sim.rb_flag
        X.x_va     = zeros(6,length(Sim.time_vec));
        X.x_qa     = zeros(4,length(Sim.time_vec));
        X.x_pa     = zeros(3,length(Sim.time_vec));
        
        X.x_va_dot = zeros(6,length(Sim.time_vec));
        X.x_qa_dot = zeros(4,length(Sim.time_vec));
        X.x_pa_dot = zeros(3,length(Sim.time_vec));
 
        X.x_va(:,1) = x0.x_va;
        X.x_qa(:,1) = x0.x_qa;
        X.x_pa(:,1) = x0.x_pa;
        
        Ct0             = Quat2Rot(x0.x_qa);
        X.x_pa_dot(:,1) = Ct0*x0.x_va(1:3);
    else
        X.x_va = x0.x_va;
        [t,QuatPos] = ode45(@(t,y) getQuatPos(t,y,x0.x_va,Sim.time_vec), Sim.time_vec, [x0.x_qa;x0.x_pa]);
        X.x_qa = QuatPos(:,1:4)';
        X.x_pa = QuatPos(:,5:7)';
    end

%     X_p1 = stepState(        X,     Matrices,Aero,Sim,1);
%     SYS  = getSystemMatrices(X_p1,  Matrices,Aero,Sim,0);  
%     X    = getInitialRates(  X_p1,X,Matrices,SYS, Sim  );
    
    KE        = zeros(length(Sim.time_vec),1);
    PE        = zeros(length(Sim.time_vec),1);
    residual  = zeros(length(Sim.time_vec),1);
    
    return
    
end

if Sim.Soln == 2.2
    X.x      = zeros(length(x0.x),length(Sim.time_vec));
    X.x(:,1) = x0.x;
    X.x_dot  = zeros(length(x0.x),length(Sim.time_vec));    
end

end