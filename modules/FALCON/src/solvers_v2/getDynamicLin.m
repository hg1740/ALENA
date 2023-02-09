function x_out = getDynamicLin(A,B,C,D,u,Sim)

if ~isfield(Sim,'rho_inf')
    Sim.rho_inf = 0.75;
end
alphm1      = 0.5*(3-Sim.rho_inf)/(1+Sim.rho_inf);
alphf1      = 1/(1+Sim.rho_inf);
Sim.gam1    = 0.5 + alphm1 - alphf1;
Sim.time_vec = Sim.time.limits(1):Sim.time.dt:Sim.time.limits(end);
Sim.tangent_update = 100;

X    = zeros(length(A),length(Sim.time_vec));
Xdot = zeros(length(A),length(Sim.time_vec));

for tt = 1:length(Sim.time_vec)-1
    fprintf(['\nSOLVING FOR TIME STEP ',num2str(Sim.time_vec(tt)),'s to ',num2str(Sim.time_vec(tt+1)),'s...\n'])
    
    X_p1     = X(:,tt) + Sim.time.dt*(1-Sim.gam1) * Xdot(:,tt);
    X_dot_p1 = zeros(length(A),1);
    
    iter_counter = 0;
    while 1
        Sim.tangent_flag = ~mod(iter_counter,Sim.tangent_update);

        r_full = X_dot_p1 - A*X_p1 - B*u(:,tt+1);
        
        % Save the tangent matrix if it isn't being updated
        if Sim.tangent_flag
            S_sparse = eye(size(A))/Sim.gam1/Sim.time.dt - A;
            D_scale  = sparse(diag(1./max(abs(S_sparse),[],2)));
        end
        
        if getResidualNorm(r_full)<Sim.eps
            X(:,tt+1) = X_p1;
            break
        end
        
        dq = - (D_scale*S_sparse)\(D_scale*r_full);
%         dq = - S_sparse\r_full;
        
        X_p1     = X_p1     + 1.0*dq;
        X_dot_p1 = X_dot_p1 + 1.0*dq/Sim.gam1/Sim.time.dt;
        
        iter_counter = iter_counter + 1;
    end
end

x_out = X;