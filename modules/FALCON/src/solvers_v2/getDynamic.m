function x_out = getDynamic(x0,Matrices,Aero,Sim,ParForVars)

parfor_switch = 0;
if exist('ParForVars','var') 
    fprintf('\nUpdating inputs for Parfor loop!\n')
    parfor_switch = 1;
    if~isempty(ParForVars) && isstruct(ParForVars)
        [x0,Matrices,Aero,Sim] = updateParForVars(x0,Matrices,Aero,Sim,ParForVars);
    end
end

Sim = initSystem(x0,Matrices,Aero,Sim);
X   = initState( x0,Matrices,Aero,Sim);

for tt = 1:length(Sim.time_vec)-1
    if ~parfor_switch
        fprintf(['\nSOLVING FOR TIME STEP ',num2str(Sim.time_vec(tt)),'s to ',num2str(Sim.time_vec(tt+1)),'s...\n'])
    end
    
    X_p1 = stepState(X,Matrices,Aero,Sim,tt);
        
    iter_counter = 0;
    while 1
        Sim.tangent_flag = ~mod(iter_counter,Sim.tangent_update);
        
        SYS            = getSystemMatrices(X_p1,Matrices,Aero,Sim,iter_counter); 
        [r_full,Q,S,D] = getResidual(X_p1,Matrices,SYS,Sim,tt); 
        
        % Save the tangent matrix if it isn't being updated
        if Sim.tangent_flag
            S_sparse = S;
            D_scale  = D;
        end
        
        if getResidualNorm(r_full)<Sim.eps
            X = updateResults(X,X_p1,SYS,Matrices,r_full,Sim,tt);
            break
        end
        
        X_p1 = updateState(X_p1,Q,Matrices,S_sparse,r_full,D_scale,Sim);
        
        iter_counter = iter_counter + 1;
    end
end

x_out = X;