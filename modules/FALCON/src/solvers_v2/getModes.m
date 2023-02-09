function x_out = getModes(x0,Matrices,Aero,Sim)

Sim   = initSystem(x0,Matrices,Aero,Sim);
X     = initState( x0,Matrices,Aero,Sim);
X_p_1 = stepState(  X,Matrices,Aero,Sim,1);

SYS            = getSystemMatrices(X_p_1,Matrices,Aero,Sim,1);
[A,~,C,~]      = getDynLinSys(     X_p_1,Matrices,SYS, Sim  );

[vecs,vals] = eig(full(A));%(1:end-5,1:end-5)
vals        = diag(vals);
[vals,ind]  = sort(vals);
vecs        = vecs(:,ind);
% scatter(real(vals),imag(vals)/2/pi,'rx')

frq         = sort(abs(imag(vals)))/2/pi;

x_out.vecs  = vecs;
x_out.vals  = vals;

x_out.frq   = frq;