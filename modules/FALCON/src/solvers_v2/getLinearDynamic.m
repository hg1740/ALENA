function x_out = getLinearDynamic(x0,Matrices,Aero,Sim)

Sim               = initSystem(x0,Matrices,Aero,Sim);
X                 = initState( x0,Matrices,Aero,Sim);
Sim.time.dt       = 0;
Sim.speed_up      = 0;
if ~isfield(Sim,'method')
    Sim.method = 1;
end
X_p1              = stepState(X,Matrices,Aero,Sim,1);

if Sim.method == 1 || Sim.method == 2 %Calculate the linearisation analytically
    SYS               = getSystemMatrices(X_p1,Matrices,Aero,Sim,1);
    [A,B1,B2,C,D1,D2] = getDynLinSys(X,Matrices,SYS,Sim);
elseif Sim.method == 3
    [A,B1,B2,C,D1,D2] = getDynLinSysNum(X_p1,Matrices,Aero,Sim);
end

x_out.A  = A;
x_out.B1 = B1;
x_out.B2 = B2;
x_out.C  = C;
x_out.D1 = D1;
x_out.D2 = D2;