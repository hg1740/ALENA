function [OPTIM,flag,output] = run_optimisation(OPTIM)


XL = [OPTIM.skin.min*ones(OPTIM.Nsec,1);...
    OPTIM.spar.min*ones(OPTIM.Nsec,1);...
    OPTIM.stringer.min*ones(OPTIM.Nsec,1)];

XU = [OPTIM.skin.max*ones(OPTIM.Nsec,1);...
    OPTIM.spar.max*ones(OPTIM.Nsec,1);...
    OPTIM.stringer.max*ones(OPTIM.Nsec,1)];


X0 = OPTIM.X0;

TolCon = 0.00005;

% I NEED TO INCLUDE THE AERODYNAMICS

Optim_opt = optimset('Algorithm','sqp','GradObj', 'off','Display', 'off',...
    'MaxFunEvals',500000,'TolCon', TolCon, 'LargeScale', 'off', 'tolfun',1e-9);%,'OutputFcn',@outfun);%,'PlotFcns',@optimplotx);

% FEED THE AERODYNAMIC SOLVER HERE??

% based on the thickness, the beam reduction function gives the beam
% properties required to run and aeroelastic analysis. The loads are then
% taken from that and fed into the optimisation routine
try
    
    
    [SOL, fun, flag, output] = fmincon(@(X)ObjFunc(X,OPTIM), X0, ...
        [],[],[],[], ...
        XL, XU, @(X)SectionConstraintsStr(X,OPTIM), Optim_opt);
    OPTIM.PlotStress = 1;
    OPTIM.CombinedPlot = 1;
    [CIN, ~] = SectionConstraintsStr(SOL, OPTIM);
    OPTIM.PlotStress = 0;
%     fprintf('\n\tOptimisation solution');
%     fprintf('\n\t- Root skin thickness = %-.2f mm',SOL(1)*1000);
%     fprintf('\n\t- Tip skin thickness = %-.2f mm',SOL(OPTIM.Nsec)*1000);
%     fprintf('\n\t- Root spar thickness = %-.2f mm',SOL(OPTIM.Nsec+1)*1000);
%     fprintf('\n\t- Tip spar thickness = %-.2f mm',SOL(2*OPTIM.Nsec)*1000);
%     fprintf('\n\t- Root Stringer Area = %-.2f m^2',SOL(2*OPTIM.Nsec+1));
%     fprintf('\n\t- Tip Stringer Area = %-.2f m^2\n',SOL(end));
    
catch
    
    flag = -1;
    fun = Inf;
    output = [];
    SOL = X0;
    CIN = [];
    error('WARNING - Optimisation was not performed\n');
    
end

OPTIM.SOL  = SOL;
OPTIM.fval = fun;
OPTIM.CIN  = CIN;
end


function stop = outfun(x,optimValues,state)
stop = false;
switch state
    case 'init'
    case 'iter'
        dlmwrite('fminConOut.out',[optimValues.iteration,optimValues.fval,optimValues.constrviolation,x'],'-append','delimiter', ',', 'precision', 6);
    case 'done'
    otherwise
end

end