%OBJFUNC: The objective function used in the optimisation routine that
%allows it to calculate the mass based on the current values of the
%solution.
%
%   Summary: I will hopefully add the aerodynamic loads to this objective
%   function that allows the solver to work out the aerodynamic loads. 

function [M] = ObjFunc(X, data)

Nsec = data.Nsec;

if data.Type == 1
    tskinroot = X(1);
    tskinttip = X(2);
    tsparroot = X(3);
    tspartip  = X(4);
    % Coefskin = X(5);
    % Coefspar = X(6);
    
    x_coord = linspace(0,data.y_ends(end),Nsec+1);
    
    C = data.c_box;
    H = data.h_spar;
    
    % Interpolate thicknesses and grab lengths to work out the mass
    % Interpolate thicknesses and grab lengths to work out the mass
    [t_skin1] = Lineardist(tskinroot,tskinttip,x_coord,data.y_ends(end));
    [t_spar1] = Lineardist(tsparroot,tspartip,x_coord,data.y_ends(end));
    
    % Is this needed?
    tskin = interp1(data.y_ends,t_skin1',data.y_mbox,'Linear');
    tspar = interp1(data.y_ends,t_spar1',data.y_mbox,'Linear');

elseif data.Type == 2
    tskin = X(1:Nsec);
    tspar = X(Nsec+1:end);
    C = data.c_box;
    H = data.h_spar;
    A = 2*(tskin.*C + tspar.*(H - 2*tskin));
    
elseif data.Type == 3
    tskin = X(1:Nsec);
    tspar = X(Nsec+1:2*Nsec);
    Astr = X(2*Nsec+1:end);
    C = data.c_box;
    H = data.h_spar;
    A = 2*(data.NS.*Astr + tskin.*C + tspar.*(H-2*tskin));
    
elseif data.Type == 4
    tskin = X(1:Nsec);
    tstr = X(2*Nsec+1:end);
    %tspar = data.Optim.geo.Wing.t_spar;
    tspar = X(Nsec+1:2*Nsec);
    C = data.c_box;
    H = data.h_spar;
    A = 2*(tstr.*C + tskin.*C + tspar.*(H-2*(tskin+tstr)));
end

% Mass
M = sum(data.Rho*A.*data.l_box);

% % Gradient
% G = [2*C; 2*H; 2*Ns];

%Breg_Range = V*(L/D)*log(FB/LGW + 1)/TSFC;

%ObjOut = data.Breguet.beta*1/Breg_Range + (1-data.Breguet.beta)*M;


end