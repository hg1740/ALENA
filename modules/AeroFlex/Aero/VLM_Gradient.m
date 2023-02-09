% The following script calculates the Aerodynamic gradient matrix for the
% VLM solver
% dFs_dq = the gradient of forces and moments of the states with respect to
% the displacements
% dFa_dq = the gradient of the forces on the panels with respect to the
% displacements
% 
% 
%   Author: Dario Calderon 

function [dFs_dq,dFa_dq,dFarot_dq] = VLM_Gradient(state,states,beam_model,lattice)

%fprintf('\n Aerodynamic Gradient Matrix ... ');
% Get a few useful counters
npanels = size(beam_model.Aero.lattice_vlm.COLLOC,1);

% Determine the Vector that runs along the 1/4chord of the panel
VLength  = size(lattice.VORTEX,2);
b1       = VLength/2;
p1(:,:)  = lattice.VORTEX(:,b1,:);		
p2(:,:)  = lattice.VORTEX(:,b1+1,:);	
LEVector = (p2-p1);

% Use the states to create a matrix that includes the aerodynamic splining
% nodes:
Downwash = sparse(npanels,6*npanels);
SpanCirc = sparse(3*npanels,npanels);
CrossWind = sparse(3*npanels,3*npanels);
for n = 1:npanels
    Downwash(n,6*(n-1)+5) = state.AS*lattice.N(n,3);
    Downwash(n,6*(n-1)+6) = -state.AS*lattice.N(n,2);
    SpanCirc(3*(n-1)+1:3*(n-1)+3,n) = LEVector(n,:)'; % This needs to be spanwise vector
    %CrossWind(3*(n-1)+1:3*(n-1)+3,3*(n-1)+1:3*(n-1)+3) = state.rho*state.AS*crossm([1,0,0]);
    CrossWind(3*(n-1)+1:3*(n-1)+3,3*(n-1)+1:3*(n-1)+3) = crossm([1,0,0]);
end

% Lift for the aerodynamic panels
[AIC, ~] = fastdw(lattice,[]);

%invAIC = inv(AIC);

%dDWN_dq = Downwash*states.InterpCol*states.DOF2Str;

%dGamma_dq = AIC\(Downwash*states.InterpCol*(states.DOF2Str+states.SkAeroArm*states.ROT2TRA));
dGamma_dq = AIC\(Downwash*states.InterpCol*(states.DOF2Str));

% dummy_alpha = zeros(768,1);
% 
% for i = 5:6:768
%    dummy_alpha(i) = pi*10/180; 
% end

dFa_dq    = state.rho*state.AS*CrossWind*SpanCirc*dGamma_dq;
dFarot_dq = states.CrossCGArm*dFa_dq;

dFs_dq = states.InterpSpline*(states.InterpFor'*dFa_dq + states.SumAeroArm*states.SkAeroArm*states.InterpFor'*dFa_dq);


end