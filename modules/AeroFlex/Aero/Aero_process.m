%AERO_PROCESS Post-process for the aerodynamic solver results
%
%   Summary: The following script takes the results and post-processes them
%   to extract the force coefficients and the moments.
%
%   [results,beam_model] = VLM_HS(state,lattice,Gust,Ic,nodevel,GustAct,AIC)
%
%   Inputs:
%       REQUIRED:
%       - FlowVector  = Inflow vector (dim:1x3)
%`      - lattice.VORTEX = Vortex coordinate points (dim:npanels x 4 x 3)
%       - SIMXZ = symmetry flag (dim:1)
%       - lattice.COLLOC
%       - lattice.N
%       OPTIONAL:
%       - AIC = Aerodynamic influence matrix
% 
% 
%   Author: Dario Calderon 

function results = Aero_process(results,lattice,rho,FlowVector,GlobalRot,S_ref,geo,ref_point)


% Recover the Global Forces from the rotation matrix
GlobalForces = (GlobalRot'*results.F')';

% Calculate normal force for each panel
LocalForceMag  = sqrt(results.F(:,1).^2 + results.F(:,2).^2 + results.F(:,3).^2);
GlobalForceMag = sqrt(GlobalForces(:,1).^2 + GlobalForces(:,2).^2 + GlobalForces(:,3).^2);

% Calculate the sum of all the panel forces
LocalForceSum  = sum(results.F,1);
GlobalForceSum = sum(GlobalForces,1);

Velocity = norm(FlowVector);

% Determine the dynamic pressure
q  = 0.5 * rho * Velocity^2;

% Non-dimensionalise the forces local to the panel
CX = LocalForceSum(1)/(q*S_ref);
CY = LocalForceSum(2)/(q*S_ref);
CZ = LocalForceSum(3)/(q*S_ref);

% Calculate the forces in the global frame
CD = GlobalForceSum(1)/(q*S_ref);
CS = GlobalForceSum(2)/(q*S_ref);
CL = GlobalForceSum(3)/(q*S_ref);

% Calculate the normal force coefficient
CN = norm(LocalForceSum)/(q*S_ref);

% Calculate the panel areas
panel_area = tarea(lattice.XYZ);

% Here I seek to work out the lift and area for each spanwise section of
% the wing.
numCaero = length(geo.ny);

sp_mid=(lattice.XYZ(:,1,:)+lattice.XYZ(:,2,:))/2;                  		

countaero   = 0;
countpanels = 0;
countnparts = 0;

spanlift = zeros(sum(geo.ny),1);
spanside = zeros(sum(geo.ny),1);
spandrag = zeros(sum(geo.ny),1);
spanarea = zeros(sum(geo.ny),1);

spanlift_control = zeros(sum(geo.ny) + sum(geo.ny.*logical(geo.fnx)),1);
spanarea_control = zeros(sum(geo.ny) + sum(geo.ny.*logical(geo.fnx)),1);

for i = 1:numCaero
    for j = 1:geo.ny(i)
        
        % Update the caerocount
        countaero   = countaero + 1;
        countnparts = countnparts + 1;
        
        % Number of chordwise panels for the spanwise section
        numcpanels = geo.nx(i) + geo.fnx(i);
        
        % Calculate the total lift for the entire spanwise section
        spanlift(countaero,1) = sum(GlobalForces(countpanels + 1:countpanels + numcpanels,3));
        spanside(countaero,1) = sum(GlobalForces(countpanels + 1:countpanels + numcpanels,2));
        spandrag(countaero,1) = sum(GlobalForces(countpanels + 1:countpanels + numcpanels,1));
        
        % Calculate the area for the spanwise location
        spanarea(countaero,1) = sum(panel_area(countpanels + 1:countpanels + numcpanels));

        % Calculate the lift for the spanwise section excluding the control
        % surface
        spanlift_control(countnparts,1) = sum(GlobalForces(countpanels + 1:countpanels + geo.nx(i),3));
       
        % Calculate the corresponding area
        spanarea_control(countnparts,1) = sum(panel_area(countpanels + 1:countpanels + geo.nx(i)));

        if geo.fnx(i) > 0
            countnparts = countnparts + 1;
            spanlift_control(countnparts,1) = ...
                sum(GlobalForces(countpanels + geo.nx(i) + 1:countpanels + geo.nx(i) + geo.fnx(i),3));
            
            spanarea_control(countnparts,1) = ...
                sum(panel_area(countpanels + geo.nx(i) + 1:countpanels + geo.nx(i) + geo.fnx(i)));
        end
        
        p_mid(countaero,1:3) = [sum(sp_mid(countpanels + 1:countpanels + numcpanels,:,1)),...
            sum(sp_mid(countpanels + 1:countpanels + numcpanels,:,2)),...
            sum(sp_mid(countpanels + 1:countpanels + numcpanels,:,3))]/numcpanels;
        
        
        %Update the panel count
        countpanels = countpanels + numcpanels;
    end
end


results.p_mid_r = sqrt(p_mid(:,2).^2);

% Panel vortex midpoint
wpsize  = size(lattice.COLLOC,1);
vlength = size(lattice.VORTEX,2);

b1 = vlength/2;

p1(:,:) = lattice.VORTEX(1:wpsize,b1,:);
p2(:,:) = lattice.VORTEX(1:wpsize,b1+1,:);

% Determine the moments about the reference point
MidVort = (p1+p2)./2;
MomArm  = MidVort-ones(size(MidVort,1),1)*ref_point;

Mlocal = cross(MomArm,results.F,2); 

% Store the results back into the results structure:
results.CX = CX;
results.CY = CY;
results.CZ = CZ;
results.CN = CN;
results.CD = CD;
results.CS = CS;
results.CL = CL;

results.spancl = spanlift./(q*spanarea);
results.spancs = spanside./(q*spanarea);
results.spancd = spandrag./(q*spanarea);

results.spanarea = spanarea;

results.spancl_control   = spanlift_control./(q*spanarea_control);
results.spanlift_control = spanlift_control;

results.spanlift = spanlift;
results.spanside = spanside;
results.spandrag = spandrag;

results.Fglobal = GlobalForces;

results.L = GlobalForceSum(3);
results.D = GlobalForceSum(1);
results.C = GlobalForceSum(2);

results.GlobalN = GlobalForceMag;
results.LocalN  = LocalForceMag;

results.M = Mlocal;

results.TotalForce  = GlobalForceSum;
results.TotalMoment = sum(Mlocal,1);

results.FORCES = GlobalForceSum;
results.MOMENT = sum(Mlocal,1);

end%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [panel_area]=tarea(XYZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tarea: Subsidary function for TORNADO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the area of each panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Department of Aeronautics	%
%				Copyright 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidaty function for TORNADO
% Called by:	coeff_create
%
% Calls:			MATLAB 5.2 std fcns
% Loads:	none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a b c]=size(XYZ);
for i=1:a
    p1=[XYZ(i,1,1) XYZ(i,1,2) XYZ(i,1,3)];	%sets up the vectors
    p2=[XYZ(i,2,1) XYZ(i,2,2) XYZ(i,2,3)];	%to the corners of the
    p3=[XYZ(i,3,1) XYZ(i,3,2) XYZ(i,3,3)];	%panel.
    p4=[XYZ(i,4,1) XYZ(i,4,2) XYZ(i,4,3)];
    
    a=p2-p1;	%sets up the edge vectors
    b=p4-p1;
    c=p2-p3;
    d=p4-p3;
    
    %ar1=norm(cross(b,a))/2;	%claculates the ctoss product of
    %ar2=norm(cross(c,d))/2;	%two diagonal corners
    
    % Added by Dario Calderon 2016
    ar1 = 0.5*((b(2).*a(3)-b(3).*a(2))^2 + (b(3).*a(1)-b(1).*a(3))^2 + (b(1).*a(2)-b(2).*a(1))^2)^0.5;
    ar2 = 0.5*((c(2).*d(3)-c(3).*d(2))^2 + (c(3).*d(1)-c(1).*d(3))^2 + (c(1).*d(2)-c(2).*d(1))^2)^0.5;
    
    panel_area(i)=ar1+ar2;	%Sums up the product to make the
end						    %Area
end% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


