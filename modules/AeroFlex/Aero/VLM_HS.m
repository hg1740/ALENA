%VLM_HS Vortex Lattice Method Solver (Horseshoe Vortices)
%
%   Summary: The following VLM method assumes that the bound vortices are
%   expressed as horseshoe vortices. The treatment of the frames is very
%   important. Two different treatments exist here. For the static analysis
%   the Velocity of the flow is rotated with respect to the aerodynamic
%   mesh, however, for the dynamic analysis, the mesh is assumed to be
%   rotated with respect to the flow. This is very annoying and needs to be
%   changed, however the rest of the framework needs to be updated
%   accordingly.
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
%
% File History
%   - 09/10/19 C.Szczyglowski
%       Removed for-loop for 'InFlow'.

function results = VLM_HS(FlowVector,rho,lattice,SIMXZ,velocities,AIC,DW,IndDrag,IndDownwash,DW_sym)
if nargin < 10
    DW_sym = [];
end
% Calculate the panel velocities
[wpsize, vor_length, ~] = size(lattice.VORTEX);

% Calculate the panel velocities
Inflow = zeros(wpsize,3);

% Check if the panel velocities are an input to the function
if isempty(velocities)
    velocities = zeros(wpsize,3);
end

% Assume that the inflow can either come from the velocities array or the
% flow vector. The latter would more appropriate for a static analysis,
% where the local velocities are constant.
% for i = 1:wpsize
%     Inflow(i,:) = FlowVector + velocities(i,:);
% end
Inflow = FlowVector + velocities;

% Do not calculate the AIC matrix if the AIC matrix already exists. To get
% around this, especially when dealing with a nonlinear problem, simply set
% it to the empty matrix outside of this function.
if isempty(AIC)
    
    [AIC, DW] = fastdw(lattice,DW);
    
    if SIMXZ
        % TODO CHECK THIS!!! - it should be the vortex points and not the
        % collocation points that should be rotated
        lattice_dummy             = lattice;
        lattice_dummy.COLLOC(:,2) = -lattice_dummy.COLLOC(:,2);
        %lattice_dummy.VORTEX(:,:,2) = - lattice_dummy.VORTEX(:,:,2);
        lattice_dummy.N(:,2)      = -lattice_dummy.N(:,2);
        [wdummy, DW_sym]          = fastdw(lattice_dummy,DW_sym);
        AIC                       = AIC + wdummy;
    end
    
end

% Calculate the newmann boundary velocity
RHS = sum(lattice.N.*Inflow,2);

% Calculate the vortex strength
gamma = AIC\RHS;

% Determine the inboard leading edge vortex point. This allows us to use
% both a horseshoe method and a Vortex ring method
% Horseshoe: b1 = 3, Vortex Ring: b1 = 4
%   - Note (C.Szczyglowski 09/10/2019) Points are numbered starting at the
%   inboard TE point and going clockwise.
b1 = vor_length/2;

% Inboard and outboard vortex points
p1(:,:) = lattice.VORTEX(:,b1,:);
p2(:,:) = lattice.VORTEX(:,b1+1,:);

% Vortex Midpoint Location
lattice_midVortex(:,:) = (p1+p2)./2;

% Calculate the induced downwash from the trailing vortices, if requested
if isempty(IndDownwash)
    if IndDrag  == 1
        
        % Move the collocation points to the vortex midpoints
        lattice.COLLOC_orig = lattice.COLLOC;
        lattice.COLLOC      = lattice_midVortex;
        
        % Calculate the induced downwash due to drag
        [~, IndDownwash] = fastdw(lattice,[]);
        
        % Take symmetry into account here by reflecting the collocation and
        % 
        if (SIMXZ)
            lattice_dummy               = lattice;
            %lattice_dummy.COLLOC(:,2)   = -lattice_dummy.COLLOC(:,2);
            lattice_dummy.VORTEX(:,:,2) = - lattice_dummy.VORTEX(:,:,2);
            lattice_dummy.N(:,2)        = -lattice_dummy.N(:,2);
            [~, DW_sym]               = fastdw(lattice_dummy,[]);
            IndDownwash                 = IndDownwash + DW_sym;
        end
        
        DWX = IndDownwash(:,:,1);
        DWY = IndDownwash(:,:,2);
        DWZ = IndDownwash(:,:,3);
        lattice.COLLOC = lattice.COLLOC_orig;
    end
else
    DWX = IndDownwash(:,:,1);
    DWY = IndDownwash(:,:,2);
    DWZ = IndDownwash(:,:,3);
end

% Calculate the leading edge vortex vector
le = (p2-p1);
%le(:,1) = 0;

% Calculate the leading edge vortex length
Lle = sqrt(sum(le.^2,2));

% Normalise the leading edge vector
lehat(:,1) = le(:,1)./Lle;
lehat(:,2) = le(:,2)./Lle;
lehat(:,3) = le(:,3)./Lle;

% Calculate the induced downwash and compensate the inflow
if IndDrag ==1
    IW(:,1) = DWX*gamma;
    IW(:,2) = DWY*gamma;
    IW(:,3) = DWZ*gamma;
    
    Inflow = Inflow - IW(:,:);
end

% Calculate the direction and magnitude of the circulation
G(:,1) = gamma.*lehat(:,1);
G(:,2) = gamma.*lehat(:,2);
G(:,3) = gamma.*lehat(:,3);

% Calculate the Force per span length
Fprim = rho*cross(Inflow,G,2);   

% Calculate the force per panel
F(:,1) = Fprim(:,1).*Lle;				
F(:,2) = Fprim(:,2).*Lle;				
F(:,3) = Fprim(:,3).*Lle;				

% Store interesting quantities into the results structure
results.F               = F;
results.gamma           = gamma;
results.gammapanel      = gamma.*Lle;
results.AIC             = AIC;
results.DW              = DW;
results.DW_sym          = DW_sym;
results.IndDownwash     = IndDownwash;

end
