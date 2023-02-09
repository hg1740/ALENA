function x_out = getTrim(x0,Matrices,Aero,Sim,TrimVars,hp,maxTrimIter)

if nargin < 6
    hp = [];
end
if nargin < 7 || isempty(maxTrimIter)
    maxTrimIter = 50;
end

%% Set up window for analysis and other graphics objects

nPart     = numel(Matrices.n_elem);
legend_in = [TrimVars.Flaps(:, 4) ; TrimVars.rb_vars(:, 3)];

%Make the axes object for drawing trim quantities
ha = setupTrimGraphics(hp, legend_in, nPart);

%Make a waitbar so the user knows something is happening.
hwb = waitbar(0, 'Trimming...please wait, or hit ''Cancel'' to interrupt', ...
    'CreateCancelBtn', 'delete(gcbf)', ...
    'Name', 'AWI Framework');
clu = onCleanup(@()delete(hwb(ishghandle(hwb))));

%% Trim parameters

%TODO - Add these values to the 'TrimOptions' object and pass in to this
%function.
if ~isfield(TrimVars,'NetForce')
    TrimVars.NetForce = [0;0;0;0;0;0];
end
if ~isfield(TrimVars,'Damping')
    TrimVars.Damping = 0.5;
end
if ~isfield(TrimVars,'eps')
    TrimVars.eps = 1e-5;
end

%Define analysis parameters for the trim solution
TrimParam = awi.methods.TrimOptions;
TrimParam.TrimVariableStepSize = 1e-3;
TrimParam.MaxIterations        = maxTrimIter;
TrimParam.NumericalDamping     = TrimVars.Damping;
TrimParam.ResidualTolerance    = TrimVars.eps;
TrimParam.ForceBalanceIndex    = TrimVars.ForceBalance;

%Define analysis parameters for static solution
AnalysisParam = awi.methods.Options;
%   - Custom parameters for static analysis solution
AnalysisParam.AnalysisType     = 'static';
AnalysisParam.NumericalDamping = 0.5;
AnalysisParam.MaxIterations    = 1e3;
%   - Define old parameters for the simulation
Sim           = initSystem(x0, Matrices, Aero, AnalysisParam);
Sim.Soln      = 1;
Sim.aero_flag = 1;
Sim.rb_flag   = 1;

Aero.inflow_switch = 1;

%% Initial conditions

%Construct vector of trim variables (control surfaces & rigid body motions) 
J = vertcat(TrimVars.Flaps{:, 3}, TrimVars.rb_vars{:, 2});

nControlSurf  = size(TrimVars.Flaps, 1);
nRigidBodyDOF = size(TrimVars.rb_vars, 1);

%Define global-to-aircraft rotaton matrix 'CGa'
rb_eul            = zeros(3,1);
trimIndex         = vertcat(TrimVars.rb_vars{:, 1});
rb_eul(trimIndex) = vertcat(TrimVars.rb_vars{:, 2});
Matrices.CGa      = Euler2Rot(rb_eul);

%Get inital coordinates & orientations of nodes
x_trim       = initSim(x0, Matrices, AnalysisParam);
[coords0, ~] = strains2coords_all(x_trim.x_f, Matrices);

%% Run Newton-Raphson trim solution

%Trim parameters
%   - Step size for incrementing the trim variables 
ds            = TrimParam.TrimVariableStepSize; 
%   - Invoke dependent parameters once!
extract_loads = TrimParam.TrimLoadsIndexMatrix;

%Trim the model
iTrim = 0;
while and(iTrim < TrimParam.MaxIterations, ishghandle(hwb))
    
    %Update the waitbar
    waitbar(iTrim / maxTrimIter, hwb);
    iTrim = iTrim + 1;    

    %% Initial static solution
    
    %Calculate the deformed shape & associated loads for the current states
    x_stat      = getStatic(x_trim, Matrices, Aero, Sim, AnalysisParam);
    globalLoads = x_stat.total_Aero_G + x_stat.total_Grav_G + TrimVars.NetForce;

    %Total trim residual and residual in each trim DOF
    res0         = extract_loads * globalLoads;
    trimResidual = norm(res0);
        
    %Update the trim graphics so the user knows the solver is working
    updateTrimGraphics(ha, iTrim, trimResidual, J, Matrices, x_stat)
    
    %Convergence reached?
    if trimResidual < TrimParam.ResidualTolerance
        break
    end
    
    %% Calculate numerical jacobian by perturbing each trim state
    dresdJ = zeros(numel(res0), numel(J));
    var_counter = 0;
    
    %Control surfaces
    for ii = 1 : nControlSurf
        var_counter = var_counter + 1;
        
        %Update the control surface deflection 
        %   - TODO : Generalise
        for jj = 1:length(TrimVars.Flaps{ii,1})
            Aero.delta_flap{TrimVars.Flaps{ii,1}(jj)}(TrimVars.Flaps{ii,2}{jj}) = J(var_counter) + ds;
        end
        
        %Run new static solution using updated control surface deflection 
        x_trim      = initSim(x0, Matrices, AnalysisParam); 
        x_stat      = getStatic(x_trim, Matrices, Aero, Sim, AnalysisParam);
        globalLoads = x_stat.total_Aero_G + x_stat.total_Grav_G + TrimVars.NetForce;
        res         = extract_loads * globalLoads;
        
        %Calculate rate-of-change in residual w.r.t this control surface
        dresdJ(:, var_counter) = (res - res0) / ds;
        
        %Restore to previous value
        for jj = 1:length(TrimVars.Flaps{ii,1})
            Aero.delta_flap{TrimVars.Flaps{ii,1}(jj)}(TrimVars.Flaps{ii,2}{jj}) = J(var_counter);
        end
        
    end
    
    %Rigid body trim DOFs
    for ii = 1 : nRigidBodyDOF
        var_counter = var_counter + 1;
        
        %Update the global-to-aircraft frame
        ds_vec                = zeros(3,1);
        ds_vec(trimIndex(ii)) = ds;
        Matrices.CGa          = Euler2Rot(rb_eul+ds_vec);
        
        %Run new static solution using updated states
        x_trim      = initSim(x0, Matrices, AnalysisParam); %Update velocities in aircraft frame
        x_stat      = getStatic(x_trim, Matrices, Aero, Sim, AnalysisParam);
        globalLoads = x_stat.total_Aero_G + x_stat.total_Grav_G + TrimVars.NetForce;
        res         = extract_loads * globalLoads;
        
        %Calculate rate-of-change in residual w.r.t this rigid DOF
        dresdJ(:, var_counter) = (res - res0)/ds;
        
        %Reset global orienation to current trim step value (unpeturbed)
        Matrices.CGa          = Euler2Rot(rb_eul);
        
    end    
        
    %% Update trim variables based on numerical jacobian
    
    %Calculate the new trim variables
    %   J_(n + 1) = J_(n) + zeta * dJ/dR * R
    J = J - (TrimParam.NumericalDamping * pinv(dresdJ) * res0);
    
    %Control surfaces
    var_counter = 0;
    for ii = 1:size(TrimVars.Flaps,1)
        var_counter = var_counter + 1;
        for jj = 1:length(TrimVars.Flaps{ii,1})
            Aero.delta_flap{TrimVars.Flaps{ii,1}(jj)}(TrimVars.Flaps{ii,2}{jj}) = J(var_counter);
        end
    end
    
    %Update the global orientation with the current rigid body states
    rb_eul = zeros(3,1);
    rb_eul(trimIndex) = J(nControlSurf + 1 : end);
    Matrices.CGa = Euler2Rot(rb_eul);    
    
    %Calculate new state vector based on updated trim variables
    x_trim = initSim(x0, Matrices, AnalysisParam);
    x_stat = getStatic(x_trim, Matrices, Aero, Sim, AnalysisParam);
    
    %Update the reference states so the jacobian is calculated correctly
    x0.x_f       = x_stat.x_f;
    
end

%% Postprocessing aero loads & displacements

% Determine the cl distribution
dyn_pres = 0.5 * Aero.rho * norm(x0.x_vg(1:3))^2;

% TODO - Matrices.s is actually the length along the beam, so the trapezoid
% rule doesn't quite work here
%   - Needs to be width of the strip
AeroLoads        = arrayfun(@(nE) zeros(nE, 3) , Matrices.n_elem, 'Unif', false);
AeroAreas        = arrayfun(@(nE) zeros(1 , nE), Matrices.n_elem, 'Unif', false);
AeroNDLoads      = arrayfun(@(nE) zeros(nE, 3) , Matrices.n_elem, 'Unif', false);
AeroNDLoadsChord = arrayfun(@(nE) zeros(nE, 3) , Matrices.n_elem, 'Unif', false);
for jj = 1 : nPart 
    AeroLoads{jj}   = [x_stat.aeroForce_G{jj}(1:6:end),x_stat.aeroForce_G{jj}(2:6:end),x_stat.aeroForce_G{jj}(3:6:end)].*diff(Matrices.s{jj})';
    AeroAreas{jj}   = Aero.chord{jj}.*diff(Matrices.s{jj});
    AeroNDLoads{jj} =AeroLoads{jj}./(dyn_pres*AeroAreas{jj}');
    AeroNDLoadsChord{jj} = AeroNDLoads{jj}.*(Aero.chord{jj}'./(Aero.chord{jj}(1)));
end

%Calculate the local displacements
[coords, ~]        = strains2coords_all(x_stat.x_f, Matrices);
LocalDisplacements = arrayfun(@(i) coords{i} - coords0{i}, 1 : nPart, 'Unif', false);

%% Stash outputs

x_out.CGa        = Matrices.CGa;
x_out.x_vg       = x0.x_vg;
x_out.x_f        = x_stat.x_f;
x_out.delta_flap = Aero.delta_flap;
x_out.x_disp     = LocalDisplacements;
x_out.Cl         = AeroNDLoads;
x_out.Cll        = AeroNDLoadsChord;

end

%Plotting trim quantities
function ha = setupTrimGraphics(hp, legend_in, nPart)
%setupTrimGraphics Makes the axes objects for plotting trim quantities
%during the trim solution. 
%
% Detailed Description:
%   - Adds dummy hg-objects to axes to improve downstrea graphics updates
%   - Autopopualtes a legend item for the trim variable plot

%Make the axes objects
if isempty(hp)
    h_res  = figure('Name','Residual, Trim Geometry and Control Effectors');
    ha     = subplot(2,2,[1 3],'Parent',h_res);
    ha(2)  = subplot(2,2,2,'Parent',h_res);
    ha(3)  = subplot(2,2,4,'Parent',h_res);        
elseif isa(hp,'uiextras.VBox')
    hv    = uiextras.HBox('Parent',hp);
    ha(1) = axes('Parent', uicontainer('Parent', hp));    
    ha(2) = axes('Parent', uicontainer('Parent', hv));
    ha(3) = axes('Parent', uicontainer('Parent', hv));
    %Adjust heights
    hp.Heights = [-2, -1];
else
    error('Bad input!')
end

%Update axes appearance
set(ha,  ...
    'XGrid'     , 'on', 'YGrid'     , 'on', 'ZGrid'     , 'on', ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on', 'ZMinorGrid', 'on');
%   - Deformed shape
xlabel(ha(1), 'X [m]');
ylabel(ha(1), 'Y [m]');
zlabel(ha(1), 'Z [m]');
set(ha(1)   , 'View', [-90, 0]);
%   - Trim Residual
xlabel(ha(2), 'Iteration')
ylabel(ha(2), 'Residual Norm')
set(ha(2)   , 'YScale', 'log');
%   - Trim DOFs
legend(ha(3), 'Location','best');
xlabel(ha(3), 'Iteration')
ylabel(ha(3), 'Control Variable Values')

%Fix the formatting (apply hold to each axes)
arrayfun(@(hax) hold(hax, 'on'), ha);

%Add dummy hg-objects
arrayfun(@(x) plot(ha(1), nan, nan, 'LineStyle', '-', ...
    'Marker', 'none', 'Color', 'k'), 1 : nPart);
plot(ha(2), nan, nan, 'LineStyle', '-', 'Color', 'k', 'Marker', '*');
hTV = cellfun(@(x) plot(ha(3), nan, nan, 'DisplayName', x, ...
    'LineStyle', '-', 'Marker', 'o', 'MarkerEdgeColor', 'k'), legend_in);
set(hTV, {'MarkerFaceColor'}, get(hTV, {'Color'}));

end

function updateTrimGraphics(ha, trimStepNumber, residual, trimStates, Matrices, x_stat)
%updateTrimGraphics Updates the trim graphics with the current states.
%
% Detailed Description:
%   - Plots the trim residual
%   - Plots the trim state values (aero controllers & rigid body terms)
%   - Plots the deformed shape


%Convert trim states to degrees 
%   * Vector is flipped so we can map the values across to the hg objects
%   which are in reverse order.
%   * TODO - assumes all trim states are angles!!
trimStates = flipud(trimStates * 180 / pi);

%Update the trim residual
hRes = ha(2).Children;
hRes.XData = [hRes.XData, trimStepNumber];
hRes.YData = [hRes.YData, residual];

%Update the trim DOF values
hTrimVar = ha(3).Children;
xd = [hTrimVar(1).XData, trimStepNumber];
yd = get(hTrimVar, {'YData'});
yd = arrayfun(@(i) [yd{i}, trimStates(i)], 1 : numel(trimStates), 'Unif', false); 
set(hTrimVar, 'XData', xd, {'YData'}, yd');

%Get current coordinates & orientations (aircraft frame)
[coords, CaB] = strains2coords_all(x_stat.x_f, Matrices);

%Rotate coordinates into gobal frame
coords_rot = i_aircraftToGlobal(coords, CaB, Matrices);

%Update the deformed shape
hBeam = ha(1).Children;
xd = cellfun(@(x) x(1, :), coords_rot, 'Unif', false)';
yd = cellfun(@(x) x(2, :), coords_rot, 'Unif', false)';
zd = cellfun(@(x) x(3, :), coords_rot, 'Unif', false)';
set(hBeam, {'XData'}, xd, {'YData'}, yd, {'ZData'}, zd);

%Propogate changes
% drawnow

    function coords_rot = i_aircraftToGlobal(coords, CaB, Matrices)
        %i_aircraftToGlobal Transforms the coordinates 'coords' into the
        %global frame from the aircraft frame.
        
        nPart = numel(Matrices.n_elem);
        coords_rot = cell(1, nPart);
        
        for jj = 1 : nPart
            
            x_(:,1) = coords{jj}(1,:,:)+Matrices.p0{jj}(1);
            y_(:,1) = coords{jj}(2,:,:)+Matrices.p0{jj}(2);
            z_(:,1) = coords{jj}(3,:,:)+Matrices.p0{jj}(3);
            
            if Matrices.Parent{jj}(1)
                coords_parent  = coords{Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)));
                CaB_parent     = CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)};
                coords_rot{jj} = Matrices.CGa*(repmat(coords_parent,1,length(x_)) + CaB_parent*[x_';y_';z_']);
            else
                coords_rot{jj} = Matrices.CGa*[x_';y_';z_'];
            end
            
            clear x_ y_ z_
            
        end
    end

end
