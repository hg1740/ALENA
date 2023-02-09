function [x_out, bSuccess] = getStatic(x0, Matrices, Aero, Sim, AnalysisParam, hp)
%getStatic Solves the static equations of motion of the system described by
%'Matrices', 'Aero' & 'Sim' using the parameters in 'AnalysisParam'.
%
% Detailed Information:
%   - Plots diagnostic information in the graphics container 'hp'.
%

bSuccess = false;

%Grab variables from analysis object 
maxIter  = AnalysisParam.MaxIterations;
epsilon  = AnalysisParam.ResidualTolerance;
zeta     = AnalysisParam.NumericalDamping;
nTangent = AnalysisParam.TangentMatrixUpdateFrequency; 

%Set up iteration counters
iter   = 1 : maxIter;
r_norm = 1e9 .* ones(1, maxIter);

%% Set up graphics objects for plotting convergence of system residual
bClose = false;
%   - Container
if nargin < 6 
    hp     = figure('Name', 'Static solution');
    bClose = true;
end
%   - Axes
hAx = axes('Parent', hp, 'NextPlot', 'add', 'XLim', [0 10], 'YScale', 'log', 'Box', 'on');
xlabel(hAx, 'Iteration [-]');
ylabel(hAx, 'Residual Norm [-]');
%   - Line
hRes = line(hAx, 'XData', nan, 'YData', nan, ...
    'Marker', 's', ...
    'Color' , 'b', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'b', ...
    'LineStyle'      , '-');
hEps = line(hAx, 'XData', hAx.XLim, 'YData', [epsilon, epsilon], ...
    'LineStyle', '-', ...
    'LineWidth', 2  , ...
    'Color'    , 'r', ...
    'Marker'   , 'none');

hAx(2) = axes('Parent', hp, 'NextPlot', 'add', ...
    'XLim'    , [-20, 30], ...
    'YLim'    , [0, 85]  , ...
    'ZLim'    , [-20, 60]  , ...
    'View'    , [150, 35], ...
    'Box'     , 'on'     , ...
    'XGrid'   , 'on'     , ...
    'YGrid'   , 'on'     , ...
    'ZGrid'   , 'on');
set(hAx, {'OuterPosition'}, {[0 0 1 0.5] ; [0 0.5 1 0.5]});
xlabel(hAx(2), 'X [m]');
ylabel(hAx(2), 'Y [m]');
zlabel(hAx(2), 'Z [m]');
[coords, ~] = strains2coords_all(x0.x_f, Matrices);
hNode = plot3(hAx(2), squeeze(coords{1}(1, :, :)), squeeze(coords{1}(2, :, :)), squeeze(coords{1}(3, :, :)), 'Marker', 'o', 'Color', 'k', 'LineStyle', '-', 'MarkerFaceColor', 'g');

%% Initialise the System and State
X    = initState(x0, Matrices, Aero, Sim, AnalysisParam);   %Initial state vector
X_p1 = stepState(X , Matrices, Aero, Sim, AnalysisParam);   %First estimate of the state vector

%% Loop until convergence is reached
ii = 1;
while ii <= maxIter 
    
    %Flag to determine whether to update tangent matrix this iteration
    bUpdateTangent = (mod((ii - 1), nTangent) == 0);
    
    % Get all the system matrices, the tangent matrix and the residual
    SYS               = getSystemMatrices(X_p1, Matrices, Aero, Sim, AnalysisParam);
    [r_full, Q, S, D] = getResidual(X_p1, Matrices, SYS, Sim, AnalysisParam, bUpdateTangent);
    
    %If the tangent matrix been updated then replace the current value with
    %the new value
    if bUpdateTangent
        S_sparse = S;
        D_scale  = D;
    end        
    
    %Update the state vector using Newton-Raphson
    dq            = - (D_scale * S_sparse)\(D_scale * r_full);
    Q.q_p1_full   = Q.q_p1_full + zeta * dq;
    
    %Reallocate the full state vector to the cell structure
    for jj = 1:length(Matrices.n_elem)
        X_p1.x_f_p1{jj} = Q.q_p1_full(Sim.Ind.x_f_ind{jj});
    end
    
    %Residual norm
    r_norm(ii) = abs(norm(r_full));
    
    %Update graphics
    set(hRes, 'XData', iter(1 : ii), 'YData', r_norm(1 : ii));    
    hAx(1).XLim = [0, ii];
    set(hEps, 'XData', get(hAx(1), 'XLim'));
    [coords, ~] = strains2coords_all(X_p1.x_f_p1, Matrices);
    set(hNode, 'XData', squeeze(coords{1}(1, :, :)), 'YData', squeeze(coords{1}(2, :, :)), 'ZData', squeeze(coords{1}(3, :, :)));
    drawnow
    
    %Was the run successful?
    if r_norm(ii) < epsilon
        bSuccess = true;
        break
    end

    %Update the iteration counter
    ii = ii + 1;
    
end

%Final states
x_out = updateResults(X, X_p1, SYS, Matrices, r_full, Sim, AnalysisParam);

if bClose %Close graphics window?
    close(hp);
end

end
