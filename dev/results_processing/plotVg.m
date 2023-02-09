function [hF, hF_Disp, hF_MC, FlutterModeData, FlutterMachVariation] = plotVg(f06File, varargin)
%plotVg Extracts the flutter data from a FO6 file and makes the V-g plot.
%   - TODO : Add mode tracking

hF                   = [];
hF_Disp              = [];
hF_MC                = [];
FlutterModeData      = [];
FlutterMachVariation = [];

%Pick a file
if nargin < 1 || isempty(f06File)
    [filename, dirOut] = uigetfile({'*.f06';}, 'Pick MSC.Nastran output file');
    if isnumeric(filename) || isnumeric(dirOut)
       return 
    end
    f06File = fullfile(dirOut, filename);
end

%Parse
p = inputParser;
addParameter(p, 'DampingLimit'       , 0.03 , @(x)validateattributes(abs(x), {'numeric'}, {'scalar', '<', 1}));
addParameter(p, 'VOperating'         , []   , @(x)validateattributes(abs(x), {'numeric'}, {'scalar', 'positive'}));
addParameter(p, 'VDive'              , []   , @(x)validateattributes(abs(x), {'numeric'}, {'scalar', 'positive'}));
addParameter(p, 'PlotUnstableModes'  , false, @(x)validateattributes(x     , {'logical'},  {'scalar'}));
addParameter(p, 'RetainUnstableModes', false, @(x)validateattributes(x     , {'logical'}, {'scalar'}));
addParameter(p, 'ModelFEM'           , []   , @(x)validateattributes(x     , {'awi.fe.FEModel'}, {'vector'}));
addParameter(p, 'PlotModeIndex'      , []   , @(x)validateattributes(x     , {'numeric'}, {'row', 'integer', 'positive'}));
addParameter(p, 'RetainAllModes'     , false, @(x)validateattributes(x     , {'logical'}, {'scalar'}));
parse(p, varargin{:});

%Define data for plotting annotations
zetaCrit = 100 * p.Results.DampingLimit;
VO       = p.Results.VOperating;
VD       = p.Results.VDive;

%% Extract flutter results

%Extract results
BaseFlutter = readF06FLUTTER(f06File);

%Any rigid body modes?
if p.Results.RetainAllModes
    idxRB = true(size(BaseFlutter.EigenSummary(2).Freq));
else
    idxRB = (BaseFlutter.EigenSummary(2).Freq ~= 0);
end
flexModeInd = find(idxRB == true, 1, 'first');
nFlexMode = nnz(idxRB);
switch BaseFlutter.Method
    case 'PK'
        nMode = size(BaseFlutter.Summary, 1);
    case 'PKNL'
        nMode = size(BaseFlutter.Summary, 2);
end

%Set up colours and labels for mode branches
str = arrayfun(@(iM) sprintf('Branch %i', iM), flexModeInd : nMode, 'Unif', false)';
clr = colorSet(nFlexMode, {'k', 'w'});
clr = num2cell(clr, 2);

%% Plot the data
switch BaseFlutter.Method
    
    case 'PK'
                
        nMach = numel(BaseFlutter.Mach);
        nVel  = numel(BaseFlutter.Summary(1).velocity);        
        zMin  = zeros(1, nMach);                
                
        %Set up hg objects
        hF    = arrayfun(@(i) figure('Name', ...
            sprintf('V-g plot - Mach = %.3f', BaseFlutter.Mach(i))), 1 : nMach);
        hAx   = arrayfun(@(i) axes('Parent', hF(ceil(i / 2)), 'NextPlot', 'add'), 1 : nMach * 2);
        hAx   = reshape(hAx, [2, nMach]);
        idx   = false(nFlexMode, nMach);
        
        %Plot V-g for each Mach
        for i = 1 : nMach
            
            %Grab velocity, freq & damping data
            FlutterData = BaseFlutter.Summary(:, i);
            vel  = vertcat(FlutterData.velocity)';
%             freq = vertcat(FlutterData.frequency)';
            freq = vertcat(FlutterData.kFreq)';
            damp = 100 .* vertcat(FlutterData.damping)'; 
            
            %Filter out rigid body modes
            vel  = vel(:, idxRB);
            freq = freq(:, idxRB);
            damp = damp(:, idxRB);
            
            %Plot frequency
            hL = plot(hAx(1, i), vel, freq, 'LineWidth', 1.5);
            if numel(hL) > 1
                set(hL, {'Color'}, clr);
                set(hL, {'DisplayName'}, str);
            end
            
            %Plot damping
            hL = plot(hAx(2, i), vel, damp, 'LineWidth', 1.5);
            if numel(hL) > 1
                set(hL, {'Color'}, clr);
                set(hL, {'DisplayName'}, str);
            end
            
            %Stash minimum damping value
            zMin(i) = min(min(damp));
            
            %Stash the unstable branches
            idx(:, i) =  any(damp > 0, 1);
            
            v_  = vel(:, idx(:, i));
            g_  = damp(:, idx(:, i));
            for ijk = 1 : nnz(idx(:, i))
                vF{i}(ijk) = interp1(g_(:, ijk), v_(:, ijk), 0);
            end
            
            %Consolidate complex eigenvalues
            re = vertcat(FlutterData.complex)';
            im = vertcat(FlutterData.eigenvalue)';
            
            %Calculate damping and frequency from complex eigenvalue
            %   - Calculated value of g matches 
            omega = im(idxRB);
            gamma = re(idxRB) ./ omega;
            g     = gamma .* 2;
            
        end
        zMin = min(zMin);
        
        %Retain any unstable branches across all plots
        if ~p.Results.RetainUnstableModes && numel(hL) > 1
            %Logical index of stable modes
            idx  = ~any(idx, 2);
            unstableInd = find(~idx) + flexModeInd - 1;
            idx_ = repmat(idx, [nMach, 1]);
            %Line handles
            hL  = flipud(findobj(hAx(1, :), 'Type', 'Line'));            
            set(hL(idx_), 'Visible', 'off');
            hL  = flipud(findobj(hAx(2, :), 'Type', 'Line'));            
            set(hL(idx_), 'Visible', 'off');
        end
        
    case 'PKNL'
        
        if numel(idxRB) ~= size(BaseFlutter.Summary, 2)
            bBadData = true;
            idxRB = true(size(BaseFlutter.Summary));
            str = arrayfun(@(iM) sprintf('Branch %i', iM), 1 : numel(idxRB), 'Unif', false)';
            clr = colorSet(numel(idxRB), {'k', 'w'});
            clr = num2cell(clr, 2);
        else
            bBadData = false;
        end
        
        % zeta = 100 .* vertcat(BaseFlutter.Summary.damping);
        vel  = vertcat(BaseFlutter.Summary.velocity)';
        
        %Consolidate damping ratios
        zeta0 = 100 .* vertcat(BaseFlutter.Summary.damping)';
        freq  = vertcat(BaseFlutter.Summary.kFreq)';
        
        %Get rid of rigid body-body
        freq   = freq(:, idxRB);
        zeta0  = zeta0(:, idxRB);
        vel    = vel(:, idxRB);
        
        %Only retain unstable modes
        idx = ~any(zeta0 >= 0, 1);

        %Make the plot
        hF = figure('Name', 'V-g plot - Matched Velocity/Density/Mach pairs');
        hAx(1, 1) = axes('Parent', hF, 'NextPlot', 'add');
        hAx(2, 1) = axes('Parent', hF, 'NextPlot', 'add');        
        
        %Plot frequency
        if size(vel, 1) == 1
            hL(:, 1) = arrayfun(@(i) plot(hAx(1, 1), vel(i), freq(i), 'Marker', 'o', 'LineStyle', 'none'), 1 : numel(vel));
        else
            hL(:, 1) = plot(hAx(1, 1), vel, freq, 'LineWidth', 1.5);
        end
        if numel(hL) > 1
            set(hL(:, 1), {'Color'}, clr);
            set(hL(:, 1), {'DisplayName'}, str);
        end
        
        %Plot damping
        if size(vel, 1) == 1
            hL(:, 2) = arrayfun(@(i) plot(hAx(2, 1), vel(i), zeta0(i), 'Marker', 'o', 'LineStyle', 'none'), 1 : numel(vel));
        else
            hL(:, 2) = plot(hAx(2, 1), vel, zeta0, 'LineWidth', 1.5);
        end
        if numel(hL) > 1
            set(hL(:, 2), {'Color'}, clr);
            set(hL(:, 2), {'DisplayName'}, str);
        end
        
        zMin = min(min(zeta0));
        
        %Switch off any stable branches
        if any(idx) && ~p.Results.RetainUnstableModes
            set(hL(idx, :), 'Visible', 'off');
            unstableInd = find(~idx) + flexModeInd - 1;
            if bBadData
                unstableInd = find(~idx);
            end
        end
        unstableInd = find(~idx) + flexModeInd - 1;
end

%% Format plots

%Set up subplot appearance
set(hAx(1, :), 'OuterPosition', [0, 0.5, 1, 0.5]);
set(hAx(2, :), 'OuterPosition', [0, 0  , 1, 0.5]);

%Take each axes so we know what is in it
set(hAx(1, :), 'Tag', 'V-k');
set(hAx(2, :), 'Tag', 'V-g');

%Format axes
set([hAx.XLabel], 'String', 'Velocity [m/s]');
set([hAx(2, :).YLabel], 'String', '$g \quad [\%]$');
set([hAx(1, :).YLabel], 'String', '$k \quad [-]$');
yl = get(hAx(2, :), {'YLim'});
yl = vertcat(yl{:});
yLim(1) = zMin;
yLim(2) = max(yl(:, 2));
set(hAx(2, :), 'YLim', yLim);

%Add annotations
arrayfun(@(h) plot(h, h.XLim, [0, 0], 'k-'), hAx(2, :));
arrayfun(@(h) plot(h, h.XLim, [zetaCrit, zetaCrit], 'r--'), hAx(2, :));
arrayfun(@(h) text(h, h.XLim(1), zetaCrit, sprintf('g = %.2g\\%%', zetaCrit), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left'), hAx(2, :));
if ~isempty(VO)
    arrayfun(@(h) plot(h, [VO, VO], h.YLim, 'k--'), hAx(2, :));
    arrayfun(@(h) text(h, VO, h.YLim(1), '$V_{C}$', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'), hAx(2, :));
end
if ~isempty(VD)
    arrayfun(@(h) plot(h, [VD, VD]        , h.YLim, 'k--'), hAx(2, :));
    arrayfun(@(h) plot(h, 1.15 .* [VD, VD], h.YLim, 'r--'), hAx(2, :));
    arrayfun(@(h) text(h, VD, h.YLim(1), '$V_{D}$', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'), hAx(2, :));
    arrayfun(@(h) text(h, 1.15 * VD, h.YLim(1), '$1.15V_{D}$', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right'), hAx(2, :));
end

%Force certain modes to be shown
index = p.Results.PlotModeIndex;
if ~isempty(index)
    %Account for rigid body modes
    index = index - flexModeInd + 1;
    %Update 'idx'
    idx        = true(nFlexMode, 1);
    idx(index) = false;
    for i = 1 : numel(hAx(2, :))
        %Grab mode branches
        hL = flipud(findobj(hAx(2, i), 'Type', 'Line'));
        hL = hL(arrayfun(@(h) ~isempty(h.DisplayName), hL));
        %Force visibility
        set(hL(~idx), 'Visible', 'on');
        set(hL(idx) , 'Visible', 'off');
    end
end

%Add a legend
for i = 1 : numel(hAx(1, :))
    %Grab mode branches
    hL = flipud(findobj(hAx(1, i), 'Type', 'Line'));
    hL = hL(arrayfun(@(h) ~isempty(h.DisplayName), hL));
    %Make the legend
    if ~isempty(hL)
        legend(hAx(1, i), hL(~idx(:, 1)), {hL(~idx(:, 1)).DisplayName}, 'Location', 'NorthEast');
    end
end

%Force the V-g curves to be on top
for i = 1 : numel(hAx(2, :))
    hAx(2, i).Children = flipud(hAx(2, i).Children);
end

%% Plot complex modeshapes

FEM = p.Results.ModelFEM;

%Get summary data for unstable modes
%   - 'searchStr' denotes the order in which the damping vector is
%   searched to find the first velocity where the damping goes unstable.
%   For the PK method the 
switch BaseFlutter.Method
    case 'PK'
        FSummary  = BaseFlutter.Summary(unstableInd, :);
    case 'PKNL'
        FSummary  = BaseFlutter.Summary(unstableInd)';
end

%Get complex eigenvectors from the h5 file
% ComplexEigData = BaseFlutter.ComplexEigenvector;
[path, name, ~] = fileparts(f06File);
h5File = fullfile(path, [name, '.h5']);
[h5Raw, ~, h5Results] = h5extract(h5File);

%Stash the complex eigenvalue - Use absolute value
if isfield(h5Results, 'EIGENVECTOR')
    re = abs([h5Results.EIGENVECTOR.TIME_FREQ_EIGR]);
    im = abs([h5Results.EIGENVECTOR.EIGI]);
else
    re = nan;
    im = nan;
end

nUSMode = numel(unstableInd);

%Preallocate
zrs    = zeros(1, nUSMode);
fNames = {'Branch' ; 'Mach' ; 'Velocity' ; 'DensityRatio' ; 'G' ; 'K' ; 'RealEig' ; 'ImagEig' ; 'ModeIndex'};
data   = [{unstableInd'} ; repmat({zrs}, [numel(fNames) - 1, 1])];
FlutterModeData = cell2struct(data, fNames);

%Find first point where each mode goes unstable - assume that is the
%flutter mode.
for iM = 1 : nUSMode
    %Grab damping data
    SummaryData = FSummary(iM, :);
    damp = vertcat(SummaryData.damping);
    vel  = vertcat(SummaryData.velocity);
    %First (e.g. lowest) Mach number and velocity to go unstable
    idxUnstable = damp > 0;
    machInd     = find(any(idxUnstable, 2), 1, 'first');
    unstableVel = min(vel(machInd, idxUnstable(machInd, :)));
    velInd      = find(vel(machInd, :) == unstableVel);
    if isempty(velInd) || isempty(machInd)
        FlutterModeData.ModeIndex(iM) = unstableInd(iM);
        continue
    end
    %Stash Mach, velocity, density and eigenvalue
    FlutterModeData.Mach(iM)         = SummaryData(machInd).mach(velInd);
    FlutterModeData.Velocity(iM)     = SummaryData(machInd).velocity(velInd);
    FlutterModeData.DensityRatio(iM) = SummaryData(machInd).density(velInd);
    FlutterModeData.G(iM)            = SummaryData(machInd).damping(velInd);
    FlutterModeData.K(iM)            = SummaryData(machInd).kFreq(velInd);
    FlutterModeData.RealEig(iM)      = SummaryData(machInd).complex(velInd);
    FlutterModeData.ImagEig(iM)      = SummaryData(machInd).eigenvalue(velInd);
    %Find the index number of the complex eigenvector that corresponds to
    %this flight point
    indRe = find(FlutterModeData.RealEig(iM) == re);
    if isempty(indRe)
        [~, indRe] = min(abs(abs(FlutterModeData.RealEig(iM)) - re));
    end
    indIm = find(FlutterModeData.ImagEig(iM) == im);
    if isempty(indIm)
        [~, indIm] = min(abs(abs(FlutterModeData.ImagEig(iM)) - im));
    end
    idx = indRe == indIm;
    if any(idx)
        if nnz(idx) > 1
            warning('Ambiguous match when matching flutter modeshapes. Setting mode index to NaN');
            FlutterModeData.ModeIndex(iM) =  nan;
        else
            FlutterModeData.ModeIndex(iM) = indIm(find(idx, 1, 'first'));
        end
    else
        FlutterModeData.ModeIndex(iM) =  nan;
    end
end

%Find variation in flutter speed as a function of Mach
nMach = size(FSummary, 2);
vF    = zeros(nUSMode, nMach);
mach = nan;
%Find flutter speed for each mach
for iM = 1 : nUSMode
    %Grab damping data
    SummaryData = FSummary(iM, :);
    damp = vertcat(SummaryData.damping);
    vel  = vertcat(SummaryData.velocity);
    mach = vertcat(SummaryData.mach);
    %First (e.g. lowest) velocity to go unstable for each Mach
    for iMach = 1 : nMach
        idxUnstable   = damp(iMach, :) > 0;
        indexUnstable = find(idxUnstable);
        [v, ind] = min(vel(iMach, idxUnstable));
        if indexUnstable(ind) == 1
            vF(iM, iMach) = v;
        else
            %Calcualte flutter velocity by interpolating
            vF(iM, iMach) = interp1(damp(iMach, :), vel(iMach, :), 0);
        end
    end
end
%Stash
FlutterMachVariation.Method          = BaseFlutter.Method;
FlutterMachVariation.Mach            = mach(:, 1)';
FlutterMachVariation.FlutterVelocity = vF;
FlutterMachVariation.UnstableModeInd = unstableInd;


%Request modeshape plots?
if ~p.Results.PlotUnstableModes || isempty(FEM) || ~isfield(h5Results, 'EIGENVECTOR')
    return
end

%Stash the modeshapes
FlutterModeData.Eigenvectors = h5Results.EIGENVECTOR(FlutterModeData.ModeIndex);

%   - Split this up into multiple functions
hF_Disp = plotDisplacementResult(FEM, h5File, 'PlotIndex', FlutterModeData.ModeIndex, 'DrawSurface', true, 'PlotScaleFactor', 100);
hAx = findobj(hF_Disp, 'Type', 'axes');

%Update the title and figure name
label = arrayfun(@(i) sprintf( ...
    'Flutter Branch %i, Velocity = %.3f m/s, Mach = %.3f, g = %.3f %%, k = %.3f', ...
    FlutterModeData.Branch(i)  , ...
    FlutterModeData.Velocity(i), ...
    FlutterModeData.Mach(i)    , ...
    FlutterModeData.G(i) * 100 , ...
    FlutterModeData.K(i)      ), 1 : nUSMode, 'Unif', false)'; 
set(hF_Disp, {'Name'}, label);
set([hAx.Title], {'String'}, label);

%Plot modal coords
label = strrep(label, 'Flutter', 'Flutter Modal Coords - ');
if isempty(BaseFlutter.ModalCoords)
    hF_MC = [];
else
    mc = BaseFlutter.ModalCoords.MPF(FlutterModeData.ModeIndex, :);
    hF_MC = arrayfun(@(i) figure('Name', label{i}), 1 : nUSMode);
    %Make a bar chart
    for i = 1 : nUSMode
        hAx_(1) = axes('Parent', hF_MC(i));
        bar(hAx_(1), real(mc(i, :)));
        hAx_(2) = axes('Parent', hF_MC(i));
        bar(hAx_(2), imag(mc(i, :)));
        set(hAx_, {'OuterPosition'}, {[0 0.5 1 0.5] ; [0 0 1 0.5]});
        set([hAx_.Title] , {'String'}, {'Real Parts' ; 'Imaginary Parts'});
        set([hAx_.YLabel], {'String'}, {'$\xi^{R} \quad [-]$' ; '$\xi^{I} \quad [-]$'});
        set([hAx_.XLabel], 'String'  , '$\textrm{Normal Mode Number} \quad [-]$');
    end
end
    
end

%% Old code
%   - Aimed at linking up every single complex eigenvector to every single
%   flutter summary point
%
% %Associate complex eigenvector results with flutter summary data
% ComplexEig = assignComplexEigenvectorData(BaseFlutter, unstableInd);
% 
% %Only retain unstable modes
% ComplexEig = ComplexEig(unstableInd, :);
% 
% %Split by modes then Mach number
% for iM = 1 : numel(unstableInd)
%     nam = ['Mode', num2str(unstableInd(iM))];
%     ComplexEigData.(nam) = reshape(ComplexEig(iM, :), [nMach, nVel]);
% end
% clear ComplexEig
% eig = complex([h5Results.EIGENVECTOR.TIME_FREQ_EIGR], [h5Results.EIGENVECTOR.EIGI]);
% eig_ = [ComplexEig.Eigenvalue];

function ComplexEig = assignComplexEigenvectorData(BaseFlutter, retentionModeIndex)
%assignComplexEigenvectorData Associates the complex eigenvector data with
%the flutter summary data.
%
% Because Nastran presents the Eigenvector data in a strange way we have no
% way of precisely matching the data presented in the "Flutter Summary" 
% with the data in the "Complex Eigenvectors"
% 
% Solution:
%   - Find the best match between the absolute values of the real and
%     imaginary parts of the complex eigenvectors presented in the Flutter
%     Summary and the Complex Eigenvector data.
%   - The complex eigenvalue provided in the Flutter summary often has the
%     opposite sign for the imaginary term as the eigenvector data. This is
%     because the eigenvalues extracted from the flutter analysis are
%     complex conjugate pairs. As a workaround we compare the absolute
%     values of the real and imaginary terms.

ModeData   = BaseFlutter.Summary;
nModes     = size(ModeData.Summary, 1);
ComplexEig = [];

if nargin < 2
    retentionModeIndex = 1 : nModes;
end

%Velocity/Mach/Density data
mach         = vertcat(ModeData.mach);
vel          = vertcat(ModeData.velocity);
densityRatio = vertcat(ModeData.density);
damp         = vertcat(ModeData.damping);
mach         = mach(:);
vel          = vel(:);
densityRatio = densityRatio(:);
damp         = damp(:);

%Grab complex eigenvalue from flutter summary
reSum   = abs(vertcat(ModeData.complex));
imSum   = abs(vertcat(ModeData.eigenvalue));
reSum   = reSum(:);
imSum   = imSum(:);

%Complex eigenvalue for all complex eigenvectors
allComplexEig = [BaseFlutter.ComplexEigenvector.Eigenvalue];
re = abs(real(allComplexEig));
im = abs(imag(allComplexEig));

%Calculate the difference in real and imaginary terms
delRe = abs(reSum - re);
delIm = abs(imSum - im);
[delRe, indRe] = min(delRe);
[delIm, indIm] = min(delIm);

%Only proceed if the real and imaginary deltas agree
idx = (indRe - indIm) == 0;

if ~any(idx)
    return
end

%Make the Complex Eigenvector data structure
fNames = fieldnames(BaseFlutter.ComplexEigenvector(1));
fNames = [fNames ; {'G' ; 'Mach' ; 'Vel' ; 'DensityRatio' ; 'DelRe' ; 'DelIm'}];
ComplexEig = repmat(cell2struct(cell(size(fNames)), fNames), [1, numel(damp)]);

%Grab complex eigenvector dat
cEig = squeeze(struct2cell(BaseFlutter.ComplexEigenvector));

%Convert to cell notation
damp  = num2cell(damp);
delRe = num2cell(delRe);
delIm = num2cell(delIm);
mach  = num2cell(mach);
vel   = num2cell(vel);
dens  = num2cell(densityRatio);

%Assign data - Only assign data where the re and im terms have matched
[ComplexEig(idx).Eigenvalue] = cEig{1, idx};
[ComplexEig(idx).GID]        = cEig{2, idx};
[ComplexEig(idx).T1]         = cEig{3, idx};
[ComplexEig(idx).T2]         = cEig{4, idx};
[ComplexEig(idx).T3]         = cEig{5, idx};
[ComplexEig(idx).R1]         = cEig{6, idx};
[ComplexEig(idx).R2]         = cEig{7, idx};
[ComplexEig(idx).R3]         = cEig{8, idx};
[ComplexEig(idx).G]          = damp{idx};
[ComplexEig.Mach]            = mach{:};
[ComplexEig.Vel]             = vel{:};
[ComplexEig.DensityRatio]    = dens{:};
[ComplexEig.DelRe]           = delRe{:};
[ComplexEig.DelIm]           = delIm{:};

%Reshape
nPoints    = numel(ComplexEig) / nModes;
ComplexEig = reshape(ComplexEig, [nModes, nPoints]);

end

