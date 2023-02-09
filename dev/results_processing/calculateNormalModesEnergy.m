function hF = calculateNormalModesEnergy(FEM, h5File, varargin)
%calculateNormalModesEnergy

%Parse
if nargin < 1
    return
else
    assert(isa(FEM, 'awi.fe.FEModel'), ['Expected the first input ', ...
        'argument to be a valid ''awi.fe.FEModel'' object.']);
end
p= inputParser;
addParameter(p, 'PlotIndex', []   , @(x)validateattributes(x, {'numeric'}, {'row', 'integer', 'positive'}));
addParameter(p, 'PlotMAC'  , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
parse(p, varargin{:});

%Get the data
[h5Raw, ~, h5Results] = h5extract(h5File);
% if ~isfield(h5Results, 'SUMMARY') && ~isfield(h5Results.SUMMARY, 'EIGENVALUE')
%     error('Expected the .h5 file to reference a SOL103 results file.');    
% end
% assert(all(h5Results.SUMMARY.EIGENVALUE.MASS == 1), ['Expected the ', ...
%     '.h5 file to reference a SOL103 results file which calculated ' , ...
%     'mass-normalised modes.']);
% if ~isfield(h5Results, 'KINETIC_ENERGY') || ~isfield(h5Results.ENERGY, 'STRAIN_ELEM')
%     error(['Expected the .h5 file to reference a SOL 103 results file '  , ...
%         'containing grid point kinetic energy (GPKE) and element strain ', ...
%         'energy (ESE)']);
% end

%Flatlist the collection
AllFEM = flatlist(FEM);

%Rotation matrix data
[rMatrix, cid] = getCoordSysRMatrix(h5Raw.NASTRAN.INPUT);

%Grab energy data
[KE, SE, AllFEM] = extractModeEnergyBreakdown(h5Results, AllFEM, rMatrix, cid);
nPart  = numel(AllFEM);

%Determine total kinetic energy from all elements
if isfield(h5Results, 'SUMMARY')
    EigSummary = h5Results.SUMMARY.EIGENVALUE;
    seTot = EigSummary.OMEGA .^ 2 ./ 2;
else
   %Try the f06?
   [path, name, ~] = fileparts(h5File);
   %Look for an .f06 in the same directory as the .h5
   f06File = fullfile(path, [name, '.f06']);
   Results = readF06FLUTTER(f06File);
   if isfield(Results, 'EigenSummary')
       data = Results.EigenSummary(1);
       assert(all(data.GenMass == 1), ['Expected the eigenanalysis to ', ...
           'use mass-normalised modes.']);
       seTot = data.Omega.^ 2 ./ 2;
   else
       warning(['Unable to compute total Strain Energy as no ', ...
           'Eigenanlysis summary data could be found in the h5 or f06 file.']);
       return
   end
end
nMode = numel(seTot);

%Determine total percentage of SE & KE across all parts per mode
seByPart = cellfun(@(x) sum(x, 1), SE, 'Unif', false);
seByPart = vertcat(seByPart{:});
keByPart = cellfun(@(x) sum(x, 1), KE, 'Unif', false);
keByPart = vertcat(keByPart{:});

%Normalise the strain energy (KE already normalised)
if numel(seTot) ~= size(seByPart, 2)
    seByPart_norm = nan;
    nMode = numel(h5Results.EIGENVECTOR);
else
    seByPart_norm = 100 .* seByPart ./ repmat(seTot', [nPart, 1]);
end

%Fracton of KE and SE of flexible nodes as a fraction of model total
partSETot  = sum(seByPart, 1);
partKETot  = squeeze(sum(keByPart, 1));
partSEFrac =  partSETot' ./ seTot;

%% Calculate the MAC

if p.Results.PlotMAC
    
    nMode = numel(h5Results.EIGENVECTOR);
    
    %Construct modeshape matrix (PHI)
    PHI = [ ...
        [h5Results.EIGENVECTOR.X]  ; ...
        [h5Results.EIGENVECTOR.Y]  ; ...
        [h5Results.EIGENVECTOR.Z]  ; ...
        [h5Results.EIGENVECTOR.RX] ; ...
        [h5Results.EIGENVECTOR.RY] ; ...
        [h5Results.EIGENVECTOR.RZ] ];
    
    
    %Construct ESE & GPKE matrix
    ESE  = vertcat(SE{:});
    GPKE = cell(size(KE));
    for i = 1 : numel(KE)
        ke = arrayfun(@(j) KE{i}(:, :, j), 1 : 6, 'Unif', false);
        GPKE{i} = vertcat(ke{:});
    end
    GPKE = vertcat(GPKE{:});
    
    %Calculate MAC matrices
    MAC      = zeros(nMode);
    ESE_MAC  = zeros(nMode);
    GPKE_MAC = zeros(nMode);
    for ii = 1 : nMode
        for jj = 1 : nMode
            MAC(ii,jj)      = i_calculateMAC(PHI(:,ii) , PHI(:,jj));
            ESE_MAC(ii,jj)  = i_calculateMAC(ESE(:,ii) , ESE(:,jj));
            GPKE_MAC(ii,jj) = i_calculateMAC(GPKE(:,ii), GPKE(:,jj));
        end
    end        
    
    %Consolidate MAC quantities into a single measure
    CombinedMAC = sqrt(MAC.^2 + ESE_MAC.^2 + GPKE_MAC.^2);
    AverageMAC  = (MAC + ESE_MAC + GPKE_MAC) ./ 3;
    
end

    function mac = i_calculateMAC(phi1, phi2)
        %i_calculateMAC Calculates the MAC between two eigenvectors.
        mac= (abs(phi1'*phi2))^2/((phi1'*phi1)*(phi2'*phi2));
    end

%% Plotting

%Get plot options
%   - Plot Index
pInd = p.Results.PlotIndex;
if isempty(pInd)
    %Plot everything
    pInd = 1 : nMode;
else
    %Check the idnexing is within bounds
    expInd = (1 : nMode)';
    pInd = pInd(and(any(pInd >= expInd), any(pInd <=  expInd)));
    if isempty(pInd)
        warning('matlab:ALENA:badResultsIndex', ['Requested results '  , ...
            'indices are outside the bounds of the results data. This ', ...
            'can happen when the user has requested a results index '  , ...
            'greater or less than the total number of results.'        , ...
            '\n\n\tUSER ACTION: Define new values of ''PlotIndex''.\n\n']);
        return
    end
end

% clr = num2cell(colorSet(nModes), 2);
% 
% %Distribution of strain & kinetic energy
% for i = 1 : nPart
%     
%     %Make the graphics object
%     figure('Name', sprintf('%s - Energy distribution', AllFEM(i).Name));
%     
%     %Eta positions
%     [~, beamR] = calculateBeamEta(AllFEM(i));
%     eta = beamR ./ beamR(end);
%     etaE = eta(1 : end - 1) + diff(eta) ./ 2;
%     
%     %Strain energy
%     hAx(1) = subplot(2, 1, 1);    
%     hL = plot(hAx(1), etaE, SE{i}');
%     set(hL, {'Color'}, clr);
%     title(hAx(1), 'Strain Energy');
%     
%     %Kinetic energy
%     hAx(2) = subplot(2, 1, 2);
%     hL = plot(hAx(2), eta, sum(KE{i}, 3)');
%     set(hL, {'Color'}, clr);
%     title(hAx(2), 'Kinetic Energy');
%     
% end

%% Plot the bar charts for each mode that has been requested
nRes = numel(pInd);
hF_ = arrayfun(@(~) figure, 1 : nRes * 2);
hF_ = reshape(hF_, [numel(hF_) / 2, 2]);
str = [ ...
    arrayfun(@(i) sprintf('Mode %i - KE Summary' , pInd(i)), 1 : nRes, 'Unif', false), ...
    arrayfun(@(i) sprintf('Mode %i - SE Summary' , pInd(i)), 1 : nRes, 'Unif', false)] ;
set(hF_, {'Name'}, str');
hAx = arrayfun(@(hf) axes('Parent', hf, 'NextPlot', 'add'), hF_);
for i = 1 : nRes 
    %KE
    bar(hAx(i, 1), squeeze(keByPart(:, pInd(i), :)));
    %SE
    if isnan(seByPart_norm)
        bar(hAx(i, 2), squeeze(seByPart(:, pInd(i))));
    else
        bar(hAx(i, 2), squeeze(seByPart_norm(:, pInd(i))));
    end
end
set([hAx.XAxis], 'TickValues', 1 : nPart);
set([hAx.XAxis], 'TickLabels', {'SW' ; 'PW' ; 'SS' ; 'PS' ; 'SJ' ; 'PJ'});

%% Plot the strain & KE energy by part
hF(1)  = figure('Name', 'Strain Energy by part');
hAx(1) = axes('Parent', hF(1));
if isnan(seByPart_norm)
    b = bar3(seByPart);    
else
    b = bar3(seByPart_norm);
end
set(b, {'CData'}, get(b, {'ZData'}));
set(b, 'FaceColor', 'interp');

%Plot the strain energy by part
hF(2)  = figure('Name', 'Kinetic Energy by part');
hAx(2) = axes('Parent', hF(2));
b = bar3(sum(keByPart, 3));
set(b, {'CData'}, get(b, {'ZData'}));
set(b, 'FaceColor', 'interp');

%Format axes
set([hAx.XLabel], 'String'  , 'Mode Number [-]');
set([hAx.ZLabel], {'String'}, {'SE [\%]' ; 'KE [\%]'});
set([hAx.YAxis] , 'TickLabels', {'SW' ; 'PW' ; 'SS' ; 'PS' ; 'SJ' ; 'PJ'});

%% Plot MAC quantities

if p.Results.PlotMAC
    
    % plot mac matrix
    hF(end + 1) = figure('Name', 'MAC');
    hAx = axes('Parent', hF(end));
    b = bar3(hAx, MAC);
    colorbar
    title('MAC')
    set(b, {'CData'}, get(b, {'ZData'}));
    view([-90 90]);
    %
    hF(end + 1) = figure('Name', 'MAC_ESE');
    hAx = axes('Parent', hF(end));
    b = bar3(hAx, ESE_MAC);
    title('ESE MAC')
    set(b, {'CData'}, get(b, {'ZData'}));
    view([-90 90]);
    %
    hF(end + 1) = figure('Name', 'MAC_GPKE');
    hAx = axes('Parent', hF(end));
    b = bar3(hAx, GPKE_MAC);
    title('GPKE MAC')
    set(b, {'CData'}, get(b, {'ZData'}));
    view([-90 90]);
    
    hF(end + 1) = figure('Name', 'MAC_RSS');
    hAx = axes('Parent', hF(end));
    b = bar3(hAx, CombinedMAC);
    title('RSS(MAC, ESE_MAC, GPKE_MAC)')
    set(b, {'CData'}, get(b, {'ZData'}));
    view([-90 90]);
    %
    hF(end + 1) = figure('Name', 'MAC_MEAN');
    hAx = axes('Parent', hF(end));
    b = bar3(hAx, AverageMAC);
    title('MEAN(MAC, ESE_MAC, GPKE_MAC)')
    set(b, {'CData'}, get(b, {'ZData'}));
    view([-90 90]);
    
end

%Collapse figure handles
hF = [hF(:) ; hF_(:)];

end

function [KE, SE, AllFEM] = extractModeEnergyBreakdown(h5Results, AllFEM, rMatrix, cid)

if isfield(h5Results, 'KINETIC_ENERGY')
    NasKE = h5Results.KINETIC_ENERGY;
else
   NasKE = []; 
end
if isfield(h5Results.ENERGY, 'STRAIN_ELEM')
    NasSE = h5Results.ENERGY.STRAIN_ELEM;
else
    NasSE = [];
end

%Manually remove non-flexible parts
GObj   = [AllFEM.GeometryObject];
idx    = ~contains({GObj.Name}, {'Fuselage', 'HTP', 'VTP', '765-095-RD'});
AllFEM = AllFEM(idx);

%How many parts in the assembly?
nPart  = numel(AllFEM);
if isfield(h5Results, 'SUMMARY') && isfield(h5Results.SUMMARY, 'EIGENVALUE')
    nModes = numel(h5Results.SUMMARY.EIGENVALUE.MODE);
else
    nModes = numel(h5Results.EIGENVECTOR);
end

%Preallocate
KE = cell(1, nPart);
SE = cell(1, nPart);

%Grab kinetic energy and strain energy for each beam element and node
for iP = 1 : nPart
    
    %Grab 'awi.fe' objects    
    if isempty(AllFEM(iP).Beams)
        continue
    end    
    [Beams, BeamNodes] = getSortedBeamData(AllFEM(iP));
    
    %Get ID numbers
    beamID    = [Beams.ID];
    nodeID    = [BeamNodes.ID];
    
    %Output coordinate system?
    csys  = [BeamNodes.CD];
    ucsys = unique(csys);
    
    %Preallocate
    nNodes    = numel(BeamNodes);
    ke = zeros(nNodes, nModes, 6);
    se = zeros(nNodes - 1, nModes);
    
    %Extract KE and SE for each mode
    for iM = 1 : nModes
        
        %Kinetic energy
        if ~isempty(NasKE)
            
            %Grab the KE related to this flexible component
            idx = ismember(NasKE(iM).ID, nodeID);
            keByDof = [ ...
                NasKE(iM).KET1(idx), NasKE(iM).KET2(idx), NasKE(iM).KET3(idx), ...
                NasKE(iM).KER1(idx), NasKE(iM).KER2(idx), NasKE(iM).KER3(idx)];
            
            %Transform into basic (global) coordinate system
            %TODO - Decide whether 'abs' term is necessary.
            for iCSys = 1 : numel(ucsys)
                %Grab roation matrix
                idx = ismember(csys, ucsys(iCSys));
                rMat = rMatrix{ismember(cid, ucsys(iCSys))};
                %KE in directions (T1,T2,T3)
                ke_ = keByDof(idx, 1 : 3);
                ke_ = ke_ * rMat;
                keByDof(idx, 1 : 3) = abs(ke_);
                %KE in directions (R1,R2,R3)
                ke_ = keByDof(idx, 4 : 6);
                ke_ = ke_ * rMat;
                keByDof(idx, 4 : 6) = abs(ke_);
            end
            
            ke(:, iM, :) = permute(keByDof, [1, 3, 2]);
            
        end
        
        %Strain energy
        idx = ismember(NasSE(iM).ID, beamID);
        idxE = ismember(beamID, NasSE(iM).ID);
%         se(idxE, iM) = NasSE(iM).ENERGY(idx);
        se(idxE, iM) = NasSE(iM).PCT(idx);
        
    end
    
    %Assign to cell
    KE{iP} = ke;
    SE{iP} = se;
    
end

%Filter out any parts that have no Beam elements
idx = ~cellfun(@isempty, KE);
KE  = KE(idx);
SE  = SE(idx);
AllFEM = AllFEM(idx);

end

function [Beams, BeamNodes] = getSortedBeamData(FEM)
%getSortedBeamData Grabs the nodes which are related to beams but makes
%sure they are sorted from their 'root' to 'tip' using the underlying
%geometry object behaviour.

GObj      = FEM.GeometryObject;
Beams     = FEM.Beams;
BeamNodes = FEM.BeamNodes;
BeamNodes = unique(BeamNodes(:));

%Sort be span vector
X      = [BeamNodes.X];
NodesA = [Beams.Nodes];
XA = [NodesA(1, :).X];
if isa(GObj, 'awi.model.LiftingSurface')
    %Lifting Surface
    if GObj.Span < 1
        direction = 'descend';
    else
        direction = 'ascend';
    end
    switch GObj.SpanVector
        case 'Y'
            [~, ind]  = sort(X(2, :), direction);
            [~, indA] = sort(XA(2, :), direction);
        case 'Z'
            [~, ind]  = sort(X(3, :) , direction);
            [~, indA] = sort(XA(3, :), direction);
    end
else
    %Bluff-Body
    [~, ind]  = sort(X(1, :) , 'descend');
    [~, indA] = sort(XA(2, :), 'descend');
end

%Sort and return
BeamNodes = BeamNodes(ind);
Beams     = Beams(indA);


end
