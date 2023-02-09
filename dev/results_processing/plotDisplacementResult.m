function hF = plotDisplacementResult(FEM, ResultsObj, varargin)
%plotDisplacementResult Plots the displacements of a generic FEM.
%
% Inputs:
%   - 'FEM'        : Finite Element Model
%   - 'DispResult' : Displacement results
%
% Parameters:
%   - 'DrawSurface':  true | false (Default = false)
%   Draws a deformed surface for each displacement result by joining the
%   leading and trailing edge nodes of any lifting surfaces in the FEM
%   collection.
%   - 'PlotScaleFactor' : scalar | real (Default = 1)
%   Defines the scale factor which all displacements will be multiplied by.
%   - 'PlotIndex' : row | integer | positive (Default = [])
%   Defines the index number of the results sets which will be plotted.

%Parse
if nargin < 1
    return
else
    assert(isa(FEM, 'awi.fe.FEModel'), ['Expected the first input ', ...
        'argument to be a valid ''awi.fe.FEModel'' object.']);
end
p= inputParser;
addParameter(p, 'DrawSurface'    , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
addParameter(p, 'PlotScaleFactor', 1    , @(x)validateattributes(x, {'numeric'}, {'scalar', 'real', 'finite'}));
addParameter(p, 'PlotIndex'      , []   , @(x)validateattributes(x, {'numeric'}, {'row', 'integer', 'positive'}));
parse(p, varargin{:});

%What have we got?
if ischar(ResultsObj)                            % Input file     
    %Is it a file?
    [~, ~, ext] = fileparts(ResultsObj);
    if isempty(ext)
        %Cannot proceed
        return
    end
    %What file?
    switch ext
        case '.h5'
            ResultsObj = getResultsDataFromH5(ResultsObj);
        case '.f06'
            error('Update code');
            ResultsObj = getResultsDataFromF06(ResultsObj);
    end
end
if isstruct(ResultsObj)                          % Structure      
    %Parse feldnames
    expNam  = {'Translation', 'Rotation', 'ID', 'TimeFreqEig', 'Subcase'};
    fNames  = fieldnames(ResultsObj(1));
    assert(all(ismember(expNam, fNames)), sprintf(['When defining the ' , ...
        'displacment results as a MATLAB structure the following fields', ...
        'are required:', repmat('\n\t- %s', [1, numel(expNam)])], expNam{:}));
end
if isa(ResultsObj, 'awi.fe.result.Displacement') % Results object 
    %Check it is valid
    error('Update code');
end

%Flatlist the collection
AllFEM = unique(flatlist(FEM));

%Grab all Nodes in the FEM
Nodes = sortNodes(AllFEM);
Beams = [AllFEM.Beams];

%Grab planform nodes
if p.Results.DrawSurface
    NodeStruct = grabPlanformNodes(AllFEM);
else
    NodeStruct = [];
end

%Get plot options
%   - Plot Index
pInd = p.Results.PlotIndex;
if isempty(pInd)
    %Plot everything
    pInd = 1 : size(ResultsObj.Translation, 3);
else
    %Check the idnexing is within bounds
    nRes = size(ResultsObj.Translation, 3);
    expInd = (1 : nRes)';
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
%   - Scale Factor
sf = p.Results.PlotScaleFactor;

%Plot the deformation results
if p.Results.DrawSurface
    hF = drawDeformedSurface(ResultsObj, pInd, sf, NodeStruct, Beams);
else
    hF = drawDeformedModel(ResultsObj, pInd, sf, Nodes, Beams);
end

%Update plot apperance
formatPlotAppearance(ResultsObj, pInd, hF);

end

%% Processing results

function DispResults = getResultsDataFromH5(h5File)
%getResultsDataFromH5 Extracts the results data from a h5File and returns a
%valid results object.

rFields = {'EIGENVECTOR', 'DISPLACEMENT'};

DispResults = [];

%Get the data
[h5Raw, ~, h5Results] = h5extract(h5File);

%Any results types?
resNames = fieldnames(h5Results);
rFields = rFields(ismember(rFields, resNames));

%Escape route
if isempty(rFields)
    return
end
if ~isfield(h5Raw.NASTRAN, 'INPUT')
    return
else
    InputData = h5Raw.NASTRAN.INPUT;
end

%Get rotation matrices
[rMatrix, cid] = getCoordSysRMatrix(InputData);

    function [rMatrix, cid] = getCoordSysRMatrix(InputData)
        %getCoordSysRMatrix Extracts the rotation matrices from the .h5
        %model data.
        
        %Grab rotation matrices
        RData     = InputData.COORDINATE_SYSTEM;
        rotData   = RData.TRANSFORMATION.RDATA.DATA;
        nData     = numel(rotData);
        nCoordSys = nData / 12;
        validateattributes(nCoordSys, {'numeric'}, {'scalar', 'integer'});
        rotData   = reshape(rotData, [12, nCoordSys]);
        cid       = RData.TRANSFORMATION.IDENTITY.CID;
        
        %Index into rotation data to extract rotation matrix and origin of coordsys
        r0 = rotData(1 : 3, :);
        rMatrix = rotData(4 : end, :);
        rMatrix = arrayfun(@(i) reshape(rMatrix(:, i), [3, 3]), 1 : nCoordSys, 'Unif', false);
        
    end

%Grab grid ID and the ID of their output coordinate systems
id_cid_map = [InputData.NODE.GRID.ID, InputData.NODE.GRID.CD];
if isfield(InputData.NODE, 'SPOINT')
    sid        = InputData.NODE.SPOINT.ID;
    id_cid_map = [id_cid_map  ; [sid, zeros(size(sid))]];
end
cid_out    = unique(id_cid_map(:, 2));
cid_out    = cid_out(cid_out ~= 0);

%What is the results data?
ResData = h5Results.(rFields{1});

%Preallocate
nNode = numel(id_cid_map) / 2;
nMode = numel(ResData);
trans = zeros(nNode, 3, nMode);
rotat = zeros(nNode, 3, nMode);
id    = zeros(nNode, 1, nMode);

%Process displacements
%   - Transform all displacements into basic coordinate system
for i = 1 : nMode
    %Stash a copy of the displacement results and the node IDs
    temp        = ResData(i);
    id(:, 1, i) = temp.ID;
    %Construct translation (T) and rotation (R) results in 3xn matrices
    trans(:, :, i) = [temp.X , temp.Y , temp.Z];
    rotat(:, :, i) = [temp.RX, temp.RY, temp.RZ];
    for iR = 1 : numel(cid_out)
        %Select results in this coord sys
        rMat = rMatrix{cid == cid_out(iR)};
        idx  = ismember(id_cid_map(:, 2), cid_out(iR));
        trans(idx, :, i) = trans(idx, :, i) * rMat;
    end
end

%Switch order of dims to match expected data in Node.GlobalTranslation
trans = permute(trans, [2, 1, 3]);
rotat = permute(rotat, [2, 1, 3]);
id    = permute(id   , [2, 1, 3]);

%Assign to output structure
DispResults.Translation = trans;
DispResults.Rotation    = rotat;
DispResults.ID          = id;
DispResults.TimeFreqEig = [ResData.TIME_FREQ_EIGR];
DispResults.Subcase     = [ResData.SUBCASE];

end

function Nodes = sortNodes(AllFEM)
%sortNodes Grabs all the 'awi.fe.Node' objects in the FEM collection and
%sorts them by ID number.

Nodes  = [AllFEM.Nodes];
[~, index] = sort([Nodes.ID], 'ascend');
Nodes  = Nodes(index);

end

function NodeStruct = grabPlanformNodes(AllFEM)
%grabPlanformNodes Retrieves the 'awi.fe.Node' objects that lie on the
%leading or trailing edge of any lifting surface planforms. Also returns
%the nodes which the FE beam elements are attached to.

NodeStruct = [];

%Dependent on availability of 'LiftingSurface' related FEMs
GObj        = [AllFEM.GeometryObject];
% LiftSurfFEM = AllFEM(arrayfun(@(o) isa(o, 'awi.model.LiftingSurface'), GObj));
LiftSurfFEM     = AllFEM(arrayfun(@(o) isa(o, 'awi.model.Beam'), GObj));
if isempty(LiftSurfFEM)
    return
end

%Preallocate
nLS       = numel(LiftSurfFEM);
BeamNodes = cell(1, nLS);
LENodes   = cell(1, nLS);
TENodes   = cell(1, nLS);

%Find beam and planform nodes for every lifting surface
for i = 1 : nLS
    BeamNode = LiftSurfFEM(i).BeamNodes;
    BeamNode = unique(BeamNode(:));
    if isempty(BeamNode)
        continue
    end
    if isa(LiftSurfFEM(i).GeometryObject, 'awi.model.LiftingSurface')
        %Assume planform nodes are connected to the beam by 'RigidBar' elements
        RigidBar = LiftSurfFEM(i).RigidBars;
        if isempty(RigidBar)
            continue
        end
        RigidBar = RigidBar(ismember(vertcat(RigidBar.NodesI), BeamNode));
        %TODO - Update this check because this assumes the only rigid elements
        %attached to the beam nodes are the planform nodes.
        PlanformNodes = horzcat(RigidBar.NodesD);
        LENodes{i} = PlanformNodes(1, :)';
        TENodes{i} = PlanformNodes(2, :)';
    end
    %Sort nodes in ascending order
    [~, ind]     = sort([BeamNode.ID], 'ascend');
    BeamNodes{i} = BeamNode(ind);
end

%Assign to structure
NodeStruct.LENodes   = LENodes;
NodeStruct.TENodes   = TENodes;
NodeStruct.BeamNodes = BeamNodes;

end

%% Plotting

function hF = drawDeformedModel(DispResults, pInd, scaleFactor, Nodes, Beams)
%drawDeformedModel Plots the deformations for every node in the model.

nRes = numel(pInd);

%Stash the nodeIDs
nodeIDs = [Nodes.ID];

%Mak graphics objects
hF  = arrayfun(@(i) figure, 1 : nRes);
hAx = arrayfun(@(i) axes('Parent', hF(i), 'NextPlot', 'add'), 1 : nRes);

for i = 1 : nRes
    
    rIndex = pInd(i);
    
    %Grab current translation results and scale
    dT = DispResults.Translation(:, :, rIndex) .* scaleFactor;
    dR = DispResults.Translation(:, :, rIndex) .* scaleFactor;
    id = DispResults.ID(1, :, rIndex);
    
    %Only assign data for nodes that have been provided
    idx = ismember(id, nodeIDs);
    
    %Assign to the Node objects
    set(Nodes, {'GlobalTranslation'}, num2cell(dT(:, idx), 1)');
    set(Nodes, 'DrawMode', 'deformed');
    
    %Draw the deformed model
    drawElement(Nodes, hAx(i));
    drawElement(Beams, hAx(i));
    
end

end

function hF = drawDeformedSurface(DispResults, pInd, scaleFactor, NodeStruct, Beams)
%drawDeformedSurface Plots a deformed surface using the displacement
%results 'DispResults'.
%
% Uses the Node handles in 'NodeStruct' to grab the displacement results
% for the LE and TE of any Lifting Surfaces in the collection. This is then
% used to create a deformed surface which can show the model twist.

%Colour and transparency for the undeformed plot
clr   = [220, 220, 200] ./ 255;
alfa  = 0.4;

nRes = numel(pInd);

%Grab node data
BeamNodes = NodeStruct.BeamNodes;
BNode     = vertcat(BeamNodes{:});
LENodes   = NodeStruct.LENodes;
TENodes   = NodeStruct.TENodes;
nLS       = numel(TENodes);

%Grab all nodes referenced by the Beam objects
AllBeamNodes = unique([Beams.Nodes]);
set([Beams.Nodes], 'DrawMode', 'deformed');


%Make graphics objects
hF  = arrayfun(@(i) figure, 1 : nRes);
hAx = arrayfun(@(i) axes('Parent', hF(i), 'NextPlot', 'add'), 1 : nRes);

for i = 1 : nRes
    
    rIndex = pInd(i);
    ht     = hAx(i);
    
    %Grab current translation results and scale
    dT = DispResults.Translation(:, :, rIndex) .* scaleFactor;
    dR = DispResults.Translation(:, :, rIndex);
    id = DispResults.ID(1, :, rIndex);
    
    %Take only the real parts - assuming complex
    dT = real(dT);
    dR = real(dR);
    
    %Grab the planform coordinates and beam twist
    [r0, rDeform, tws, rBeamDeform]  = grabPlanformData(dT, dR, id, BeamNodes, LENodes, TENodes);
    
    %Plot the patch (deformed & undeformed)    
    arrayfun(@(i) patch(ht, r0{i}(1, :)     , r0{i}(2, :)     ,  r0{i}(3, :)     , ...
        clr, 'FaceAlpha', alfa), 1 : numel(rDeform));
    arrayfun(@(i) patch(ht, rDeform{i}(1, :), rDeform{i}(2, :),  rDeform{i}(3, :), ...
        tws{i}), 1 : numel(rDeform));
    
    %Set up deformation data for the beam nodes
    dTBeam = horzcat(rBeamDeform{:});
    set(BNode, {'GlobalTranslation'}, num2cell(dTBeam, 1)');
    
    %Plot the displaced beams
    drawElement(AllBeamNodes, hAx(i), 'deformed');
    drawElement(Beams, hAx(i), 'deformed');     
    
%     colorbar(ht);
    
    %Add column of nan data to allow vectorisation
    %coords = cellfun(@(x) [x, x(:, 1), nan(3, 1)], coords, 'Unif', false);
    %tws    = cellfun(@(x) [x, x(1, 1), nan(1, 1)], tws   , 'Unif', false);
    
    %Can plot all in one but cannot fill becauase of nan data
    %coords = horzcat(coords{:});
    %tws    = horzcat(tws{:});
    %hF = figure;
    %hAx = axes('Parent', hF);
    %hP_ = patch(coords(1, :), coords(2, :), coords(3, :), 'b');
    
end

%Local functions

    function [r0, rDeform, tws, rBeamDeform] = grabPlanformData(trans, rotat, id, BeamNodes, LENodes, TENodes)
        
        %Preallocate
        %         coords = zeros(3, nCoord);
        %         tws    = zeros(1, nCoord);
        nLS         = numel(BeamNodes);        
        r0          = cell(1, nLS);
        rDeform     = cell(1, nLS);
        rBeamDeform = cell(1, nLS);
        tws         = cell(1, nLS);
        
        %Extract planform deformed coords and twist for each lifting surf
        for ii = 1 : nLS
            
            %Grab beam data
            bmID = get(BeamNodes{ii}, 'ID');
            bmID = horzcat(bmID{:});
            drBM = trans(:, ismember(id, bmID));
            rBeamDeform{ii} = drBM;
            
            if isempty(LENodes{ii})
                continue
            end
            
            %Undeformed coordinates            
            rLE = [LENodes{ii}.X];
            rTE = [TENodes{ii}.X];
            
            %Stash the undeformed shape
            r0{ii} = [rLE, fliplr(rTE)];  
            
            %Get ID numbers of LE/TE nodes            
            leID = get(LENodes{ii}  , 'ID');
            teID = get(TENodes{ii}  , 'ID');            
            leID = horzcat(leID{:});
            teID = horzcat(teID{:});
            
            %Get translations of these nodes            
            drLE = trans(:, ismember(id, leID));
            drTE = trans(:, ismember(id, teID));
            
            %Get displaced coords
            rLE  = rLE + drLE;
            rTE  = rTE + drTE;
            
            %Get rotation (twist) of beam nodes about beam axis
            twist   = rad2deg(rotat(1, ismember(id, bmID)));
            tws{ii} = [twist, fliplr(twist)];
%             aoa     = atand(rLE(3, :) - rTE(3, :) ./ abs(rLE(1, :) - rTE(1, :)));
%             tws{ii} = [aoa, fliplr(aoa)];

            %Construct path coords            
            rDeform{ii}     = [rLE, fliplr(rTE)];
            
        end
        
        %Remove empties
        tws     = tws(~cellfun(@isempty, tws));
        r0      = r0(~cellfun(@isempty, r0));
        rDeform = rDeform(~cellfun(@isempty, rDeform));
                        
    end

end

function formatPlotAppearance(DispResults, pInd, hF)

%Get axes handles
hAx = findobj(hF, 'Type', 'axes');

%Format
axis(hAx, 'equal');
set([hAx.XLabel], 'String', 'X [m]');
set([hAx.YLabel], 'String', 'Y [m]');
set([hAx.ZLabel], 'String', 'Z [m]');
set(hAx         , 'View'  , [-60, 20]);

%Assume we are dealing with mode data
lambda = DispResults.TimeFreqEig(pInd);
omega  = sqrt(lambda);
freq   = omega ./ 2 ./ pi;
tit    = arrayfun(@(i) sprintf('Mode %i - Freq. = %.3fHz', pInd(i), freq(i)), 1 : numel(pInd), 'Unif', false);
set([hAx.Title], {'String'}, tit');

%Rename the figure
set(hF, {'Name'}, tit');

end
