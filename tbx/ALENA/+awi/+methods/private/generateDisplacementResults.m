function DispResults = generateDisplacementResults(H5Data, FEM)
%generateDisplacementResults Defines the 'awi.results.Displacement' objects
%using the Nastran results data in 'H5Data'.

rFields = {'eigenvector', 'displacement'};

DispResults = [];

%Any results types?
resNames = fieldnames(H5Data.ResultSets);
% rFields = rFields(ismember(rFields, resNames));
resNames = resNames(ismember(lower(resNames), rFields));

%Escape route
if isempty(rFields)
    return
end
if ~isfield(H5Data.Raw.NASTRAN, 'INPUT')
    return
else
    InputData = H5Data.Raw.NASTRAN.INPUT;
end

%Get the displacements in the Nastran 'basic' coordinate system
DispData = global2BasicDisplacements(H5Data, InputData, resNames{1});

%Get the 'awi.fe.FEModel'
FEM = flatlist(FEM);

%Sort the 'awi.fe.Node' objects in ascending order to match Nastran output
Nodes    = [FEM.Nodes];
nodeIDs  = [Nodes.ID];
if isempty(nodeIDs)
    return
end
[nodeIDs, ind] = sort(nodeIDs, 'ascend');
Nodes = Nodes(ind);

%Make the results objects
DispResults = arrayfun(@(~) awi.results.NodalResults, 1 : numel(DispData));
set(DispResults, 'ResultsType', 'Displacement');

%Assign data
for iRes = 1 : numel(DispResults)
    
    %Match the results to the associated Node
    idx     = nodeIDs' == DispData(iRes).ID(1, :, 1);
    indNode = arrayfun(@(iN) find(idx(iN, :)), 1 : numel(Nodes));
    
    %Assign the data 
    set(DispResults(iRes), 'Nodes'      , Nodes);   
    set(DispResults(iRes), 'Translation', DispData(iRes).Translation(:, indNode, :));
    set(DispResults(iRes), 'Rotation'   , DispData(iRes).Rotation(:, indNode, :));
    set(DispResults(iRes), 'TimeVector' , DispData(iRes).TimeFreqEig);
    
end

end

function DispResults = global2BasicDisplacements(H5Data, InputData, resNam)
%global2BasicDisplacements Transforms the Nastran nodal data from the
%global (local) coordinate system into the basic (global) coordinate
%system using the coordinate systems defined in 'InputData'.

%Get rotation matrices
[rMatrix, cid] = getCoordSysRMatrix(InputData);

%Grab grid ID and the ID of their output coordinate systems
id_cid_map = [InputData.NODE.GRID.ID, InputData.NODE.GRID.CD];
if isfield(InputData.NODE, 'SPOINT')
    sid        = InputData.NODE.SPOINT.ID;
    id_cid_map = [id_cid_map  ; [sid, zeros(size(sid))]];
end
cid_out    = unique(id_cid_map(:, 2));
cid_out    = cid_out(cid_out ~= 0);

%What is the results data?
ResData = H5Data.ResultSets.(resNam);

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