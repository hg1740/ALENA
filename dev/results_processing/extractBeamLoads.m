function BeamLoads = extractBeamLoads(AllFEM, h5Results)
%extractBeamLoads Extracts the beam loads for each FEM object

BeamLoads = [];

%Parse
if nargin < 2
   return 
end
if ~isfield(h5Results, 'ELEMENT_FORCE') || ~isfield(h5Results.ELEMENT_FORCE, 'BEAM')
    return
else
    ElementForce = h5Results.ELEMENT_FORCE;
    BeamForce    = ElementForce.BEAM;
end

%Collapse beam force values into a matrix to allow indexing
forceNames = {'BM1', 'BM2', 'TS1', 'TS2', 'AF', 'TTRQ'};
ForceData  = cell2struct(cell(size(forceNames)), forceNames, 2);
for iQ = 1 : numel(forceNames)
    ForceData.(forceNames{iQ}) = vertcat(BeamForce.(forceNames{iQ}));    
end

%Preallocate
BeamLoads = cell2struct(cell(1, numel(forceNames) + 1), [{'Beams'}, forceNames], 2);

%Extract beam loads for each FEM component
for iF = 1 : numel(AllFEM)
    %Grab beam objects and sort in ascending order of ID
    Beams      = AllFEM(iF).Beams;
    if isempty(Beams)
        continue
    end
    [eid, ind] = sort([Beams.ID], 'ascend');
    Beams      = Beams(ind);
    
    %Assume subcase results are defined for the same set of beams 
    resID = BeamForce(1).EID;
    
    %Grab loads for this set of beams
    %   - Loads are returned at END-A and END-B of the beam
    idx  = ismember(resID, eid);    
    temp = structfun(@(f) f(:, idx, [1, end]), ForceData, 'Unif', false);
    temp = struct2cell(temp);
    
    %Return to structure
    BeamLoads(iF) = cell2struct([{Beams} ; temp], [{'Beams'} ; forceNames']);
    
end

%Only return populated data structures
BeamLoads = BeamLoads(arrayfun(@(s) ~isempty(s.Beams), BeamLoads));


end

function x = padCoordsWithNaN(x)
%padCoordsWithNaN Accepts a matrix of [2, N] sets of
%coordinates which represent the coordinate of a series of
%lines from end-A to end-B and returns a single vector with all
%of the coordinates padded by NaN terms.
%
% This function enables the plotting of line objects to be
% vectorised.

%Convert to cell so we retain the pairs of coordinates in the
%correct order
x  = num2cell(x, 1);

%Preallocate
x_ = cell(1, 2 * numel(x));

%Assign the data and NaN terms
x_(1 : 2 : end - 1) = x;
x_(2 : 2 : end)     = {nan};

%Return a column vector
x = vertcat(x_{:});

end

