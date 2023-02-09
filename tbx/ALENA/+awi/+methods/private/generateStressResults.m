function StressResults = generateStressResults(H5Data, FEM)
%generateDisplacementResults Defines the 'awi.results.Stress' objects
%using the Nastran results data in 'H5Data'.

rFields = {'eigenvector', 'stress'};

StressResults = [];

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

%Get the stresses in the Nastran
StressData = H5Data.ResultSets.STRESS.BEAM;
StressData_NL = H5Data.ResultSets.STRESS.BEAM_NL; % Nonlinear stress data

%Get the 'awi.fe.FEModel'
FEM = flatlist(FEM);

%Sort the elements in ascending order to match Nastran output
Elements    = [FEM.Beams];
ElementIDs  = [Elements.ID];
if isempty(ElementIDs)
    return
end
[ElementIDs, ind] = sort(ElementIDs, 'ascend');
Elements = Elements(ind);

%Make the results objects
StressResults = arrayfun(@(~) awi.results.ElementResults, 1 : numel(StressData));
set(StressResults, 'ResultsType', 'Stress');

%Assign data
for iRes = 1 : numel(StressResults)
    
    %Match the results to the associated Element
    idx     = ElementIDs' == StressData(iRes).EID';
    indNode = arrayfun(@(iN) find(idx(iN, :)), 1 : numel(Elements));
    
    %Assign the data 
    set(StressResults(iRes), 'Elements'      , Elements);
    set(StressResults(iRes), 'EID', StressData(iRes).EID(indNode, :));
    set(StressResults(iRes), 'GRID', StressData(iRes).GRID(:, indNode, :));
    set(StressResults(iRes), 'SD', StressData(iRes).SD(:, indNode, :));
    set(StressResults(iRes), 'XC', StressData(iRes).XC(:, indNode, :));
    set(StressResults(iRes), 'XD', StressData(iRes).XD(:, indNode, :));
    set(StressResults(iRes), 'XE', StressData(iRes).XE(:, indNode, :));
    set(StressResults(iRes), 'XF', StressData(iRes).XF(:, indNode, :));
    set(StressResults(iRes), 'MAX', StressData(iRes).MAX(:, indNode, :));
    set(StressResults(iRes), 'MIN', StressData(iRes).MIN(:, indNode, :));
    
    set(StressResults(iRes), 'GRIDA', StressData_NL(iRes).GRIDA(indNode, :));
    
    set(StressResults(iRes), 'NSXCA', StressData_NL(iRes).NSXCA(indNode, :));
    set(StressResults(iRes), 'NSECA', StressData_NL(iRes).NSECA(indNode, :));
    set(StressResults(iRes), 'TECA', StressData_NL(iRes).TECA(indNode, :));
    set(StressResults(iRes), 'EPECA', StressData_NL(iRes).EPECA(indNode, :));
    set(StressResults(iRes), 'ECECA', StressData_NL(iRes).ECECA(indNode, :));
    
    set(StressResults(iRes), 'NSXDA', StressData_NL(iRes).NSXDA(indNode, :));
    set(StressResults(iRes), 'NSEDA', StressData_NL(iRes).NSEDA(indNode, :));
    set(StressResults(iRes), 'TEDA', StressData_NL(iRes).TEDA(indNode, :));
    set(StressResults(iRes), 'EPEDA', StressData_NL(iRes).EPEDA(indNode, :));
    set(StressResults(iRes), 'ECEDA', StressData_NL(iRes).ECEDA(indNode, :));
    
    set(StressResults(iRes), 'NSXEA', StressData_NL(iRes).NSXEA(indNode, :));
    set(StressResults(iRes), 'NSEEA', StressData_NL(iRes).NSEEA(indNode, :));
    set(StressResults(iRes), 'TEEA', StressData_NL(iRes).TEEA(indNode, :));
    set(StressResults(iRes), 'EPEEA', StressData_NL(iRes).EPEEA(indNode, :));
    set(StressResults(iRes), 'ECEEA', StressData_NL(iRes).ECEEA(indNode, :));
    
    set(StressResults(iRes), 'NSXFA', StressData_NL(iRes).NSXFA(indNode, :));
    set(StressResults(iRes), 'NSEFA', StressData_NL(iRes).NSEFA(indNode, :));
    set(StressResults(iRes), 'TEFA', StressData_NL(iRes).TEFA(indNode, :));
    set(StressResults(iRes), 'EPEFA', StressData_NL(iRes).EPEFA(indNode, :));
    set(StressResults(iRes), 'ECEFA', StressData_NL(iRes).ECEFA(indNode, :));
    
    set(StressResults(iRes), 'GRIDB', StressData_NL(iRes).GRIDB(indNode, :));

    set(StressResults(iRes), 'NSXCB', StressData_NL(iRes).NSXCB(indNode, :));
    set(StressResults(iRes), 'NSECB', StressData_NL(iRes).NSECB(indNode, :));
    set(StressResults(iRes), 'TECB', StressData_NL(iRes).TECB(indNode, :));
    set(StressResults(iRes), 'EPECB', StressData_NL(iRes).EPECB(indNode, :));
    set(StressResults(iRes), 'ECECB', StressData_NL(iRes).ECECB(indNode, :));

    set(StressResults(iRes), 'NSXDB', StressData_NL(iRes).NSXDB(indNode, :));
    set(StressResults(iRes), 'NSEDB', StressData_NL(iRes).NSEDB(indNode, :));
    set(StressResults(iRes), 'TEDB', StressData_NL(iRes).TEDB(indNode, :));
    set(StressResults(iRes), 'EPEDB', StressData_NL(iRes).EPEDB(indNode, :));
    set(StressResults(iRes), 'ECEDB', StressData_NL(iRes).ECEDB(indNode, :));

    set(StressResults(iRes), 'NSXEB', StressData_NL(iRes).NSXEB(indNode, :));
    set(StressResults(iRes), 'NSEEB', StressData_NL(iRes).NSEEB(indNode, :));
    set(StressResults(iRes), 'TEEB', StressData_NL(iRes).TEEB(indNode, :));
    set(StressResults(iRes), 'EPEEB', StressData_NL(iRes).EPEEB(indNode, :));
    set(StressResults(iRes), 'ECEEB', StressData_NL(iRes).ECEEB(indNode, :));

    set(StressResults(iRes), 'NSXFB', StressData_NL(iRes).NSXFB(indNode, :));
    set(StressResults(iRes), 'NSEFB', StressData_NL(iRes).NSEFB(indNode, :));
    set(StressResults(iRes), 'TEFB', StressData_NL(iRes).TEFB(indNode, :));
    set(StressResults(iRes), 'EPEFB', StressData_NL(iRes).EPEFB(indNode, :));
    set(StressResults(iRes), 'ECEFB', StressData_NL(iRes).ECEFB(indNode, :));
    
end

end