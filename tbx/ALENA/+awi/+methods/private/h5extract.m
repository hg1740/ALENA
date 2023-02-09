function [h5Data, varargout] = h5extract(filename)
%h5extract Extracts the data from a .h5 file into a Matlab structure.
%
% h5Data = h5extract() will open up a dialogue box prompting the user to
% select a .h5 file. The data will then be extracted from that file.
%
% h5Data = h5extract(filename) will extract all the data from the .h5 file
% given by the path 'filename'.
%
% [h5Data, Meta] = h5extract(...) will return the meta data for the .h5
% file. This is the same as the output from h5info(filename).
%
% Detailed Explanation:
%   - The Matlab structure 'h5Data' preserves the hierachy of the .h5 file
%     being read.
%   - If any quantity is terminated by '_CPLX' then it is assumed to
%     contain complex data. The real and imaginary parts are selected by
%     identifying repeat fields that end in 'R' and 'I'.
%
% References:
%   [1]. MSC.Nastran Reference Manual - "The Nastran HDF5 Result Database
%   (NH5RDB)"
%
% Author    : Christopher Szczyglowski
% Username  : cs1414
% Email     : cs1414@bristol.ac.uk
% Timestamp : 30-Nov-2017 14:55:25
%
% Copyright (c) 2017 Christopher Szczyglowski
% All Rights Reserved
%
%
% Revision: 1.0 30-Nov-2017 14:55:25
%   - Initial function :
%

h5Data = [];

%Prompt user if no file is provided
if nargin == 0
    filename = getfile({'*.h5', 'HDF5 File'}, 'Select a file');
    if isempty(filename), return; end %Escape route
end

%Check extension
[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.h5'), '''Filename'' must be the name of a .h5 file');

%Grab meta data
Meta = h5info(filename);

%Extract data from the file
h5Data = extractData(filename, Meta, struct());

%Additional outputs
if nargout > 1 % h5 meta data   
    varargout{1} = Meta;
end
if nargout > 2 % Results set (structure)
    %Consolidate results sets
    %   - N.B. The behaviour of this function is not guaranteed. Recommend
    %    commenting this out.
    varargout{2} = generateResultsSet(h5Data);
end

end

function Data = extractData(filename, MetaGroup, varargin)
%findDataSet Extracts all the data from the .h5 file 'filename' into a
%Matlab structure which preserves the .h5 data hierachy.
%
% Assumptions:
%   - Assume that a structure with a non-empty 'Datasets' field is a root
%     structure. The structure above this one is the parent containing the
%     group name.

p = inputParser;
addRequired(p, 'filename' , @ischar);
addRequired(p, 'MetaGroup', @validateMetaGroup);
addOptional(p, 'Data'     , struct(), @isstruct);
parse(p, filename, MetaGroup, varargin{:});

%Grab current data
Data = p.Results.Data;

%If MetaGroup.DataSet is populated then extract data
if ~isempty(MetaGroup.Datasets)
    names = {MetaGroup.Datasets.Name};
    for i = 1 : numel(names)
        %Formulate 'groupset'
        groupset = strcat([MetaGroup.Name, '/'], names{i});
        %Use built-in 'h5read' to extract data
        h5Data   = h5read(filename, groupset);
        %Check for complex data using the token '_CPLX'
        if contains(names{i}, '_CPLX')
            %Remove '_CPLX' from the results set name
            rName = names{i}(1 : strfind(names{i}, '_CPLX') - 1);
            %Check to make sure we don't override an existing field
            if isfield(Data, rName)
                rName = [rName, '_CPLX'];
            end
            Data.(rName) = formatComplexData(h5Data);
        else
            Data.(names{i}) = h5Data;
        end
    end
end

%If MetaGroup.Groups is populated then ...
%... recursively loop through groups and find Datasets
if ~isempty(MetaGroup.Groups)
    gNames = formatGroupName({MetaGroup.Groups.Name});
    for iG = 1 : numel(MetaGroup.Groups)
        Data.(gNames{iG}) = extractData(filename, MetaGroup.Groups(iG));
    end
end

    function tf = validateMetaGroup(mg)
        
        validMetaGroupFields = {'Filename', 'Name', 'Groups', 'Datasets', 'Datatypes', 'Links', 'Attributes'};
        tf = true;
        
        validateattributes(mg, {'struct'}, {'nonempty'}, 'findDataSet', 'Meta');
        index = ~ismember(fieldnames(mg), validMetaGroupFields);
        
        if any(index)
            ME = MException('matlab:mni:NH5RDB:badFormat', ...
                ['The fieldnames of ''MetaGroup'' do not match their ', ...
                'expected values.\n''MetaGroup'' should have the '    , ...
                'following fields:\n\n\t%s\n\nInstead it had fields:' , ...
                '\n\n\t%s\n\n'], strjoin(validMetaGroupFields), ...
                strjoin(fieldnames(mg)));
            throwAsCaller(ME);
        end
        
    end

    function groupName = formatGroupName(groupString)
        %formatGroupName Grabs the final group name from the group set
        %string.
        
        groupName = cellfun(@(x) x(find(x == '/', 1, 'last') + 1 : end), ...
            groupString, 'Unif', false);
        
    end

    function ComplexData = formatComplexData(Data)
        %formatComplexData Returns the data in 'Data' in a complex format
        %where appropriate by selecting the real and imaginary data from
        %the fields in 'Data'.
        %
        % Real and imaginary data is selected by searching for repeat
        % fieldnames that are appended by 'R' or 'I'.
        
        %Grab field names
        fNames = fieldnames(Data);
        
        %Remove final character and search for unique names
        cropNames   = cellfun(@(x) x(1 : end - 1), fNames, 'Unif', false);
        uniqueNames = unique(cropNames);
        
        %If no. unique names = no. field names then no formatting needed
        if numel(uniqueNames) == numel(fNames)
            ComplexData = Data;
        else
            %Convert structure to cell for easier indexing
            dataCell = struct2cell(Data);
            %Formulate real ('R') and imaginary ('I') fieldnames
            realNames = strcat(uniqueNames, 'R');
            imagNames = strcat(uniqueNames, 'I');
            %Logical indexing
            iR = ismember(fNames, realNames);
            iI = ismember(fNames, imagNames);
            iO = ~or(iR, iI);
            %Grab real and imaginary data & seperate other data
            rData = dataCell(iR);
            iData = dataCell(iI);
            other      = dataCell(iO);
            otherNames = fNames(iO);
            dataNames  = cropNames(iR);
            %Convert to type 'complex double' - Loop is quicker than
            %vectorised use of 'complex' in this instance
            data = cell(size(rData));
            for iii = 1 : numel(dataNames)
                data{iii} = complex(rData{iii}, iData{iii});
            end
            ComplexData = cell2struct([other; data], [otherNames; dataNames]);
            %%Attempt vectorisation (Slower than loop)
            %r = cat(2, rData{:});
            %j = cat(2, iData{:});
            %data = num2cell(complex(r, j), 1)';
        end
    end

end

function ResultsSet = generateResultsSet(h5Data)
%generateResultsSet Formats the results data in the structure 'h5Data'
%into a results set.

%Just return nothing for now
ResultsSet = struct();

%Check if the 'RESULT' field is present
if ~isfield(h5Data.NASTRAN, 'RESULT') || ~isfield(h5Data, 'INDEX')
    ResultsSet = [];
    return
end

%Get name of results groups
rgNames = fieldnames(h5Data.NASTRAN.RESULT);

%SUMMARY maps straight across
if isfield(h5Data.NASTRAN.RESULT, 'SUMMARY')
    ResultsSet.SUMMARY = h5Data.NASTRAN.RESULT.SUMMARY;
    rgNames(ismember(rgNames, 'SUMMARY')) = [];
end

%Ignore 'DOMAINS' & 'ENERGY' & 'OPTIMIZATION' for now
rgNames(ismember(rgNames, 'DOMAINS')) = [];
% rgNames(ismember(rgNames, 'ENERGY'))  = [];

%OPTIMIZATION results have a special method
if any(ismember(rgNames, 'OPTIMIZATION'))
    
    %Remove from the default list
    rgNames(ismember(rgNames, 'OPTIMIZATION')) = [];
    
    %Grab results for easier indexing
    Opt        = h5Data.NASTRAN.RESULT.OPTIMIZATION;
    if isfield(h5Data.INDEX.NASTRAN.RESULT, 'OPTIMIZATION')
        AllIndex   = h5Data.INDEX.NASTRAN.RESULT.OPTIMIZATION;
    else
        AllIndex = [];
    end
    DomainData = h5Data.NASTRAN.RESULT.DOMAINS;
    
    %Extract the OPTIMIZATION results
    ResultsSet = formatOptimisationResults(Opt, AllIndex, DomainData, ResultsSet);
    
end

for iRG = 1 : numel(rgNames)
    DomainData   = h5Data.NASTRAN.RESULT.DOMAINS;
    ResultsGroup = h5Data.NASTRAN.RESULT.(rgNames{iRG});
    ResultsIndex = h5Data.INDEX.NASTRAN.RESULT.(rgNames{iRG});
    ResultsSet   = splitIntoResultsSet(ResultsGroup, ResultsIndex, DomainData, ResultsSet);
end

    function ResultsSet = formatOptimisationResults(Opt, OptIndex, DomainData, ResultsSet)
        %formatOptimisationResults Returns a structure detailing the design
        %cycle history of the objective function and all design variables.
        %If the sensitivities have been requested then these will also be
        %included in the 'ResultsSet'.        
        
        %% Deal with objective function (DESOBJ)
        %   - Simple, everything maps across.
        %   - Just need to index the design cycle!
        
        DESOBJ.EXACT = Opt.OBJECTIVE.EXACT;
        DESOBJ.APPRX = Opt.OBJECTIVE.APPRX;
        DESOBJ.MAXIM = Opt.OBJECTIVE.MAXIM;
        
        idx = ismember(DomainData.ID, Opt.OBJECTIVE.DOMAIN_ID);
        DESOBJ.DesignCycle = DomainData.DESIGN_CYCLE(idx);
        
        %% Deal with Design Variables (DESVAR)
        %   - Return a structure with one field per design variable.
        %   - Name each field after the design variable.
        %   - For each design variable provide the following information:
        %       + IDs
        %       + Bounds
        %       + Initial conditions (X0)
        %       + Design Cycle history (Xi)
        
        %Design variable data
        dvarData  = struct2cell(Opt.LABEL);
        dvarData_ = fieldnames(Opt.LABEL);
        
        %Grab design variable names
        dvarNames = cellstr(Opt.LABEL.LABEL');
        
        %Throw away label data
        idxL = ismember(dvarData_, 'LABEL');
        dvarData(idxL)  = [];
        dvarData_(idxL) = [];
        
        %Return a structure with one field per design variable
        DESVAR = struct();
        for iD = 1 : numel(dvarNames)
            
            %Assign basic data
            for ii = 1 : numel(dvarData)
                DESVAR.(dvarNames{iD}).(dvarData_{ii}) = dvarData{ii}(iD);
            end
            
            %Grab initial value
            idx = ismember(Opt.VARIABLE.VID, DESVAR.(dvarNames{iD}).IDVID);
            DESVAR.(dvarNames{iD}).X0 = Opt.VARIABLE.INIT(idx);
            
            %Grab design variable history for each design cycle
            idx = ismember(Opt.HISTORY.VID, DESVAR.(dvarNames{iD}).IDVID);
            DESVAR.(dvarNames{iD}).Xi = Opt.HISTORY.VALUE(idx);
            
            %Assign the correct design cycle using the 'DOMAINS' data
            domainID = Opt.HISTORY.DOMAIN_ID(idx);
            idxD = ismember(DomainData.ID, domainID);
            DESVAR.(dvarNames{iD}).DesignCycle = DomainData.DESIGN_CYCLE(idxD);
            
        end        
       
        %% Deal with Design Sensitivities (DESSENS) 
        %   - This will be more involved...
        %   - Need to look into the 'OptIndex' structure
        
        if ~isempty(OptIndex) && isfield(OptIndex, 'SENSITIVITY')
            DESSENS = struct(); %Do nothing for now...
        else            
            DESSENS = struct();
        end
        
        %% Assign to structure
        
        ResultsSet.OPTIMISATION.DESOBJ = DESOBJ;
        ResultsSet.OPTIMISATION.DESVAR = DESVAR;
        ResultsSet.OPTIMISATION.DESSENS = DESSENS;
                 
    end

    function ResultsSet = splitIntoResultsSet(ResultsGroup, ResultsIndex, DomainData, ResultsSet)
        
        %Get names of all results
        rNames = fieldnames(ResultsGroup);
        rNames(ismember(rNames, 'IDENT')) = [];  %IGNORE IDENT FOR NOW --> Part of the ENERGY field
        rNames(ismember(rNames, 'STRAIN_ELEM')) = [];
        for iR = 1 : numel(rNames)
            %Grab result & index data
            result = ResultsGroup.(rNames{iR});
            index  = ResultsIndex.(rNames{iR});
            %Have we reached the bottom level?
            if ~isfield(index, 'POSITION')
                ResultsSet.(rNames{iR}) = struct(); %APPEND 'ResultsSet'
                ResultsSet.(rNames{iR}) = splitIntoResultsSet(result, index, DomainData, ResultsSet.(rNames{iR}));
            else
                %Split by subcase
                %scID = DomainData.SUBCASE(ResultsGroup.(rNames{iR}).DOMAIN_ID);
                %[C, ia, ic] = unique(scID);
                %Define the structure arrays
                nSet    = numel(index.POSITION);
                lb      = index.POSITION + 1;
                ub      = index.POSITION + index.LENGTH;
                rqNames = fieldnames(result);
                qNames  = [rqNames ; {'SUBCASE' ; 'DESIGN_CYCLE' ; 'TIME_FREQ_EIGR'}];
                %Preallocate
                ResultsSet.(rNames{iR})(nSet, 1) = cell2struct(cell(size(qNames)), qNames);
                %Loop through subcases and defined the ResultsSets (1 per subcase)
                for iSet = 1 : nSet
                    for iQ = 1 : numel(rqNames)    
                        %Determine which dimension contains the data
                        %relating to 'index.POSITION'
                        [m, n] = size(result.(rqNames{iQ}));
                        %Index the data
                        if mod(m, nSet) == 0 && mod(n, nSet) == 0
                            try %Try both
                                data = result.(rqNames{iQ})(lb(iSet) : ub(iSet), :);
                            catch
                                data = result.(rqNames{iQ})(:, lb(iSet) : ub(iSet));
                            end
                        elseif mod(m, nSet) == 0
                            %Index the rows & take all columns
                            data = result.(rqNames{iQ})(lb(iSet) : ub(iSet), :);
                        elseif mod(n, nSet) == 0
                            %Index the columns & take all rows
                            data = result.(rqNames{iQ})(:, lb(iSet) : ub(iSet));
                        else
                            data = result.(rqNames{iQ})(lb(iSet) : ub(iSet), :);
                        end
                        %Check dimensions are correct 
                        %   - This is based on the first entry which is
                        %     typically the ID entry.
                        if iQ > 1
                            if size(data, 1) ~= size(ResultsSet.(rNames{iR})(iSet).(rqNames{1}), 1)
                               data = permute(data, [3, 2, 1]);
                            end
                        end
                        %Assign to 'ResultsSet'
                        ResultsSet.(rNames{iR})(iSet).(rqNames{iQ}) = data;
                    end
                    
                end
                %Assign domain data
                scID  = DomainData.SUBCASE(horzcat(ResultsSet.(rNames{iR}).DOMAIN_ID));
                dcID  = DomainData.DESIGN_CYCLE(horzcat(ResultsSet.(rNames{iR}).DOMAIN_ID));
                tfID  = DomainData.TIME_FREQ_EIGR(horzcat(ResultsSet.(rNames{iR}).DOMAIN_ID));
                scNum = num2cell(scID(1, :));
                dcNum = num2cell(dcID(1, :)); 
                tfNum = num2cell(tfID(1, :));                 
                [ResultsSet.(rNames{iR})(:).SUBCASE]        = deal(scNum{:});
                [ResultsSet.(rNames{iR})(:).DESIGN_CYCLE]   = deal(dcNum{:});
                [ResultsSet.(rNames{iR})(:).TIME_FREQ_EIGR] = deal(tfNum{:});                
                %Filter any content that does not appear to match subcase format
                SubcaseData = DomainData.SUBCASE;
                if nSet ~= numel(SubcaseData)
                    continue
                end
                %Split each 'ResultsSet' by subcase
                %   - Assume that each subcase has an equal number of entries
                %       + This is highly unlikely and needs to be improved upon.
                %       Simple enough to do propertly - Just need to use logical
                %       indexing and a for loop however I need a quick, dirty
                %       solution right now and a simple 'RESHAPE' fits the bill
                %   - TODO: Get rid of this.
                subcaseID = unique(SubcaseData);
                nSubcase  = numel(subcaseID);
                ResultsSet.(rNames{iR}) = reshape(ResultsSet.(rNames{iR}), ...
                    [nSet / nSubcase, nSubcase]);
            end            
        end
        
    end

end
