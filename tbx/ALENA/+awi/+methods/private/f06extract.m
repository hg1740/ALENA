function F06Data = f06extract(filename)
%f06extract Extracts the data from a .f06 file into a Matlab structure.

F06Data = [];

%Prompt user if no file is provided
if nargin == 0
    filename = getfile({'*.f06', 'F06 file'}, 'Select a file');
    if isempty(filename), return; end %Escape route
end

%Check extension
[~, ~, ext] = fileparts(filename);
assert(strcmp(ext, '.f06'), '''Filename'' must be the name of a .f06 file');

%Read every line from the file (assume all in one block) & discard comments
fid  = fopen(filename, 'r');
data = textscan(fid, '%s', 'delimiter', '\n','CommentStyle', '$', 'whitespace', '');
fclose(fid);
f06Data = data{1};
clear data

%Remove preamble
ind = findstring(f06Data, 'IFP OPERATIONS COMPLETE');
f06Data = f06Data(ind + 1 : end);

%Remove empty lines
f06Data(cellfun(@(x) isempty(x), strtrim(f06Data))) = [];

%Find the page index numbers & split into pages
pageIndex = findstring(f06Data, 'PAGE');
lb = [1 ; pageIndex(1 : end - 1) + 1];
ub = pageIndex;
pageData = arrayfun(@(i) f06Data(lb(i) : ub(i)), 1 : numel(lb), 'Unif', false);

%Grab the format data
FormatData = defineFormatData;

%Extract the aerodynamic forces
[F06Data.AeroForce, pageData]    = readOutputData(pageData, FormatData.AeroForce);
[F06Data.AeroPressure, pageData] = readOutputData(pageData, FormatData.AeroPressure);
[F06Data.VehicleCoeff, ~]        = readOutputData(pageData, FormatData.VehicleCoeff);
%charles 
[F06Data.Bendingmoment, pageData] = readOutputData(pageData, FormatData.Bendingmoment);

%Bespoke extraction functions
[F06Data.TrimAngle, pageData]      = extractTrimAngle(pageData, FormatData.TrimAngle);

% comment out for fixed aoa analysis 
% [F06Data.StabilityDeriv, pageData] = extractDerivatives(pageData, FormatData.StabilityDeriv);


%charles commented out following
%More formatting
% if ~isempty(F06Data.VehicleCoeff)
%     F06Data.VehicleCoeff = formatVehicleCoeff(F06Data.VehicleCoeff);
% end
    
end

%% Reading data from the f06

function FormatData = defineFormatData
%defineFormatData : Defines the format data for the output quantities in
%the .f06 file
%
% Detailed Explanation :
%   - Below is a description of each field of the output quantity:
%       + 'target'  : String that identifies the output quantity
%       + 'header'  : Header lines immediately above the output data
%       + 'format'  : Format string for sscanf to extract the numeric data
%       + 'colName' : Fieldname for each column of the output data
%   - The following output quantities are available:
%       + Aerodynamic Forces

%Aerodynamic Forces
AeroForce.target     = 'AERODYNAMIC FORCES ON THE AERODYNAMIC ELEMENTS';
AeroForce.header{1}  = 'GROUP  GRID ID  LABEL        T1                T2                T3                R1                R2                R3';
AeroForce.format     = '%i %i %*4c %f %f %f %f %f %f';
AeroForce.colName    = {'GROUP', 'PANELID', 'T1', 'T2', 'T3', 'R1', 'R2', 'R3'};
FormatData.AeroForce = AeroForce;

%Aerodynamic Pressures
AeroPressure.target     = 'AERODYNAMIC PRESSURES ON THE AERODYNAMIC ELEMENTS';
AeroPressure.header{1}  = 'GRID   LABEL          COEFFICIENTS           PRESSURES';
AeroPressure.format     = '%i %*2c %f %f';
AeroPressure.colName    = {'AeroGrid', 'PressCoeff', 'Press'};
FormatData.AeroPressure = AeroPressure;

%Vehicle Coefficients
VehicleCoeff.target     = 'AERODYNAMIC MONITOR POINT TOTAL VEHICLE COEFFICIENTS';
VehicleCoeff.header{1}  = 'AXIS         RIGID AIR  +  RESTRAINED INCR.  -  INERTIAL  +   RIGID-APPLIED + RESTRAINED INCR.  =  BALANCE';
VehicleCoeff.header{2}  = '----         ---------     ----------------     --------      -------------   ----------------     -------';
VehicleCoeff.format     = '%*10c %*10c %*f %*f %*10c %*10c %*f %f';
VehicleCoeff.colName    = {'Balance'}; %'RigidAir', 'RestrainedIncr', 'RestrainedInc2',
FormatData.VehicleCoeff = VehicleCoeff;

%Trim Angles
TrimAngle.target     = 'TRANSFORMATION FROM REFERENCE TO WIND AXES';
TrimAngle.header{1}  = 'TRANSFORMATION FROM REFERENCE TO WIND AXES';
TrimAngle.format     = '%*20c %f %*9c %f';
TrimAngle.colName    = {'Radians', 'Degrees'};
FormatData.TrimAngle = TrimAngle;

%Stability Derivatives
StabilityDeriv.target     = 'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S';
StabilityDeriv.header{1}  = 'TRIM VARIABLE   COEFFICIENT              RIGID                         ELASTIC                          INERTIAL';
StabilityDeriv.header{2}  = 'UNSPLINED        SPLINED       RESTRAINED      UNRESTRAINED     RESTRAINED      UNRESTRAINED';
StabilityDeriv.format     = '%*6c %f %f %f %f %f %f';
StabilityDeriv.colName    = {'RigidUnSplined', 'RigidSplines', 'ElasticRestrained', 'ElasticUnrestrained', 'InertialRestrained', 'InertialUnrestrained'};
FormatData.StabilityDeriv = StabilityDeriv;

%charles Bending moment
Bendingmoment.target = 'F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )';
Bendingmoment.header{1}  = 'STAT DIST/   - BENDING MOMENTS -            - WEB  SHEARS -           AXIAL          TOTAL          WARPING';
Bendingmoment.header{2}  = 'ELEMENT-ID  GRID   LENGTH    PLANE 1       PLANE 2        PLANE 1       PLANE 2        FORCE          TORQUE         TORQUE';
Bendingmoment.format = '%i %i %f %f %f %f %f %f %f %f %i %f %f %f %f %f %f %f %f';
Bendingmoment.colName =  {'ELEMENTID', 'UGRID','ULENGTH','UMPLN1','UMPLN2','USPLN1','USPLN2','UFORCE','UTORQUE1','UTORQUE2'...
    'LGRID','LLENGTH','LMPLN1','LMPLN2','LSPLN1','LSPLN2','LFORCE','LTORQUE1','LTORQUE2'};
FormatData.Bendingmoment=Bendingmoment;



end

function idx = bContains(data, str)
%bContains Returns a logical index indicating which elements of the
%cell-str array 'data' contain the text 'str'.
%
%   - TODO : Update once confirmed for Matlab 2017 onwards

idx = cellfun(@(x) ~isempty(strfind(x, str)), data);

end

function index = findstring(data, str)
%findstring Returns the index of the cell-str array 'data' which contains
%the text 'str'.

index = find(bContains(data, str));

end

function [ResultsData, pageData] = readOutputData(pageData, DataFormat)
%readOutputData : Reads the output data from the cell array f06Data
%
% Detailed Explanation :
%   - It is assumed that the 'format' field in 'DataFormat' is correct for
%     each set of results
%   - Only numeric data will be extracted using 'sscanf'. The format string
%     should be defined so that characters are skipped
%

ResultsData = [];

%Find the page with the data on it
idx      = cellfun(@(x) any(bContains(x, DataFormat.target)), pageData);
data     = pageData(idx);
pageData = pageData(~idx);

if isempty(data) %Escape route
    return
end

nSubcase = numel(data);

%Preallocate
fNames      = [DataFormat.colName, {'Subcase'}];
temp        = cell2struct(cell(size(fNames)), fNames, 2);
ResultsData = repmat(temp, [1, nSubcase]);
clear temp

%Assume each entry is grouped by subcase 
%   - Only valid if the input file specified enough lines in the Case
%   Control command "LINES = ..." to force all subcase data to be on a
%   single page.
for ii = 1 : numel(data)
   ResultsData(ii) = i_extractResultsData(data{ii}, DataFormat); 
end

    function ResultsData = i_extractResultsData(data, DataFormat)
        
        %Extract the subcase number
        subcaseLine = data(bContains(data, 'SUBCASE'));
        if isempty(subcaseLine)
            subcaseID = nan;
        else
            ind       = strfind(subcaseLine{1}, 'SUBCASE');
            subcaseID = sscanf(subcaseLine{1}(ind: end), 'SUBCASE %f');
        end
        
        %Find the header line that sits immediately above the numerical data
        ind = findstring(data, DataFormat.header{1});
        
        %Grab the table data
        tableData = data(ind(end) + numel(DataFormat.header) : end);
        
        %Combine into one long string for scanning
        str = cat(2, tableData{:});
        
        %Extract data into a column vector (vectorised part)
        [numericData, nData] = sscanf(str(2:end), DataFormat.format);
        nCol = numel(DataFormat.colName);
        
        %Reshape into a [nCol, nData] matrix
        sortedData = reshape(numericData, [nCol, nData/nCol]);
        
        %Convert to a structure
        ResultsData = cell2struct( ...
            [num2cell(sortedData, 2) ; {subcaseID}], ...
            [DataFormat.colName' ; {'Subcase'}]);
        
    end

end

%% Bespoke formatting

function VehicleCoeff = formatVehicleCoeff(VehicleCoeffRaw)
%formatVehicleCoeff Labels each of the vehicle coefficients so that the
%data makes more sense.

%Group coefficients into body and wind reference frames
balance = vertcat(VehicleCoeffRaw.Balance);
body    = balance(:, 1 : 2 : end);
wind    = balance(:, 2 : 2 : end);

labels = { ... 
    'CX', 'CY', 'CZ', 'CMX'    , 'CMY'     , 'CMZ'    ; ...
    'CD', 'CY', 'CL', 'CM_ROLL', 'CM_PITCH', 'CM_YAW'};

VehicleCoeff.Body    = cell2struct(num2cell(body, 1), labels(1, :), 2);
VehicleCoeff.Wind    = cell2struct(num2cell(wind, 1), labels(2, :), 2);
VehicleCoeff.Subcase = [VehicleCoeffRaw.Subcase];

end

%% Bespoke extraction 

function [TrimAngle, pageData] = extractTrimAngle(pageData, DataFormat)
%extractTrimAngle Extracts the trim angle (AoA & side-slip) and wind to
%reference axis rotation matrix for each subcase.

TrimAngle = [];

%Find the page with the data on it
idx      = cellfun(@(x) any(bContains(x, DataFormat.target)), pageData);
data     = pageData(idx);
pageData = pageData(~idx);

if isempty(data) %Escape route
    return
end

nSubcase = numel(data);

%Preallocate
TrimAngle = repmat(struct('AoA', [], 'SideSlip', [], 'RMatrix', []), [1, nSubcase]); 
aoa       = cell(1, nSubcase);
sideslip  = cell(1, nSubcase);
subcaseID = cell(1, nSubcase);
rMat      = cell(1, nSubcase);

%Extract the data
for i = 1 : nSubcase
     [aoa{i}, sideslip{i}, rMat{i}, subcaseID{i}] = i_extractTrimData(data{i}, DataFormat);   
end

%Replace empties with nan
aoa(cellfun(@isempty, aoa))           = {nan};
sideslip(cellfun(@isempty, sideslip)) = {nan};

%Convert to degrees
aoa      = cellfun(@rad2deg, aoa, 'Unif', false); 
sideslip = cellfun(@rad2deg, sideslip, 'Unif', false); 

%Assign to structure
[TrimAngle.AoA]      = deal(aoa{:});
[TrimAngle.SideSlip] = deal(sideslip{:});
[TrimAngle.RMatrix]  = deal(rMat{:});

    function [aoa, sideslip, rMat, sc] = i_extractTrimData(data, DataFormat)
        %i_extractTrimData Extracts the trim data and subcase number
        
        sc = [];
        
        %Extract the subcase number
        subcaseLine = data(bContains(data, 'SUBCASE'));
        if ~isempty(subcaseLine)
            ind = strfind(subcaseLine{1}, 'SUBCASE');
            sc  = sscanf(subcaseLine{1}(ind: end), 'SUBCASE %f');
        end
        
        %Find the header line that sits immediately above the numerical data
        ind = findstring(data, DataFormat.header{1});
        
        %Extract exactly the number of lines we need
        tableData = strtrim(data(ind + 1 : ind + 5));
        
        %Scan data
        aoa      = sscanf(tableData{1}, 'ANGLE OF ATTACK   = %f');
        sideslip = sscanf(tableData{2}, 'ANGLE OF SIDESLIP = %f');
        x = sscanf(tableData{3}, '%*16c %f %f %f');
        y = sscanf(tableData{4}, '%*16c %f %f %f');
        z = sscanf(tableData{5}, '%*16c %f %f %f');
        rMat = [x, y, z]';
        
    end
end

function [StabilityDeriv, pageData] = extractDerivatives(pageData, DataFormat)
%extractDerivatives Extracts the stabilitiy derivatives for each subcase.

StabilityDeriv = [];

%Find the page with the data on it
idx      = cellfun(@(x) any(bContains(x, DataFormat.target)), pageData);
data     = pageData(idx);
pageData = pageData(~idx);

if isempty(data) %Escape route
    return
end

nSubcase = numel(data);

for i = 1 : nSubcase
    
    %Extract the subcase number
    subcaseLine = data{i}(bContains(data{i}, 'SUBCASE'));
    if ~isempty(subcaseLine)
        ind       = strfind(subcaseLine{1}, 'SUBCASE');
        subcaseID = sscanf(subcaseLine{1}(ind: end), 'SUBCASE %f');
    end
    
    %Find the header line that sits immediately above the numerical data
    ind = findstring(data{i}, DataFormat.header{1});
    
    %Extract exactly the number of lines we need
    tableData = data{i}(ind + numel(DataFormat.header) : end);
    
    %Strip trim variables and stability data
    trimVarData = strtrim(cellfun(@(x) x(1 : 20), tableData, 'Unif', false));
    tableData   = cellfun(@(x) x(21 : end), tableData, 'Unif', false);
    
    %Grab trim variables     
    trimVarData = trimVarData(and( ...
        isnan(str2double(trimVarData)),  ...
        ~cellfun(@isempty, trimVarData)));
    
    %Make into valid fieldnames
    trimVarData(bContains(trimVarData, 'REF. COEFF.')) = {'RefCoeff'};    
        
    if i == 1
        fNames = [trimVarData ; {'Subcase'}];
        StabilityDeriv = repmat(cell2struct(cell(size(fNames)), fNames), [1, nSubcase]);
    end
    
    StabilityDeriv(i).Subcase  = subcaseID;
    
    %Set bounds for indexing
    nTrimVar = numel(trimVarData);
    ub       = cumsum(repmat(6, [1, nTrimVar]));
    lb       = [1, ub(1 : end - 1) + 1];
    
    %Extract the stability derivatives
    for j = 1 : nTrimVar
        
        trimVar = trimVarData{j};
        nCol    = numel(DataFormat.colName);
        
        %Extract data using 'sscanf'
        stabStrData = tableData(lb(j) : ub(j));
        str         = cat(2, stabStrData{:});
        stabData    = sscanf(str, DataFormat.format);
        
        %Reshape into a [nCol, nData] matrix
        sortedData = reshape(stabData, [nCol, numel(stabData)/nCol]);
        
        %Return to a structure
        StabilityDeriv(i).(trimVar) = cell2struct( ...
            num2cell(sortedData, 1), ...
            {'CX', 'CY', 'CZ', 'CMX', 'CMY', 'CMZ'}, 2);
        
        
    end

end

end
