%% DEPLOYAEROFLEX:
% The following function takes in p = profile('info') and copies the files
% listed in the profile over to a specifed folder.
% 
%   Author: Dario Calderon (University of Bristol 2018)

function deployAeroFlex(p,DestFolder,Str)

% Allow the user not to specify the destination folder
if nargin == 1
    error('Destination folder required');
end

if nargin == 2
    NoStr = 1;
end

if nargin == 3
    NoStr = 0;
end

% Extract the list of functions identifed by the profiler
Sizing_f   = {p.FunctionTable.FileName};

% Isolate the functions which contain the string AeroFlex
if NoStr == 0
    AeroFlex_f = find(~cellfun(@isempty,regexp(Sizing_f,Str)));
else
    AeroFlex_f = 1:numel(Sizing_f);
end

% Create the destination folder
if ~exist(DestFolder,'dir')
    mkdir(DestFolder);
end

% Copy the files across
fprintf('\n\t...Copying files ...');
for i = 1:numel(AeroFlex_f)
    path = p.FunctionTable(AeroFlex_f(i)).FileName;
    [~,fname,ext] = fileparts(path);
    copyfile(path,[DestFolder filesep fname ext]);
end
fprintf('\n\t...Copy complete.\n');

end