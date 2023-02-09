%% DEPLOYAEROFLEX:
% The following function takes in p = profile('info') and copies the files
% listed in the profile over to a specifed folder.
% 
%   Author: Dario Calderon (University of Bristol 2018)

function deployAeroFlex(p,DestFolder)

% Allow the user not to specify the destination folder
if nargin < 2
    DestFolder = 'AeroFlexSizing';
end

% Extract the list of functions identifed by the profiler
Sizing_f   = {p.FunctionTable.FileName};

% Isolate the functions which contain the string AeroFlex
AeroFlex_f = find(~cellfun(@isempty,regexp(Sizing_f,'AeroFlex')));

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