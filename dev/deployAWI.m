function deployAWI(deployTarget, thingsToDeploy)
%DeployDevVersion Deploys a user-appropriate version of the code.
%
% A couple of things to note:
%   1. The code will be published using p-files.
%   2. It is asssumed that this function is saved in a Git repository.
%      Therefore, when deploying the code base the new code is placed in
%      the parent directory of the Git repository.

if nargin < 1
    [~, deployTarget, ~] = fileparts(pwd);
end
if nargin < 2
    thingsToDeploy = { ...
        fullfile(pwd, '+awi'), ...
        fullfile(pwd, '+mvc'), ...
        fullfile(pwd, 'UoB') , ...
        fullfile(pwd, 'Models'), ...
        fullfile(pwd, 'ALENA.m')};
end

%SET UP
currentDir = pwd;
addpath(genpath(currentDir));

%Get path to the parent directory of the Git Repo
parentDir = navigateToParentDir(pwd);

%Create a new folder - Rename after this folder
tf           = mkdir(parentDir, deployTarget);
deployTarget = fullfile(parentDir, deployTarget);

if ~tf %Tell user if we fail to make the new folder
    error(['Failed to create a new folder in ''%s'' for the deployed ', ...
        'code.'], deployTarget);
end

for iT = 1 : numel(thingsToDeploy) %Attempt to p-code the item
    obfuscateTheThing(thingsToDeploy{iT}, deployTarget)
end

%TEAR DOWN
cd(currentDir);     %Return to original workspace
rmpath(currentDir); %Remove from the path

end

function parentDir = navigateToParentDir(currentDir)
%navigateToParentDir Finds the parent directory of the .Git repo

%Get current directory data
dirData  = dir(currentDir);

%Want folders only (including hidden)
dirNames = {dirData([dirData.isdir]).name};

%Any folders containg .git?
idx = contains(dirNames, '.git');

if any(idx)
    %Go one level up - Recurse function
    [parentDir, ~, ~] = fileparts(currentDir);
else
    %Recurse function
    [parentPath, ~, ~] = fileparts(currentDir);
    parentDir = navigateToParentDir(parentPath);
end

end

function obfuscateTheThing(thing, parentDir)
%obfuscateTheThing Provides the functionality for p-coding the .m files.

%Change directory so that the p-files end up in this folder
cd(parentDir);

%What have we got?
[path, subDir, ext] = fileparts(thing);

if isdir(thing)
    %Make a new directory so we preserve the file structure
    newDir = fullfile(parentDir, subDir);
    tf     = mkdir(parentDir, subDir);
    if ~tf %Tell user if we fail to make a folder
        error('Failed to make folder ''%s'' in the deployed directory', ...
            newDir);
    end
    %What is in the directory?
    contents = dir(thing);
    subThings = {contents.name};
    subThings(ismember(subThings, {'.', '..'})) = [];
    %Recurse
    for iST = 1 : numel(subThings)
        theThing = fullfile(path, subDir, subThings{iST});
        obfuscateTheThing(theThing, newDir);
    end
elseif strcmp(ext, '.m')
    %Turn thing into p-code
    pcode(thing)
else
    %Copy the thing
    copyfile(thing, pwd);
end

end
