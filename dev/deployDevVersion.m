function deployDevVersion(deployTarget, thingsToDeploy)
%DeployDevVersion Deploys a development version of the code which can be
%used for sharing with other collaborators. 
%
% A couple of things to note:
%   1. The code will be published using m-files.
%   2. All the directories and files in 'thingsToDeploy' will be copied to
%   'deployTarget'.
%   3. If 'deployTarget' is not provided then the user will be prompted for
%   a directory to copy the files/folders to.

if nargin < 1 %Ask the user for a directory to copy the files to
    deployTarget = uigetdir(pwd, 'Select a folder location for the AWI Framework code base to be deployed to.');
end
if nargin < 2 %Maintain a manual list of directories and files to copy
    thingsToDeploy = { ...
        fullfile(pwd, '+awi'), ...
        fullfile(pwd, '+mvc'), ...
        fullfile(pwd, 'UoB') , ...
        fullfile(pwd, 'Models'), ...
        fullfile(pwd, 'ALENA.m') , ...
        fullfile('\\ads\filestore\Engineering\Research\Projects\AWI\Analysis_Tools\UoB_Framework\GLT')};
end


%How many things are we deploying?
nThings = numel(thingsToDeploy);

hWB = waitbar(0, 'Copying files ...');

for iDir = 1 : nThings %Copy all items across
    
    %What is the name of the subdirectory?
    [~, subDir, ext] = fileparts(thingsToDeploy{iDir});
    
    %Define full path of the subdirectory
    subDirPath = fullfile(deployTarget, [subDir, ext]);
        
    %Copy the contents into the new subdirectory
    tf = copyfile(thingsToDeploy{iDir}, subDirPath);
        
    if ~tf %Tell user if we fail to copy a folder
        
        %Delete the waitbar
%         delete(hWB);
        
        %Tell the user what happened
        warning('Failed to copy folder ''%s'' to the deployed directory', ...
            thingsToDeploy{iDir});
        
    end
   
    %Update the waitbar
    waitbar((iDir / nThings), hWB, sprintf( ...
        'Successfully copied file/directory %i/%i', iDir, nThings));
    
end

%Get rid of the waitbar
delete(hWB);

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
