%%GETFAMEVLM: A function that reads a FAME output VLM file and stores the
%%results into a matlab structure.

function results = getFAMEVLM(folderName,PlotMesh)

% Deal with the inputs
if nargin > 0
    pathName = [folderName filesep,'__mon',filesep,'l3_fame-w',filesep,'aeromeshes',filesep,'VL_panels'];
    fileName = 'fame-w_VL_panels__report.txt';
else
    [fileName,pathName,~] = uigetfile('*.txt','Pick the VLM file');
    fprintf('\n%s',[pathName,fileName]);
end

if nargin < 2
    PlotMesh = 0 ;
end

% Read the file
blk = textread(fullfile(pathName,fileName),'%s','delimiter','\n');

% Locate the different load cases
idx    = find(~cellfun(@isempty,regexp(blk,'Load case no.')));

for i = 1:numel(idx)-1
    
    if i == numel(idx)-1
        bulk{i} = blk(idx(i+1)+7:numel(blk)-2);
    else
        bulk{i} = blk(idx(i+1)+7:idx(i+2)-5);
    end
    
    tab = char(bulk{i});
    tab = str2num(tab);
    
    % Extract the combined cs deflection
    CSDeflection = sum(tab(:,32:38),2);
    
    % Extract the camber
    Camber = tab(:,28);
    
    % Extract the unit normal of the panel
    N      = tab(:,24:26);
    
    XYZ(:,:,1) = [tab(:,5),tab(:,8),tab(:,11),tab(:,14)]/1000;
    XYZ(:,:,2) = [tab(:,6),tab(:,9),tab(:,12),tab(:,15)]/1000;
    XYZ(:,:,3) = [tab(:,7),tab(:,10),tab(:,13),tab(:,16)]/1000;
    
    % Store the results into a structure
    results(i).XYZ           = XYZ;
    results(i).N             = N;
    results(i).Camber        = Camber;
    results(i).CSDelflection = CSDeflection;
    
    % Plot the mesh if requested by the user
    if PlotMesh == 1
        figure;
        
        fill3(XYZ(:,:,1)', XYZ(:,:,2)', ...
            XYZ(:,:,3)', CSDeflection,'FaceColor','flat')
        
        title(['LoadCase #' num2str(i)]);
        axis equal;
        colorbar
        colormap(jet);
        drawnow;
    end
    
end

end