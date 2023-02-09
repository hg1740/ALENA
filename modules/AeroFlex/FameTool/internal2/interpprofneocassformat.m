function obj = interpprofneocassformat(obj)
% Use this routine to interpolate profile data at specific eta values in
% NEOCASS format

% set eta values here.
% This should not be uniformly distributed!!!!!
%eta = 0:1/(ny-1):1;

% Extract the airfoils

FolderFlag = 0;

if isempty(obj.Dirs.airfoilFolder)
    warning('No airfoil folder exists in the "Dirs" object');
    %                 error('FAME2MAT:DirNotExist','Fame directory not defined. Please add a directory in the "Dirs" object');
    FolderFlag = 1;
    
elseif ~exist(obj.Dirs.airfoilFolder)
    
    warning('The airfoil in the "Dirs" object does not exist');
    FolderFlag = 1;
    
end

if FolderFlag == 1
    
    AIRFOIL_folder = uigetdir([],'Choose airfoil directory');
    set(obj.Dirs.airfoilFolder,AIRFOIL_folder);
    
end


prof = [];

% Filter through the airfoils found in the FAME input file
for i = 1:length(obj.Fame.fm4Input.WING.GLOBAL_GEOMETRY.PROFILE_DISTRIBUTION.NAME)
    
    prof(i).eta  = str2double(obj.Fame.fm4Input.WING.GLOBAL_GEOMETRY.PROFILE_DISTRIBUTION.ETA{i});
    
    check_var = strtrim(regexp(obj.Fame.fm4Input.WING.GLOBAL_GEOMETRY.PROFILE_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = str2double(obj.Fame.fm4Input.USER_VARIABLES.' check_var{1} '{1});']);
        prof(i).eta = Eta;
    end
    
    prof(i).name = obj.Fame.fm4Input.WING.GLOBAL_GEOMETRY.PROFILE_DISTRIBUTION.NAME{i};
    
    ProfileStructure = getFameInput([obj.Dirs.airfoilFolder filesep prof(i).name '.prf' ]);
    
    prof_x = [];
    prof_z = [];
    for j = 1:length(ProfileStructure.PROFILE.GEOMETRY.COORDINATES.X_REL)
        prof_x = [prof_x;str2double(ProfileStructure.PROFILE.GEOMETRY.COORDINATES.X_REL{j})];
        prof_z = [prof_z;str2double(ProfileStructure.PROFILE.GEOMETRY.COORDINATES.Z_REL{j})];
    end
    
    prof(i).coords.x = prof_x;
    prof(i).coords.z = prof_z;
    
end

etavec = [prof.eta];

plot_airfoil = 0;

% Determine where the airfoils need to be interpolated
WingCaero  = findobj(obj.Mdl.Caero,'part','Wing');
eta1 = sort(unique([WingCaero.startY])/max([WingCaero.startY] + [WingCaero.b]));
eta2 = sort(unique([WingCaero.startY]+[WingCaero.b])/max([WingCaero.startY] + [WingCaero.b]));

eta = uniquetol([eta1,eta2],1e-4);

% Perform an interpolation between the airfoils instead

x_coord = [];
z_coord = [];

for i = 1:numel(etavec)
    
    x_coord = [x_coord,prof(i).coords.x];
    z_coord = [z_coord,prof(i).coords.z];
    
    % Calculate the number of points that define the lower/upper surface
    numcoord = ceil(numel(prof(i).coords.x)/2);
    [prof(i).LScoords.x,c_order] = sort(prof(i).coords.x(1:numcoord),'ascend');
    
    prof(i).LScoords.z = prof(i).coords.z(1:numcoord);
    prof(i).LScoords.z = prof(i).LScoords.z(c_order);
    
    prof(i).UScoords.x = prof(i).coords.x(end-numcoord+1:end);
    prof(i).UScoords.z = prof(i).coords.z(end-numcoord+1:end);
end

x_interp = [];
z_interp = [];

for i = 1:size(x_coord,1)
    
    x_interp(i,:) = interp1(etavec,x_coord(i,:),eta);
    z_interp(i,:) = interp1(etavec,z_coord(i,:),eta);
    
end

for i = 1:size(x_interp,2)
   
    interpprof(i).eta = eta(i);
    interpprof(i).coords.x = x_interp(:,i);
    interpprof(i).coords.z = z_interp(:,i);
    
     % Calculate the number of points that define the lower/upper surface
    numcoord = ceil(numel(interpprof(i).coords.x)/2);
    [interpprof(i).LScoords.x,c_order] = sort(interpprof(i).coords.x(1:numcoord),'ascend');
    
    interpprof(i).LScoords.z = interpprof(i).coords.z(1:numcoord);
    interpprof(i).LScoords.z = interpprof(i).LScoords.z(c_order);
    
    interpprof(i).UScoords.x = interpprof(i).coords.x(end-numcoord+1:end);
    interpprof(i).UScoords.z = interpprof(i).coords.z(end-numcoord+1:end);
end

%% CHECK TO MAKE SURE THAT INTERPOLATION HAS WORKED

h_airfoil = figure;

h_prof    = hggroup;
h_interp  = hggroup;
for i = 1:numel(etavec)
    plot3(etavec(i)*ones(size(x_coord(:,i))),prof(i).coords.x,prof(i).coords.z,'ro','Parent',h_prof);
    hold on;
end

for i = 1:numel(eta)
    plot3(eta(i)*ones(size(x_interp(:,i))),interpprof(i).coords.x,interpprof(i).coords.z,'bo','Parent',h_interp);
    hold on;
end

legend([h_prof,h_interp],'Original Profiles','Interpolated Profiles','Location','Best');
xlabel('x/c'); ylabel('eta'); zlabel('z/c');

saveas(h_airfoil,[obj.Dirs.writeFolder filesep 'ModelFigures' filesep 'AirfoilInterpolation.fig']);

for i = 1:numel(obj.Mdl.Caero)
    
    if strcmp(obj.Mdl.Caero(i).part,'Wing')
        
        Caero_eta_in = obj.Mdl.Caero(i).startY/max([WingCaero.startY] + [WingCaero.b]);
        Caero_eta_out = (obj.Mdl.Caero(i).startY+obj.Mdl.Caero(i).b)/max([WingCaero.startY] + [WingCaero.b]);
        
        idx_in  = find(abs([interpprof.eta] - Caero_eta_in)<0.0001);
        idx_out = find(abs([interpprof.eta] - Caero_eta_out)<0.0001);
        
        obj.Mdl.Caero(i).rootAirfoil = ['fame' num2str(idx_in)];
        obj.Mdl.Caero(i).tipAirfoil  = ['fame' num2str(idx_out)];
    end
end


if ~exist([obj.Dirs.writeFolder filesep 'FAMEAirfoils'],'dir')
    mkdir([obj.Dirs.writeFolder filesep 'FAMEAirfoils']);
else
    rmdir([obj.Dirs.writeFolder filesep 'FAMEAirfoils'],'s');
    mkdir([obj.Dirs.writeFolder filesep 'FAMEAirfoils']);
end

for i = 1:numel(prof)
      
    fid = fopen([obj.Dirs.writeFolder filesep 'FAMEAirfoils' filesep prof(i).name '.dat'],'w');
    
    %fprintf(fid,'$NACA 0012 AIRFOIL WITH ADJUSTED CAMBER, ETA : %4.3f\n',interpprof(i).eta);
    
    fprintf(fid,'%i  %i \n',numcoord,numcoord);
     
%     c = {};
%     f = {};
%     
%     c{1,1} = {[prof(i).UScoords.x;prof(i).LScoords.x]};
%     c{2,1} = {[prof(i).UScoords.z;prof(i).LScoords.z]};
%     
%     f{1,1} = '%10.7f  %10.7f \n';
%     
%     dataToWrite = cat(2, c{:});
%     format      = cat(2, f{:});
%     fprintf(fid, format, dataToWrite{:});
    
    %fprintf(fid,'$ Upper Surface');
    %fprintf(fid,'$      x        z');
    for j = 1:numel(prof(i).UScoords.x)
        fprintf(fid,'\n%10.7f  %10.7f',prof(i).UScoords.x(j),prof(i).UScoords.z(j));
    end
    %fprintf(fid,'\n$ Lower Surface');
    %fprintf(fid,'$      x        z');
    for j = 1:numel(prof(i).LScoords.x)
        fprintf(fid,'\n%10.7f  %10.7f',prof(i).LScoords.x(j),prof(i).LScoords.z(j));
    end
    
    fclose(fid);
    
     % Extract the Camber information
    [~, ~, ~, ~, angle] = airfoil_mean_line({[obj.Dirs.writeFolder filesep 'FAMEAirfoils' filesep ,prof(i).name,'.dat']},0,1,0);
    
    prof(i).camber = angle';
end

obj.Fame.Geometry.Wing.Airfoils = prof;

if ~exist([obj.Dirs.writeFolder filesep 'NeoCASSAirfoils'],'dir')
    mkdir([obj.Dirs.writeFolder filesep 'NeoCASSAirfoils']);
else
    rmdir([obj.Dirs.writeFolder filesep 'NeoCASSAirfoils'],'s');
    mkdir([obj.Dirs.writeFolder filesep 'NeoCASSAirfoils']);
end
% 
% dataToWrite = cat(1, c{:});
% format      = cat(2, f{:});
% fprintf(fid, format, dataToWrite{:});

for i = 1:numel(interpprof)
    
    name = strrep(sprintf('%i',i),'.','');
    
    
    fid = fopen([obj.Dirs.writeFolder filesep 'NeoCASSAirfoils' filesep 'fame' name  '.dat'],'w');
    
    %fprintf(fid,'$NACA 0012 AIRFOIL WITH ADJUSTED CAMBER, ETA : %4.3f\n',interpprof(i).eta);
    
    fprintf(fid,'%i  %i \n',numcoord,numcoord);
    
    
%     c = {};
%     f = {};
%     
%     c{1,1} = {interpprof(i).UScoords.x};
%     c{2,1} = {interpprof(i).UScoords.z};
%     
%     f{1,1} = '%10.7f  %10.7f \n';
%     f{2,1} = '%10.7f  %10.7f \n';
%     
%     dataToWrite = cat(1, c{:});
%     format      = cat(2, f{:});
%     fprintf(fid, format, dataToWrite{:});
%     
    %fprintf(fid,'$Upper Surface');
    for j = 1:numel(interpprof(i).UScoords.x)
        fprintf(fid,'\n%10.7f  %10.7f',interpprof(i).UScoords.x(j),interpprof(i).UScoords.z(j));
    end
    %fprintf(fid,'\n$Lower Surface');
    for j = 1:numel(interpprof(i).LScoords.x)
        fprintf(fid,'\n%10.7f  %10.7f',interpprof(i).LScoords.x(j),interpprof(i).LScoords.z(j));
    end
    
    fclose(fid);
    
         % Extract the Camber information
    [~, ~, ~, ~, angle] = airfoil_mean_line({[obj.Dirs.writeFolder filesep 'NeoCASSAirfoils' filesep 'fame',name,'.dat']},0,1,0);
    
    interpprof(i).camber = angle';
end