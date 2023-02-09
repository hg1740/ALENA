function release()
%RELEASE An example release function to build, test and package a MATLAB
%toolbox. The release version is assumed to be up to date in the
%'Contents.m' file. The release in the prj file is automatically updated. 
%
%  Copyright 2016-2017 The MathWorks, Inc.

%% Set toolbox name
tbxname = get_toolbox_name;

%% Get release script directory
cfdir = fileparts( mfilename( 'fullpath' ) );
tbxDir = fullfile( cfdir, 'tbx');

%% Add sandbox to MATLAB path
run( fullfile( cfdir, 'addsandbox.m' ) );

%% Check MATLAB and related tools, e.g.:
assert( ~verLessThan( 'MATLAB', '9.5' ), 'MATLAB R2018b or higher is required.' )

%% Check installation
fprintf( 1, 'Checking installation...' );
v = ver( tbxname );
switch numel( v )
    case 0
        fprintf( 1, ' failed.\n' );
        error( '%s not found.', tbxname );
    case 1
        % OK so far
        fprintf( 1, ' Done.\n' );
    otherwise
        fprintf( 1, ' failed.\n' );
        error( 'There are multiple copies of ''%s'' on the MATLAB path.', tbxname );
end

%% Build documentation & examples
fprintf( 1, 'Generating documentation & examples...' );
try
    % Do something;
    fprintf( 1, ' Done.\n' );
catch e
    fprintf( 1, ' failed.\n' );
    e.rethrow()
end

%Build doc search database
builddocsearchdb( fullfile( tbxDir, 'doc' ));

%% Run tests
fprintf( 1, 'Running tests...' );
tdir = fullfile( fileparts( tbxDir ), 'tests' );
[log, results] = evalc( sprintf( 'runtests( ''%s'' )', tdir ));
if ~any( [results.Failed] )
    fprintf( 1, ' Done.\n' );
else
    fprintf( 1, ' failed.\n' );
    error( '%s', log )
end

%%  Package and rename.
fprintf( 1, 'Packaging...' );
try
    % Update prj file version to match the version in the Contents file
    prj = fullfile( cfdir, [ tbxname, '.prj'] );
    matlab.addons.toolbox.toolboxVersion(prj,v.Version);
    
    if verLessThan( 'matlab', '9.5' )
        % Undocumented API.
        name = char( com.mathworks.toolbox_packaging.services.ToolboxPackagingService.openProject( prj ) );
        com.mathworks.toolbox_packaging.services.ToolboxPackagingService.packageProject( name )
        com.mathworks.toolbox_packaging.services.ToolboxPackagingService.closeProject( name )
    else
        matlab.addons.toolbox.packageToolbox( prj );
    end
    oldMltbx = which( [tbxname '.mltbx'] );
    newMltbx = fullfile( fileparts( tbxDir ), 'releases', [tbxname ' v' v.Version '.mltbx'] );
    movefile( oldMltbx, newMltbx )
    fprintf( 1, ' Done.\n' );
catch e
    fprintf( 1, ' failed.\n' );
    e.rethrow()
end

%% Check package
fprintf( 1, 'Checking package...' );
if verLessThan( 'matlab', '9.5' )
    info = mlAddonGetProperties( newMltbx );
    tver = info.version;
else
    tver = matlab.addons.toolbox.toolboxVersion( newMltbx );
end 

if strcmp( tver, v.Version )
    fprintf( 1, ' Done.\n' );
else
    fprintf( 1, ' failed.\n' );
    error( 'Package version ''%s'' does not match code version ''%s''.', tver, v.Version )
end

%% Show message
fprintf( 1, 'Created package ''%s''.\n', newMltbx );

end %release
