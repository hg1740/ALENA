function publish_code(Options)
%publish_code Publishes the ALENA documentation using the comments in the
%code files. 
%
% Detailed Description:
%
%
%

if nargin < 1
    Options = struct( ...
        'Doc_Destination', pwd, ...
        'Doc_Format'     , {{'html'}}, ...
        'TopLevelPackage', 'awi');      %Must be on the path
end

mp = meta.package.fromName(Options.TopLevelPackage);
if isempty(mp)
    error(['Unable to generate documentation as the ''%s'' package ', ...
        'could not be found on the path.']);
end

publish_package(mp, Options);

end

%% Publish packages

function publish_package(mp, Options)
%publish_package

valid_format = {'html', 'latex'};

validateattributes(mp, {'meta.package'}, {'scalar'}, 'mp', 'publish_package');
validateattributes(Options.Doc_Format, {'cell'}, {'vector'}, 'format', 'publish_package');
i_validate_format(Options.Doc_Format, valid_format);

    function i_validate_format(format, valid_format)
        
        assert(iscellstr(format), ['Expected ''format'' to be a cell-array ', ...
            'of strings.']);
        assert(all(contains(format, valid_format)), sprintf(['Expected ', ...
            '''format'' to contain the following tokens\n\t%s\n\n']     , ...
            strjoin(valid_format, ', ')));
        
    end

%Create a summary of the package 

%Differentiate by package


%Publish by file-type
for ii = 1 : numel(mp.ClassList) %Classes
    publish_class(mp.ClassList(ii));    
end
for ii = 1 : numel(mp.FunctionList) %Functions 
    
end
for ii = 1 : numel(mp.PackageList) %Packages
    publish_package(mp.PackageList(ii), Options)
end

%Differentiate by file type (script, function, class)

%

end

function publish_class(mc, Options)

%fullfilepath = which(mc.Name);

h1 = help(mc.Name);

%Interpret markup

%Write to file
switch Options.Format
    
end
fid = fopen(fullfile(pwd, 'MyTestFile.txt'), 'w');
fprintf(fid, '%s\n', h1);
fclose(fid);

end

function MarkupMap = getMarkupMapping

MarkupMap = cell( ...
    'reStructuredText', 'LaTeX', 'HTML' ; ...
    '* *'
end

