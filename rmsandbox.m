function rmsandbox()
%rmsandbox  Uninstall sandbox
%
%  See also: rmsandbox, 
%            modify_sandbox_path
%
%  Copyright 2016 The MathWorks, Inc.
%
% Edits by Christopher Szczyglowski, University of Bristol, 2020
%   - Refactored most of the code into 'modify_sandbox_path'

sub_directory_to_remove = {'tbx' ; 'tests' ; 'modules' ; 'examples' ; 'dev'};

modify_sandbox_path(sub_directory_to_remove, 'remove');

% call submodules
run(fullfile('modules','Matran','rm_sandbox.m'))

end