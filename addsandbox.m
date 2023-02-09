function addsandbox()
%addsandbox  Install sandbox
%
%  See also: rmsandbox, 
%            modify_sandbox_path
%
%  Copyright 2016 The MathWorks, Inc.
%
% Edits by Christopher Szczyglowski, University of Bristol, 2020
%   - Refactored most of the code into 'modify_sandbox_path'

sub_directory_to_add = {'tbx' ; 'tests' ; 'modules' ; 'examples' ; 'dev'};

modify_sandbox_path(sub_directory_to_add, 'add');

% call submodules
run(fullfile('modules','Matran','add_sandbox.m'))

end 