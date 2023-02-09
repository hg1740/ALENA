%% Development notes 16/06/12018


% 1. We want to be able to import from CPACS files
%   - The problem is, the files are massive!
%   - Need a bespoke xml import routine.
%   - Make a list of the tokens we care about:
%       + (YES) cpacs>vehicles
%       + (YES) cpacs>vehicles>aircraft
%       + (YES) cpacs>vehicles>profiles  (Grab everything)
%       + (YES) cpacs>vehicles>materials (Grab everything)
%       + (YES) cpacs>vehicles>engines
%       + (NO)  cpacs>airports
%       + (NO)  cpacs>missions
%
%
%
% 2. The log function in UITools is spitting a logical out to the command
% window under certain circumstances. See the logfcn call after adding 
% control surfaces in the "import_fm4" method.  (logfcn('Control surface
% deflections assigned!')
%
% 3. Update the name of the 'awi.model.PointMasses' collector so it
% reflects how much mass the children point masses contain.