%
% 1. Update properties of LiftingSurface, BLuffBody, etc. so that we can
% initiate a blank object and call 'get(obj)' on the object and not have
% the code fail.
%
% 2. Finish off the 'Beam' class and investigate a flexible way of
% implemented a beam that has properties.
%   * Can I implement a Structure called 'BeamProperties' structure or use
%   Dynamic properties to implement different beam forumulations?
%
% 3. Implement a 'view.BeamAnalysis' class that allows us to view any
% object of class Beamable and view the distribution of all or none of the
% beam properties across a specified axes...