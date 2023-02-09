%% update 05/18 - 06/18
%
%   This file details the major updates that have been made to the AWI
%   Framework codebase in the period May 2018 - June 2018.
%
%   Author(s): Chris Szczyglowski, Dario Calderon
%
%
%   1. Added methods for updating the 'Area Set' (aSet). This calculates
%      the (projected) area, taper ratio and aspect ratio of each section
%      along a lifting surface.
%
%   2. Significantly updated the "Beamable" and "BeamProperty" classes. 
%
%       a. If a "Beamable" object is also of class "awi.model.Stick" then
%       additional listeners are added so that each beam property stores
%       the 'XData', 'YData' & 'ZData' of the beam and listens to any
%       changes. This allows different beam properties to be defined over 
%       the global axes (X,Y,Z) as well as along the beam axis (R) and the
%       interpolation will automatically combine and handle the multiple
%       axis definitions and eta positions.
%
%       b. If a "Beamable" is of class "awi.model.LiftingSurface" then all
%       beam properties that are added by the LiftingSurface constructor
%       will automatically have their '_flag' variable set to the
%       SpanVector of the lifting surface.
%   
%       c. New methods and properties were added to "Beamable" and
%       "BeamProperty" to handle this change.
%
%       c. This now means that for all beam properties that are added as
%       dynamic properties during the "addBeamProperty" method the '_flag'
%       property will default to 'R' (along the beam axis) which is
%       contrary to much of the default logic of LiftingSurfaces.
%       Therefore, it is best practice for all objects inheriting from
%       "Beamable" to explicitly define the '_flag' property in the input
%       file. Or in the case of FAME files, during the conversion script.
%
%   3. Added new mixin class "FEable". This class handles converting a 
%      'geometry model' to an 'analysis model'. For now, this relies on a 
%      finite element logic to discretise the model into constituent parts.
%      The logic is that all analyses will go through this block in the
%      future.

%% Wish list
%
%   1. Allow objects to be added to the collection as children but hidden
%      from the tree view. This is required because sometimes additional
%      objects are added to the collection during the build process (e.g.
%      CoordSys for LiftingSurfaces) and this causes the build process to
%      slow down significantly. Furthermore, as these are 'utility' objects
%      we don't really want the use to see them and interact with them.
%      These objects are basically dependent and will be removed/added each
%      time the model is built so want to minimise disruption. 
%       a. Would it be possible to define a new type of collector for these
%       objects? The collector would be a leaf node (in the sense of a tree
%       view) but would still have child objects. In the tree view it would
%       display the number of children. For the case of point masses (of
%       which there could be 100s/1000s) we would have "Point Masses(1345)"
%       where 1345 is the number of point masses in that collection.
%
%   2. Add groups to the layered drawing. E.g. Geometry and Beam/FE
%      properties. Could have the option to toggle all members of a group
%      at once.
%
%   3. Ability to relate properties to other properties in the xml file.
%      i.e. treat properties as variables.
%
%   4. "FAME_FEM" is appearing in the results set as it is of type
%      "BeamModel" which really should be moved as a peer of the "Aircraft"
%      & LoadCases containers in the uitree. Action: Make the BeamModel an
%      Entity instead of a results set.
%
%   5. Need the ability to multi-select trim load cases to run in the trim
%      analysis.
%
%   6. All results sets (e.g. TrimResult, SizeResults, etc.) should inherit
%      from a common superclass. To do this, the analyses will need to 
%      output to a commn superclass in the first place!
%

