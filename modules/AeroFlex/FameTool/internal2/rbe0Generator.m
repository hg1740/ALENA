%% RBE0GENERATOR   Create aero connection elements
%
%  RBE0GENERATOR(CLASSNAME) adds connection elements between aerodynamic 
%  and structure nodes for NEOCASS model
%
%  CLASSNAME is a Fame2mat object
%
%  This is a private function and should be accessed from the |Fame2mat|
%  object.
%
%   Example:
%      obj = fame2mat;
%      obj = getFameResults(obj);
%      obj = getFameStructure(obj);
%      obj = remeshStructure(obj);
%      obj = getFameAero(obj);
%      obj = setRbe0Nodes(obj);   % calls this function
%
%   Restrictions:
%       NOTE: FAME model must be loaded before this function can be used
%
%   See also |Rbe0| definition in NEOCASS manual

%   Copyright 2016 University of Bristol
%   Private function.
function obj = rbe0Generator(obj,part)

PartID = [part,'StartID'];
if strcmpi(part,'vtp')
    PartName = part;
else
    PartName = ['Stbd',part];
end

Node      = obj.Mdl.Grid;

NodeIdx   = obj.Opts.Struct.(PartID); 

Outline   = obj.Fame.Geometry.(part).Outline;
rbeOrient = obj.Opts.Struct.rbeOrient;

%% Create the leading / trailing edge representative vectors
% Count the number of kinks on the wing - assume the same number on LE and TE
nKinks    = size(Outline,1) / 2;
leKinks   = Outline(1:nKinks,:);
teKinks   = Outline(nKinks + 1:end,:);
nSections = size(leKinks,1) - 1;

if strcmpi(part,'vtp')
    [~,idx] = sort(leKinks(:,3),'ascend');
    leKinks = leKinks(idx,:);
    [~,idx] = sort(teKinks(:,3),'ascend');
    teKinks = teKinks(idx,:);
else
    [~,idx] = sort(leKinks(:,2),'ascend');
    leKinks = leKinks(idx,:);
    [~,idx] = sort(teKinks(:,2),'ascend');
    teKinks = teKinks(idx,:);
end

if leKinks ~= teKinks
    warning('The number of leading edge kinks do not match the trailing edge')
end

leVec = zeros(nSections,3);
teVec = zeros(nSections,3);

for i = 1:nSections
    leVec(i,1:3) = leKinks(i + 1,:) - leKinks(i,:);
    teVec(i,1:3) = teKinks(i + 1,:) - teKinks(i,:);
end

%% Create the beam representative vectors
% Assume that the beam tangent is equal to the vector connecting the
% inboard grid node
WingBeamNodes = findobj(Node,'part',PartName,'-and','type','Beam');

beamVec = zeros(numel(WingBeamNodes),3);

count = 0;

for i = 1:numel(WingBeamNodes)
    count = count + 1;
    if i == 1
        beamVec(count,1:3) = WingBeamNodes(i + 1).coord - WingBeamNodes(i).coord;
    else
        beamVec(count,1:3) = WingBeamNodes(i).coord - WingBeamNodes(i - 1).coord;
    end
end

%% Create the RBE2 connecting vector
% Assume that they do not follow the twist (this may need to be ammended)
ribVec = zeros(numel(WingBeamNodes),3);

for i = 1:numel(WingBeamNodes)
    if i == 1
        ribVec(i,1:3) = [ -1 0 0];
    elseif i == numel(WingBeamNodes)
        ribVec(i,1:3) = [ -1 0 0];
    elseif rbeOrient == 1
        ribVec(i,1:3) = cross([0 0 1],beamVec(i,1:3) / norm(beamVec(i,1:3)));
    else
        ribVec(i,1:3) = [ -1 0 0];
    end
end

%% Find the intersection of the Rib and leading / trailing edge vectors

lePoints = zeros(numel(WingBeamNodes),3);
tePoints = zeros(numel(WingBeamNodes),3);
count = 0;
for i = 1:numel(WingBeamNodes)
    solTl =  -1;
    solTt =  -1;
    j = 0;
    count = count +1;
    while( solTl > 1 || solTl < 0) % work out the solution as long as no solution exists
        j = j + 1; % Cycles through the different sections
        if ~strcmpi(part,'vtp')
            solTl = ((WingBeamNodes(i).coord(2) - leKinks(j,2)) - ((WingBeamNodes(i).coord(1) - leKinks(j,1)) * ribVec(count,2)) / ribVec(count,1)) / ...
                (leVec(j,2) - (leVec(j,1) * ribVec(count,2)) / ribVec(count,1));
        else
            solTl = ((WingBeamNodes(i).coord(3) - leKinks(j,3)) - ((WingBeamNodes(i).coord(1) - leKinks(j,1)) * ribVec(count,3)) / ribVec(count,1)) / ...
                (leVec(j,3) - (leVec(j,1) * ribVec(count,3)) / ribVec(count,1));
        end
        solTr = ( - (WingBeamNodes(i).coord(1) - leKinks(j,1)) + solTl * leVec(j,1)) / ribVec(count,1);
    end
    
    if solTl > 1 && solTl < 0
        warning('No solution found for the intersection of the rib and the leading edge')
    end
    
    %LEpoints(i,1:3) = NODE.coord(i,:) + 0.625 * solTR * RibVec(i,:);
    lePoints(count,1:3) = WingBeamNodes(i).coord + 1 * solTr * ribVec(count,:);
    j = 0;
    
    while (solTt > 1 || solTt < 0) % work out the solution as long as no solution exists
        j = j + 1; % Cycles through the different sections
        if ~strcmpi(part,'vtp')
            solTt = ((WingBeamNodes(i).coord(2) - teKinks(j,2)) - ((WingBeamNodes(i).coord(1) - teKinks(j,1)) * ribVec(count,2)) / ribVec(count,1)) / ...
                (teVec(j,2) - (teVec(j,1) * ribVec(count,2)) / ribVec(count,1));
        else
            solTt = ((WingBeamNodes(i).coord(3) - teKinks(j,3)) - ((WingBeamNodes(i).coord(1) - teKinks(j,1)) * ribVec(count,3)) / ribVec(count,1)) / ...
                (teVec(j,3) - (teVec(j,1) * ribVec(count,3)) / ribVec(count,1));
        end
        solTr = ( - (WingBeamNodes(i).coord(1) - teKinks(j,1)) + solTt * teVec(j,1)) / ribVec(count,1);
    end
    
    if solTt > 1 && solTt < 0
        warning('No solution found for the intersection of the rib and the trailing edge')
    end
    
    %TEpoints(i,1:3) = NODE.coord(i,:) + 0.4167 * solTR * RibVec(i,:);
    tePoints(count,1:3) = WingBeamNodes(i).coord + 1 * solTr * ribVec(count,:);
end

%% ROTATE POINTS ACCORDING TO TWIST
% twist = interp1(span * Fame.Results.Aerodynamic(1).eta,Fame.Results.Aerodynamic(1).jig_tw,struct2mat(Node,'coord',2));
% leTwistpoints(:,1) = lePoints(:,1) + (struct2mat(Node,'coord',1) - lePoints(:,1)) .* (1 - cos(pi * twist(:,1) / 180));
% leTwistpoints(:,2) = lePoints(:,2);
% leTwistpoints(:,3) = lePoints(:,3) + (struct2mat(Node,'coord',1) - lePoints(:,1)) .* sin(pi * twist(:,1) / 180);
% teTwistpoints(:,1) = tePoints(:,1) + (struct2mat(Node,'coord',1) - tePoints(:,1)) .* (1 - cos(pi * twist(:,1) / 180));
% teTwistpoints(:,2) = tePoints(:,2);
% teTwistpoints(:,3) = tePoints(:,3) + (struct2mat(Node,'coord',1) - tePoints(:,1)) .* sin(pi * twist(:,1) / 180);
% leTwistpoints = lePoints;
% teTwistpoints = tePoints;

%% Create new grid points for the RBE2s
% idx = find(find([Node.id]<(NodeIdx+100000) && [Node.Id]>NodeIdx));
% structNode = Node(idx);

for i = 1:size(lePoints,1)
    j = numel(Node);
    Node(j + 1).id       = NodeIdx + obj.Opts.Struct.LEOffsetID + i;
    Node(j + 1).cs       = 0;
    Node(j + 1).coord(1) = lePoints(i,1);
    Node(j + 1).coord(2) = lePoints(i,2);
    Node(j + 1).coord(3) = lePoints(i,3);
    Node(j + 1).cd       = 0;
    Node(j + 1).ps       = 0;
    Node(j + 1).seid     = 0;
    Node(j + 1).part     = PartName;
    Node(j + 1).type     = 'Aerodynamic';
end

for i = 1:size(tePoints,1)
    j = numel(Node);
    Node(j + 1).id       = NodeIdx + obj.Opts.Struct.TEOffsetID + i;
    Node(j + 1).cs       = 0;
    Node(j + 1).coord(1) = tePoints(i,1);
    Node(j + 1).coord(2) = tePoints(i,2);
    Node(j + 1).coord(3) = tePoints(i,3);
    Node(j + 1).cd       = 0;
    Node(j + 1).ps       = 0;
    Node(j + 1).seid     = 0;
    Node(j + 1).part     = PartName;
    Node(j + 1).type     = 'Aerodynamic';
end

%% Create RBE0 data array
for i = 1:size(lePoints,1)
    Rbe(i)         = Rbe0;
    Rbe(i).id      = NodeIdx + i;
    Rbe(i).grid    = WingBeamNodes(i).id;
    Rbe(i).data    = [Rbe(i).id + obj.Opts.Struct.LEOffsetID;Rbe(i).id + obj.Opts.Struct.TEOffsetID];
    Rbe(i).outline = Outline;
    Rbe(i).le      = leKinks;
    Rbe(i).te      = teKinks;
    Rbe(i).part    = PartName;
    
end

if isempty(obj.Mdl.Rbe0(1).id)
    obj.Mdl.Rbe0 = setFields(obj.Mdl.Rbe0,Rbe);
else
    obj.Mdl.Rbe0 = cat(2,obj.Mdl.Rbe0,Rbe);
end

obj.Mdl.Grid = setFields(obj.Mdl.Grid,Node);

end

