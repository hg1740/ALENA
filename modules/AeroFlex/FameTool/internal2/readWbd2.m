%% READWBD2   Read unsorted FAME structural mass data file
%
%  READWBD2(CLASSNAME) reads mass data from FAME WBD2 file. Thousands of
%  small mass elements are created. The coordinates and mass values are
%  defined in the WBD2 file.
%
%  CLASSNAME is a Fame2mat object
%
%  This is a private function and should be accessed from the |Fame2mat|
%  object.
%
%   Example:
%      obj = fame2mat;
%      obj = getFameResults(obj); % calls this function
%
%

%   Copyright 2016 University of Bristol
%   Private function.
function d = readWbd2( varargin )

if ~isempty(varargin)
    fileName = varargin{1};
    pathName = [];
else
    [fileName,pathName] = uigetfile(' * .205','FAME Mass File');
end

if fileName == 0
    return
end

blk = textread(fullfile(pathName,fileName),'%s','delimiter','\n');


idx            = find(~cellfun(@isempty,regexp(blk,'^C')));
blk(idx)       = [];
blk            = char(blk);
idx            = findstr(blk(1,:),'!');
blk(:,idx:end) = [];
blk            = cellstr(blk);

blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');

blk = char(blk);

blk = str2num(blk);

d.x    = blk(:,4);
d.y    = blk(:,5);
d.z    = blk(:,6);
d.Mass = blk(:,7);

end

