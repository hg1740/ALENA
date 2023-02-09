function [rMatrix, cid] = getCoordSysRMatrix(InputData)
%getCoordSysRMatrix Extracts the rotation matrices from the .h5
%model data.

%Grab rotation matrices
RData     = InputData.COORDINATE_SYSTEM;
rotData   = RData.TRANSFORMATION.RDATA.DATA;
nData     = numel(rotData);
nCoordSys = nData / 12;
validateattributes(nCoordSys, {'numeric'}, {'scalar', 'integer'});
rotData   = reshape(rotData, [12, nCoordSys]);
cid       = RData.TRANSFORMATION.IDENTITY.CID;

%Index into rotation data to extract rotation matrix and origin of coordsys
r0 = rotData(1 : 3, :);
rMatrix = rotData(4 : end, :);
rMatrix = arrayfun(@(i) reshape(rMatrix(:, i), [3, 3]), 1 : nCoordSys, 'Unif', false);

end

