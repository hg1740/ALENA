function R = generateRotationMatrix(eX, eY_temp)
%generateRotationMatrix Forms a [3, 3, N] rotation matrix R = [eX, eY, eZ]
%using the vector 'eX' and a vector eY_temp' that lies in the desired XY
%plane.
%
% Accepts a matrix of eX and matrix of eY points, i.e. eX = eY = [3, N].

eZ = cross(eX, eY_temp);
eY = cross(eX, eZ);

eX = eX ./ vecnorm(eX);
eY = eY ./ vecnorm(eY);
eZ = eZ ./ vecnorm(eZ);

R  = permute(cat(3, eX, eY, eZ), [1, 3, 2]);

end