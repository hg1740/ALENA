
function [iel, jel, vel] = nodal_mass_matrix(NMAT, dof, NODER)

% transform node mass matrix into global reference frame
% involved dofs
index = find(dof); % look for constraints
R     = [NODER, zeros(3,3); zeros(3,3), NODER];
M     = R * NMAT * R';

Mel   = M(index, index);
nr    = length(index);

% sparse mat row index		
iel = double(repmat(dof(index)', nr, 1));

% TODO:  CHECK THIS CHANGE!
jel = repmat(double(dof(index)),nr,1);
jel = jel(:);
vel = Mel(:);

end
