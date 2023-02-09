function lattice_defo = rotate_control_norm(ref, state, geo, lattice, ROT_VALUE, HINGE, varargin)

PLOT_hinge = 1;
nc = length(ROT_VALUE);
lattice_defo = lattice;
h1 = zeros(1,3);
h2 = zeros(1,3);

for n=1:nc
    
    if ROT_VALUE(n)
        
        % Store the inboard and outboard HINGE coordinates
        h1(1, 1:3) = HINGE(n, 1, 1:3);
        h2(1, 1:3) = HINGE(n, 2, 1:3);
        
        patch = lattice.Control.Patch(n);
        partition = lattice.Control.Part(n);
        
        DOF = (lattice.DOF(patch, partition, 1) : lattice.DOF(patch, partition, 2)); % patch dof
        lattice_defo = set_control_rot(patch, partition, ROT_VALUE(n), lattice_defo, geo, h1, h2, ...
            lattice.Control.DOF(n).data);
        
        
    end
end

% update wake position
lattice_defo = defo_wakesetup(lattice_defo, state, ref);

end

function lattice = set_control_rot(wing, division, deflection, lattice, geo, a1, b1,control_dof)

ny = geo.ny(wing, division);
fnx2 = length(control_dof)/ny;
fnx = fnx2;
h = zeros(1,3);
hinge = zeros(1,3);

% Define hinge vector
h = b1 - a1;
h(1) = 0;

% Normalise hinge vector
hinge = h ./ norm(h); 

for n=1:ny
    offset = (n-1) * fnx;
    control_row(n, :) = control_dof((1 + offset : fnx + offset));
end

% Rotate control normal
RotMatrix = RotVecHinge(hinge(1),hinge(2),hinge(3),deflection);
for n = 1:ny
    for k=1:fnx
        lattice.N(control_row(n, k), 1:3) = (RotMatrix * lattice.N(control_row(n, k), 1:3)')';
        lattice.CSRM(:,:,control_row(n, k)) = RotMatrix;
    end   
end

% plot_vlm(lattice)
% hold on;plot3([a1(1),b1(1)],[a1(2),b1(2)],[a1(3),b1(3)],'b-')
% hold on;quiver3(mean([a1(1),b1(1)]),mean([a1(2),b1(2)]),mean([a1(3),b1(3)]),-lattice.N(control_row(n, k), 1),-lattice.N(control_row(n, k), 2),-lattice.N(control_row(n, k), 3),'r')
% axis equal
% rotate_N = (RotMatrix * lattice.N(control_row(n, k), 1:3)')';hold on;
% quiver3(mean([a1(1),b1(1)]),mean([a1(2),b1(2)]),mean([a1(3),b1(3)]),-rotate_N(1, 1),-rotate_N(1, 2),-rotate_N(1, 3),'b')


end