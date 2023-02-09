function lattice_defo = update_hinge_lines(geo, lattice_defo)

% Loop through all the control surfaces
for nc = 1:geo.nc
    
    pind = lattice_defo.Control.Patch(nc);
    
    % Recall the element number of the first element on the control surface
    n1 = lattice_defo.Control.DOF(nc).data(1);
    
    % Recall the number of elements that make up the control surface
    ncptot = length(lattice_defo.Control.DOF(nc).data);
    
    %  Recall the element number of the last LE element on the control surface
    n2 = lattice_defo.Control.DOF(nc).data(ncptot-geo.fnx(pind)+1);
    
    i1 = 1;
    i2 = 2;
    
    if (geo.b(pind)<0)
        i1 = 2;
        i2 = 1;
    end
    
    lattice_defo.Control.Hinge(nc,1,:) = lattice_defo.XYZ(n1,i1,:);
    lattice_defo.Control.Hinge(nc,2,:) = lattice_defo.XYZ(n2,i2,:);
    
end

end
