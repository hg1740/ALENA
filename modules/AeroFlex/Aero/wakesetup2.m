function [lattice]=wakesetup2(lattice,state,ref)

infdist = 6 * ref.b_ref;

[a b c] = size(lattice.VORTEX);

switch (b)
    
    case 6
        
        c=[1 6];
        
    case 8
        
        c=[2 7];
        
    otherwise
        
        error('Wrong vortex struct size.');
        
end

V2 = lattice.VORTEX(:, c(1):c(2), :);

infx = infdist * cos(state.alpha) * cos(state.betha);
infy = infdist * sin(state.betha);
infz = infdist * sin(state.alpha) * cos(state.betha);

for t = 1:a
    for s = 1:2
        
        x = infx + lattice.VORTEX(t,c(s),1);
        y = infy + lattice.VORTEX(t,c(s),2);
        z = infz + lattice.VORTEX(t,c(s),3);
        
        psi = state.P/state.AS*x;
        theta = state.Q/state.AS*x;
        fi = state.R/state.AS*x;
        
        dx(t,s) = -x*(2-cos(theta) - cos(fi));
        dy(t,s) = sin(psi)*z-sin(fi) * x + (1-cos(psi)) * y;
        dz(t,s) = sin(theta) * x - sin(psi) * y + (1-cos(psi)) * z;
        
        
    end
end

for i=1:a
    
    INF1(i,1,1) = lattice.VORTEX(i,c(1),1) + infx + dx(i,1);
    INF1(i,1,2) = lattice.VORTEX(i,c(1),2) + infy + dy(i,1);
    INF1(i,1,3) = lattice.VORTEX(i,c(1),3) + infz + dz(i,1);
    
    INF2(i,1,1) = lattice.VORTEX(i,c(2),1) + infx + dx(i,2);
    INF2(i,1,2) = lattice.VORTEX(i,c(2),2) + infy + dy(i,2);
    INF2(i,1,3) = lattice.VORTEX(i,c(2),3) + infz + dz(i,2);
    
end

lattice.VORTEX = [INF1(:,1,:) V2 INF2(:,1,:)];

end