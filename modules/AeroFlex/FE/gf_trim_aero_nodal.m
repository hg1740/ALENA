function Fa = gf_trim_aero_nodal(INFO, DOF, NODE, AERO, results, AeroLinear)

% generalized aero forces vector on structural master nodes
Fa = zeros(INFO.ndof, 1);
% check
ngrid = INFO.ngrid;

if AeroLinear == 1
    Fxa = results.F(:,1);
    Fya = results.F(:,2);
    Fza = results.F(:,3);
else
    Fxa = results.F(:,1);
    Fya = results.F(:,2);
    Fza = results.F(:,3);
end

% get structural nodal forces
Fxs = zeros(ngrid, 1);	     Fys = zeros(ngrid, 1);       Fzs = zeros(ngrid, 1);
Imat = (AERO.Interp.Imv)';
Fxs = Imat * Fxa;            Fys = Imat * Fya;            Fzs = Imat * Fza;

% master node forces
Fxm = zeros(ngrid, 1); 	     Fym = zeros(ngrid, 1); 	    Fzm = zeros(ngrid, 1);
Mxm = zeros(ngrid, 1); 	     Mym = zeros(ngrid, 1); 	    Mzm = zeros(ngrid, 1);

for n = 1:ngrid
    if (NODE.Index(n)) % if master node
        Fxm(n) = Fxs(n);
        Fym(n) = Fys(n);
        Fzm(n) = Fzs(n);
    end
end

% get slave forces
for n = 1:ngrid
    
    if (~isempty(NODE.Aero.Coord(n).data)) % master
        
        % forces on slaves
        Ndata = NODE.Aero.Index(n).data;
        fxs = Fxs(Ndata);
        fys = Fys(Ndata);
        fzs = Fzs(Ndata);
        
        Fxm(n) = sum(fxs) + Fxm(n);
        Fym(n) = sum(fys) + Fym(n);
        Fzm(n) = sum(fzs) + Fzm(n);
        
        % This needs to be expanded to make faster!
        %M = sum( cross(NODE.Aero.Coord(n).data', [fxs, fys, fzs], 2), 1);
        Ndata = NODE.Aero.Coord(n).data';
        M = sum([Ndata(:,2).*fzs - Ndata(:,3).*fys,-Ndata(:,1).*fzs + Ndata(:,3).*fxs,Ndata(:,1).*fys - Ndata(:,2).*fxs],1);
        
        Mxm(n) = Mxm(n) + M(1);
        Mym(n) = Mym(n) + M(2);
        Mzm(n) = Mzm(n) + M(3);
    end
    
end % end node loop

F = zeros(ngrid, 3);
M = zeros(ngrid, 3);

F = [Fxm, Fym, Fzm];
M = [Mxm, Mym, Mzm];

% store master node forces and moments in the correct DOF position

for n = 1:ngrid
    
    if (NODE.Index(n)) % if master node
        
        index = find(DOF(n, 1:3)); % get free dofs
        pos = DOF(n, index);
        Fa(pos, 1) = F(n, index);  % assembly
        
        index = find(DOF(n, 4:6)); % get free dofs
        if ~isempty(index)
            
            indexoff = index + 3;
            pos = DOF(n, indexoff);
            Fa(pos, 1) = M(n, index);  % assembly
            
        end
        
    end
end

end
