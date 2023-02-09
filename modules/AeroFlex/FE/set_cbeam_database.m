
function BEAM = set_cbeam_database(ncbeam, BEAM, PBEAM, MAT, NODE, COORD)

for n=1:ncbeam
    
    % assemble stiffnes matrix
    p = BEAM.PID(n);
    i = PBEAM.Mat(p);
    
    rho1node = MAT.Rho(i) + PBEAM.RhoNS(p).data(1); % structural density
    rho2node = MAT.Rho(i) + PBEAM.RhoNS(p).data(2); % structural density
    rho1col = MAT.Rho(i) + PBEAM.DRhoNS(1,1,p); % structural density
    rho2col = MAT.Rho(i) + PBEAM.DRhoNS(1,2,p); % structural density
    
    L = zeros(1,3); % body length
    L(1) = norm(NODE.Coord(BEAM.Conn(n,1),:) + BEAM.Offset(n, 1:3) - BEAM.Colloc(1, :, n));
    L(2) = norm(BEAM.Colloc(1, :, n) - BEAM.Colloc(2, :, n));
    L(3) = norm(NODE.Coord(BEAM.Conn(n,3),:) + BEAM.Offset(n, 7:9) - BEAM.Colloc(2, :, n));
    %
    
    EA1  =  PBEAM.DA(1,1,p) * MAT.E(i);
    EA2  =  PBEAM.DA(1,2,p) * MAT.E(i);
    
    Kshear = PBEAM.Kshear(p,:);
    
    if Kshear(1) == 0
        % Bernoulli Beam Theory
        GAy1 =  EA1;
        GAy2 =  EA2;
    else
        % Timoshenko Beam Theory
        GAy1 =  PBEAM.DA(1,1,p) * MAT.G(i) * Kshear(1);
        GAy2 =  PBEAM.DA(1,2,p) * MAT.G(i) * Kshear(1);
        
    end
    
    if Kshear(2) == 0
        % Bernoulli Beam Theory
        GAz1 =  EA1;
        GAz2 =  EA2;
    else
        % Timoshenko Beam Theory
        GAz1 =  PBEAM.DA(1,1,p) * MAT.G(i) * Kshear(2);
        GAz2 =  PBEAM.DA(1,2,p) * MAT.G(i) * Kshear(2);
        
    end
    
    GJ1  =  PBEAM.DJ(1,1,p) * MAT.G(i);
    GJ2  =  PBEAM.DJ(1,2,p) * MAT.G(i);
    
    EJy1 =  PBEAM.DI(2,1,p) * MAT.E(i);
    EJz1 =  PBEAM.DI(1,1,p) * MAT.E(i);
    EJyz1 = PBEAM.DI(3,1,p) * MAT.E(i);
    
    EJy2 =  PBEAM.DI(2,2,p) * MAT.E(i);
    EJz2 =  PBEAM.DI(1,2,p) * MAT.E(i);
    EJyz2 = PBEAM.DI(3,2,p) * MAT.E(i);
    
    D1 = diag([EA1, GAy1, GAz1, GJ1, EJy1, EJz1]);
    
    % set moment
    D1(5,6) = -EJyz1;
    D1(6,5) = D1(5,6);
    
    D2 = diag([EA2, GAy2, GAz2, GJ2, EJy2, EJz2]);
    % set moment
    D2(5,6) = -EJyz2;
    D2(6,5) = D2(5,6);
    
    % store section stiffness matrices in the two collocation points
    BEAM.D(:, :,1, n) = D1;
    BEAM.D(:, :,2, n) = D2;
    
    % set lumped nodes mass matrix
    M  = zeros(1,3);
    J  = zeros(3,3);
    MJ = zeros(1,3);
    v  = zeros(1,3);
    Jx = zeros(1,3);
    
    % Assumes a linear distribution between points
    M(1) = (PBEAM.A(1,p).data(1) * rho1node  +  PBEAM.RhoNS(1,p).data(1) + ...
        PBEAM.DA(1,1,p) * rho1col  +  PBEAM.DRhoNS(1,1,p)) * L(1)/2;
    M(2) = (PBEAM.DA(1,1,p) * rho1col  +  PBEAM.DRhoNS(1,1,p) + ...
        PBEAM.DA(1,2,p) * rho2col  +  PBEAM.DRhoNS(1,2,p)) * L(2)/2;
    M(3) = (PBEAM.A(1,p).data(2) * rho2node  +  PBEAM.RhoNS(1,p).data(2) + ...
        PBEAM.DA(1,2,p) * rho2col  +  PBEAM.DRhoNS(1,2,p)) * L(1)/2;
    
    MJ(1) = 1/12 * (M(1) * L(1)^2);  % rotational inertia respect to CG1
    MJ(2) = 1/12 * (M(2) * L(2)^2);  % rotational inertia respect to CG2
    MJ(3) = 1/12 * (M(3) * L(1)^2);  % rotational inertia respect to CG3
    
    % inertia along neutral axis
    Jx(1) = ((rho1node + PBEAM.RhoNS(1,p).data(1)/PBEAM.A(1,p).data(1))*(PBEAM.I(1,p).data(1,1) + PBEAM.I(1,p).data(2,1)) + ...
        (rho1col + PBEAM.DRhoNS(1,1,p)/PBEAM.DA(1,1,p))*(PBEAM.DI(1,1,p) + PBEAM.DI(2,1,p))) * L(1)/2;
    Jx(2) = ((rho1col + PBEAM.DRhoNS(1,1,p)/PBEAM.DA(1,1,p))*(PBEAM.DI(1,1,p) + PBEAM.DI(2,1,p)) + ...
        (rho2col + PBEAM.DRhoNS(1,2,p)/PBEAM.DA(1,2,p))*(PBEAM.DI(1,2,p) + PBEAM.DI(2,2,p))) * L(2)/2;
    Jx(3) = ((rho2node + PBEAM.RhoNS(1,p).data(2)/PBEAM.A(1,p).data(2))*(PBEAM.I(1,p).data(1,2) + PBEAM.I(1,p).data(2,2)) + ...
        (rho2col + PBEAM.DRhoNS(1,2,p)/PBEAM.DA(1,2,p))*(PBEAM.DI(1,2,p) + PBEAM.DI(2,2,p))) * L(1)/2;
    
    BEAM.M(1:3,1:3, 1, n) = diag([M(1), M(1), M(1)]);
    BEAM.M(1:3,1:3, 2, n) = diag([M(2), M(2), M(2)]);
    BEAM.M(1:3,1:3, 3, n) = diag([M(3), M(3), M(3)]);
    
    % body 1 MASS MATRIX in NODE reference frame
    S = zeros(3,3);
    v = L(1)/2 .* BEAM.R(:,1,1,n)' + BEAM.Offset(n, 1:3); % cg body 1
    S(2,1) =  M(1) * v(3); S(3,1) = -M(1) * v(2); S(1,2) = -S(2,1); S(3,2) =  M(1) * v(1); S(1,3) = -S(3,1); S(2,3) = -S(3,2);
    BEAM.M(4:6,1:3,1,n) = S; BEAM.M(1:3,4:6,1,n) = S';
    
    J = diag([Jx(1), MJ(1), MJ(1)]);
    % avoid round-off errors
    RJ = BEAM.R(:, :, 1, n) * sqrt(J);
    BEAM.M(4:6,4:6, 1, n) = RJ * RJ' - M(1) .* (crossm(v) * crossm(v));
    
    % body 2 MASS MATRIX in NODE reference frame
    S = zeros(3,3);
    v = BEAM.Offset(n, 1:3); % cg body 2
    S(2,1) =  M(2) * v(3); S(3,1) = -M(2) * v(2); S(1,2) = -S(2,1); S(3,2) =  M(2) * v(1); S(1,3) = -S(3,1); S(2,3) = -S(3,2);
    BEAM.M(4:6,1:3,2,n) = S; BEAM.M(1:3,4:6,2,n) = S';
    
    J = diag([Jx(2), MJ(2), MJ(2)]);
    RJ = BEAM.R(:, :, 2, n) * sqrt(J);
    BEAM.M(4:6,4:6,2, n) = RJ * RJ' - M(2) .* (crossm(v) * crossm(v));
    
    % body 3 MASS MATRIX in NODE reference frame
    S = zeros(3,3);
    v = -L(3)/2 .* BEAM.R(:,1,3,n)' + BEAM.Offset(n, 7:9); % cg body 3
    S(2,1) =  M(3) * v(3); S(3,1) = -M(3) * v(2); S(1,2) = -S(2,1); S(3,2) =  M(3) * v(1); S(1,3) = -S(3,1); S(2,3) = -S(3,2);
    BEAM.M(4:6,1:3,3,n) = S; BEAM.M(1:3,4:6,3,n) = S';
    
    J = diag([Jx(3), MJ(3), MJ(3)]);
    RJ = BEAM.R(:, :, 3, n) * sqrt(J);
    BEAM.M(4:6,4:6,3, n) = RJ * RJ' - M(3) .* (crossm(v) * crossm(v));
    
end % end CBEAM loop

end
