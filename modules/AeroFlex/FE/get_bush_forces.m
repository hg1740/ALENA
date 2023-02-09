function F = get_bush_forces(ngrid,CBUSH,SOL,Node,NODER,ROT,BUSHROT)

F = zeros(ngrid, 6);

for k = 1:length(CBUSH.ID)
    
    n1 = CBUSH.Node(k,1);
    n2 = CBUSH.Node(k,2);
    
    c1 = Node.Coord(n1,:)' + SOL(n1,1:3)';
    c2 = Node.Coord(n2,:)' + SOL(n2,1:3)';
    
    r1 = SOL(n1,4:6)';
    r2 = SOL(n2,4:6)';
    
    %         R1 = NODER(:,:,n1);
    R1 = BUSHROT(:,:,1,k);
    R2 = NODER(:,:,n2);
    
    ROT1 = BUSHROT(:,:,1,k)'*ROT(n1,1:3)';
    ROT2 = BUSHROT(:,:,2,k)'*ROT(n2,1:3)';
%     ROT1 = ROT(n1,1:3)';
%     ROT2 = ROT(n2,1:3)';
    
    rdiff = ROT2 - ROT1;
    %rdiff = fliplr(rotm2eul(R2*R1'))';
    
    G1 = Gmat(ROT1);
    
    u = R1'*(c2-c1);
    psi = G1*(rdiff);
    
    K1 = diag(CBUSH.K(k).data(1:3)');
    K2 = diag(CBUSH.K(k).data(4:6)');
    
    S1 = R1*K1;
    S2 = G1'*K2;
    
    Fel(1:6) = [-S1*u;-S2*psi]; 
    Fel(7:12) = [S1*u;S2*psi];
    
    F(n1, 1:6) = F(n1, 1:6) + Fel(1:6);
    F(n2, 1:6) = F(n2, 1:6) + Fel(7:12);
    
end