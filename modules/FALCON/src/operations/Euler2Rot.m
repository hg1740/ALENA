function ROT = Euler2Rot(euler)

ROT1 = [1  0             0            ;
        0  cos(euler(1)) sin(euler(1));
        0 -sin(euler(1)) cos(euler(1))];
    
ROT2 = [ cos(euler(2)) 0 sin(euler(2));
         0             1 0            ;
        -sin(euler(2)) 0 cos(euler(2))];
    
ROT3 = [ cos(euler(3)) sin(euler(3)) 0;
        -sin(euler(3)) cos(euler(3)) 0;
         0             0             1];
    
ROT  = ROT1*ROT2*ROT3;