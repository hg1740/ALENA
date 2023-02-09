function dydt = getQuatPos(t,y,v,t_vec)
for ii = 1:6
    v_in(ii,:) = interp1(t_vec,v(ii,:),t);
end
% pa    = v(1:3);
va    = v_in(1:3);
omega = v_in(4:6);
quat  = y(1:4);

OmegaQuat_a  = [ 0          omega(1)  omega(2)  omega(3);
                -omega(1)   0        -omega(3)  omega(2);
                -omega(2)   omega(3)  0        -omega(1);
                -omega(3)  -omega(2)  omega(1)  0       ];

CGa_e   = Quat2Rot(quat);

qa_dot = - 0.5*OmegaQuat_a*quat;
pa_dot =       CGa_e      *va;

dydt = [qa_dot;pa_dot];