function vtilde = v2vtilde(v)
%v2vtilde Returns a block matrix of skew-symmetric terms related to the 3
%local beam linear and angular velocities.
%
% If 'v' is a 6 x 1 vector of velocities i.e. {v} = [V ; O] then 'v_tilde'
% is [O_tilde, 0 ; V_tilde, O_tilde], where 'O_tilde' and 'V_tilde' are
% 3 x 3 skew-symmetric matrices. 
%
% See Equation 11 of Ref. [1] for more details.
%
% References: 
%   [1]. "Aeroelastic Modelling of Highly Flexible Wings", Howcroft et. al
%   SciTech 2016.

%Linear velocites
Vtilde = skew(v(1:3));

%Angualr velocities
Otilde = skew(v(4:6));

%Matrix of skew-symmetric terms
vtilde = [Otilde, zeros(3); ...
          Vtilde, Otilde];
            
end