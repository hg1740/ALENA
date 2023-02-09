function ftilde = f2ftilde(f)
%f2ftilde Returns a block matrix of skew-symmetric terms related to the 3
%local beam loads and 3 local beam moments. 
%
% If 'f' is a 6 x 1 vector of loads i.e. {f} = [F ; M] then 'f_tilde'
% is [0, F_tilde ; F_tilde, M_tilde], where 'F_tilde' and 'M_tilde' are
% 3 x 3 skew-symmetric matrices. 
%
% See Equation 12 of Ref. [1] for more details.
%
% References: 
%   [1]. "Aeroelastic Modelling of Highly Flexible Wings", Howcroft et. al
%   SciTech 2016.

%Beam forces in local (intrinsic) frame
Ftilde = skew(f(1 : 3));

%Beam moments in local (intrinsic) frame
Mtilde = skew(f(4 : 6));
    
%Matrix of skew-symmetric terms
ftilde = [zeros(3), Ftilde ; ...
          Ftilde  , Mtilde];
      
end
