function [cg,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(mass,length,radius,eff_factor)

switch nargin
    case 3
        eff_factor = 1.0;
end

Eff_length = eff_factor*length;

Ixx = 0.5*mass*radius^2;
Iyy = 0.25*mass*radius^2 + (1/12)*mass*Eff_length^2;
Izz = 0.25*mass*radius^2 + (1/12)*mass*Eff_length^2;
cg  = [length/2,0,0];
end