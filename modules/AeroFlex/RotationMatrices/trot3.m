function[p2]=trot3(hinge,p,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TROT: Auxillary rotation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotates point p around hinge alpha rads.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ref: 	R�de, Westergren, BETA 4th ed,
%			studentlitteratur, 1998
%			pp:107-108
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: 	Tomas Melin, KTH,Department of%
% 				aeronautics, Copyright 2000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Context:	Auxillary function for
%				TORNADO.
% Called by: setrudder, normals
% Calls:		norm (MATLAB std fcn)
%				sin			"
%				cos			"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELP:		Hinge=vector around rotation
%						takes place.
%				p=point to be rotated
%				alpha=radians of rotation
%				3D-workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=hinge(1);
b=hinge(2);
c=hinge(3);

rho=sqrt(a^2+b^2);
r=sqrt(a^2+b^2+c^2);

if r==0
    cost=0;
    sint=1;
else
    cost=c/r;
    sint=rho/r;
end

if rho==0
    cosf=0;
    sinf=1;
else
    cosf=a/rho;
    sinf=b/rho;
end

cosa=cos(alpha);
sina=sin(alpha);

RZF=[[cosf -sinf 0];[sinf cosf 0];[0 0 1]];
RYT=[[cost 0 sint];[0 1 0];[-sint 0 cost]];
RZA=[[cosa -sina 0];[sina cosa 0];[0 0 1]];
RYMT=[[cost 0 -sint];[0 1 0];[sint 0 cost]];
RZMF=[[cosf sinf 0];[-sinf cosf 0];[0 0 1]];

P=RZF*RYT*RZA*RYMT*RZMF;
p2=P*p';
end