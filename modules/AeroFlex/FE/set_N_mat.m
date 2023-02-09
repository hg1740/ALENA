%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 - 2011 
% 
% Sergio Ricci (sergio.ricci@polimi.it)
%
% Politecnico di Milano, Dipartimento di Ingegneria Aerospaziale
% Via La Masa 34, 20156 Milano - ITALY
% 
% This file is part of NeoCASS Software (www.neocass.org)
%
% NeoCASS is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation;
% either version 2, or (at your option) any later version.
%
% NeoCASS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public
% License along with NeoCASS; see the file GNU GENERAL 
% PUBLIC LICENSE.TXT.  If not, write to the Free Software 
% Foundation, 59 Temple Place -Suite 330, Boston, MA
% 02111-1307, USA.
%

%
%***********************************************************************************************************************
%  SimSAC Project
%
%  SMARTCAD
%  Simplified Models for Aeroelasticity in Conceptual Aircraft Design  
%
%                      Sergio Ricci         <ricci@aero.polimi.it>
%                      Luca Cavagna         <cavagna@aero.polimi.it>
%                      Alessandro Degaspari <degaspari@aero.polimi.it>
%                      Luca Riccobene       <riccobene@aero.polimi.it>
%
%  Department of Aerospace Engineering - Politecnico di Milano (DIAPM)
%  Warning: This code is released only to be used by SimSAC partners.
%  Any usage without an explicit authorization may be persecuted.
%
%***********************************************************************************************************************
%	
%   Author: Luca Cavagna, Pierangelo Masarati, DIAPM
%***********************************************************************************************************************

function N = set_N_mat(c01, c02, c03, c1, c2, c3, f1, f2, f3)

N = zeros(12,18);

x = [c01(1) c02(1) c03(1)];
y = [c01(2) c02(2) c03(2)];
z = [c01(3) c02(3) c03(3)];

% COLLOC 1
eta = -1/sqrt(3);
NI = [(0.5 * eta * (eta-1)) (1-eta^2) (0.5 * eta * (eta+1))];
NIp = [(0.5 * (2 * eta-1)) (-2 * eta) (0.5 * (2 * eta+1))]; % shape functions derivative evaluated in COLLOC1

%JI=NIp*[x;y;z]';
JI = [NIp(1)*x(1)+NIp(2)*x(2)+NIp(3)*x(3),NIp(1)*y(1)+NIp(2)*y(2)+NIp(3)*y(3),NIp(1)*z(1)+NIp(2)*z(2)+NIp(3)*z(3)];
NIp = NIp ./ norm(JI);

% COLLOC 2
eta = +1/sqrt(3);
NII = [(0.5 * eta * (eta-1)) (1-eta^2) (0.5 * eta * (eta+1))];
NIIp = [(0.5 * (2 * eta-1)) (-2 * eta) (0.5 * (2 * eta+1))]; % shape functions derivative evaluated in COLLOC2 

%JII=NIIp*[x;y;z]';
JII = [NIIp(1)*x(1)+NIIp(2)*x(2)+NIIp(3)*x(3),NIIp(1)*y(1)+NIIp(2)*y(2)+NIIp(3)*y(3),NIIp(1)*z(1)+NIIp(2)*z(2)+NIIp(3)*z(3)];
NIIp = NIIp ./ norm(JII);

x = [c1(1) c2(1) c3(1)];
y = [c1(2) c2(2) c3(2)];
z = [c1(3) c2(3) c3(3)];

%PpI = [x;y;z]*NIp';
PpI = [NIp(1)*x(1)+NIp(2)*x(2)+NIp(3)*x(3);NIp(1)*y(1)+NIp(2)*y(2)+NIp(3)*y(3);NIp(1)*z(1)+NIp(2)*z(2)+NIp(3)*z(3)];

%PpII = [x;y;z]*NIIp';
PpII = [NIIp(1)*x(1)+NIIp(2)*x(2)+NIIp(3)*x(3);NIIp(1)*y(1)+NIIp(2)*y(2)+NIIp(3)*y(3);NIIp(1)*z(1)+NIIp(2)*z(2)+NIIp(3)*z(3)];

I6 = eye(6);
N(1:6,1:6)    = NIp(1)*I6;
N(1:6,7:12)   = NIp(2)*I6;
N(1:6,13:18)  = NIp(3)*I6;
N(7:12,1:6)   = NIIp(1)*I6;
N(7:12,7:12)  = NIIp(2)*I6;
N(7:12,13:18) = NIIp(3)*I6;

% FIRST section
% node 1
%N(1:3,4:6) = -crossm(NIp(1) .* f1 - PpI .* NI(1)); 
N(1:3,4:6)   = -crossm([NIp(1)*f1(1) - PpI(1)*NI(1);NIp(1)*f1(2)- PpI(2)*NI(1);NIp(1)*f1(3)- PpI(3)*NI(1)]); 
N(1:3,10:12) = -crossm([NIp(2)*f2(1) - PpI(1)*NI(2);NIp(2)*f2(2)- PpI(2)*NI(2);NIp(2)*f2(3)- PpI(3)*NI(2)]); 
N(1:3,16:18) = -crossm([NIp(3)*f3(1) - PpI(1)*NI(3);NIp(3)*f3(2)- PpI(2)*NI(3);NIp(3)*f3(3)- PpI(3)*NI(3)]);  

% SECOND section
% node 1
%N(7:9,4:6) = -crossm(NIIp(1) .* f1 - PpII .* NII(1)); 
N(7:9,4:6)   = -crossm([NIIp(1)*f1(1) - PpII(1)*NII(1);NIIp(1)*f1(2)- PpII(2)*NII(1);NIIp(1)*f1(3)- PpII(3)*NII(1)]); 
N(7:9,10:12) = -crossm([NIIp(2)*f2(1) - PpII(1)*NII(2);NIIp(2)*f2(2)- PpII(2)*NII(2);NIIp(2)*f2(3)- PpII(3)*NII(2)]); 
N(7:9,16:18) = -crossm([NIIp(3)*f3(1) - PpII(1)*NII(3);NIIp(3)*f3(2)- PpII(2)*NII(3);NIIp(3)*f3(3)- PpII(3)*NII(3)]);   
end
