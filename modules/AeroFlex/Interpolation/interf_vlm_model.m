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

function [Ic, In, Iv, Imv] = interf_vlm_model(fid, ncaero, NODE, AERO, MeshType)

% interface collocation points
%fprintf(fid, '\n\t - Assemblying collocation points interpolation matrix...'); 
Ic = assembly_colloc_interp_mat(ncaero, NODE, AERO);    
%fprintf(fid, 'done.'); 
% interface nodes
%fprintf(fid, '\n\t - Assemblying nodes interpolation matrix...'); 
In = assembly_node_interp_mat(ncaero, NODE, AERO);
%fprintf(fid, 'done.'); 
% interface vorticies
%fprintf(fid, '\n\t - Assemblying vorticies interpolation matrix...'); 

if MeshType > 3
    Iv = assembly_vortex_interp_mat_uvlm(ncaero, NODE, AERO);
else
    Iv = assembly_vortex_interp_mat(ncaero, NODE, AERO);
end

%fprintf(fid, 'done.'); 
% interface mid vorticies for force transfer to structure
%fprintf(fid, '\n\t - Assemblying vortices midpoint interpolation matrix...'); 
Imv = assembly_midv_interp_mat(ncaero, NODE, AERO);
%fprintf(fid, 'done.\n'); 

end
