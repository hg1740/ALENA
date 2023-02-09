function MASS = WBnCG(NODE,PARAM,CONM,MASS,BAR,BEAM,INFO)


% set lumped mass CG
[MASS.CG, MASS.MCG, MASS.MRP] = ...
    wb_set_conm_mass(INFO.nconm, NODE.Index, NODE.Coord, NODE.R, PARAM.GRDPNT, CONM);
% set bar mass CG
[MASS.CG, MASS.MCG, MASS.MRP] =...
    wb_add_bar_mass(INFO.nbar, NODE.Coord, NODE.R, MASS.CG, PARAM.GRDPNT, MASS.MCG, MASS.MRP, BAR);
% set beam mass CG
[MASS.CG, MASS.MCG, MASS.MRP] =...
    wb_add_bar_mass(INFO.nbeam, NODE.Coord, NODE.R, MASS.CG, PARAM.GRDPNT, MASS.MCG, MASS.MRP, BEAM);
% get principal axes
%[MASS.MCG_pa, MASS.R_pa] = wb_principal_axis(MASS.MCG);

end