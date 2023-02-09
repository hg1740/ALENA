function [MAT,INFO] = MAT1Creator(MAT,PARAM,INFO,ID)

nmat = INFO.nmat;
nmat = nmat + 1;
MAT.ID(nmat) = ID;
MAT.E(nmat) = PARAM.E;
MAT.nu(nmat) = PARAM.nu;
MAT.G(nmat) = MAT.E(nmat) / (2 * (1 + MAT.nu(nmat)));
MAT.Rho(nmat) = PARAM.Rho;
MAT.ST(nmat) = 0;
MAT.SC(nmat) = 0;
MAT.SS(nmat) = 0;
INFO.nmat = nmat;
end