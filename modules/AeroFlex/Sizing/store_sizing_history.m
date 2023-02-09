function Sizing = store_sizing_history(OPTIM)


Sizing.I11(:,1)      = OPTIM.I1;
Sizing.I22(:,1)      = OPTIM.I2;
Sizing.I12(:,1)      = 0*OPTIM.I2;
Sizing.J(:,1)        = OPTIM.J;
Sizing.t_skin(:,1)   = OPTIM.t_skin;
Sizing.t_spar(:,1)   = OPTIM.t_spar;
if OPTIM.Type == 3
    Sizing.A_str(:,1)    = OPTIM.A_str;
elseif OPTIM.Type == 4
    Sizing.A_str(:,1)    = OPTIM.t_str;
end
Sizing.z_disp(:,1)   = zeros(OPTIM.Nsec + 1,1);
Sizing.twist(:,1)    = zeros(OPTIM.Nsec + 1,1);
Sizing.cl(:,1)       = zeros(OPTIM.Nsec,1);

end