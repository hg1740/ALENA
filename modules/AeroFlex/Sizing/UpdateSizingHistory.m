function Hist = UpdateSizingHistory(Hist,optim)


Hist.I11        = [Hist.I11,optim.I1];
Hist.I22        = [Hist.I22,optim.I2];
Hist.I12        = [Hist.I12,0*optim.I2];
Hist.J          = [Hist.J,optim.J];
Hist.t_skin     = [Hist.t_skin,optim.t_skin];
Hist.t_spar     = [Hist.t_spar,optim.t_spar];
if optim.Type == 3
    Hist.A_str  = [Hist.A_str,optim.A_str];
elseif optim.Type == 4
    Hist.A_str  = [Hist.A_str,optim.t_str];
end

end