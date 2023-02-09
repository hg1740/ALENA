function [OPTIM] = RecoverBoxProp(OPTIM)

% Recover section parameters if necessary

OPTIM.t_skin = OPTIM.SOL(1:OPTIM.Nsec);
OPTIM.t_spar = OPTIM.SOL(OPTIM.Nsec+1:2*OPTIM.Nsec);
OPTIM.A_str  = OPTIM.SOL(2*OPTIM.Nsec+1:end);

OPTIM = BeamReduction(OPTIM);

end