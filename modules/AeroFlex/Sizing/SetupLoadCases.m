%% Load Case definiton file

function data = SetupLoadCases(data,idx,LoadStruc)

M     = LoadStruc.M;
ALT   = LoadStruc.Alt;
URDD3 = LoadStruc.URDD3;

if isfield(LoadStruc,'CS')
    nc = length(LoadStruc.CS.Label);
else
    nc = 0;
end

% Define some conditions
data.Aero.state.SIMXZ = 1;  %Symmetry flag in XZ plane
data.Aero.state.SIMXY = 0;  %Symmetry flag in XY plane
data.Aero.state.AS = 1;     %Air speed
data.Aero.state.alpha = 0;  %AoA
data.Aero.state.betha = 0;  %Side-slip angle
data.Aero.state.P = 0;      %Roll rate
data.Aero.state.Q = 0;      %Pitch rate
data.Aero.state.R = 0;      %Yaw rate
data.Aero.state.pgcorr = 0;
data.Aero.state.M = M;      %Mach number
data.Aero.state.ALT  = ALT; %Altitude
[data.Aero.state.rho, data.Aero.state.p, ...
    data.Aero.state.T, data.Aero.state.a, ...
    data.Aero.state.mu] = ISA_h(data.Aero.state.ALT);   %Flight point data

% List some loadcases?
trimID = idx;
data.Aero.Trim.Select(idx) = idx;
data.Aero.Trim.ID(trimID) = trimID;
data.Aero.Trim.Mach(trimID) = M;
data.Aero.Trim.ALT(trimID)  = ALT;
data.Aero.Trim.Param(trimID).data(1,1)  = {'THRUST'};
data.Aero.Trim.Value(trimID).data(1,1)  = 0;
data.Aero.Trim.Param(trimID).data(2,1)  = {'SIDES'};
data.Aero.Trim.Value(trimID).data(2,1)  = 0;
data.Aero.Trim.Param(trimID).data(3,1)  = {'ROLL'};
data.Aero.Trim.Value(trimID).data(3,1)  = 0;
data.Aero.Trim.Param(trimID).data(4,1)  = {'PITCH'};
data.Aero.Trim.Value(trimID).data(4,1)  = 0;
data.Aero.Trim.Param(trimID).data(5,1)  = {'YAW'};
data.Aero.Trim.Value(trimID).data(5,1)  = 0;
data.Aero.Trim.Param(trimID).data(6,1)  = {'URDD1'};
data.Aero.Trim.Value(trimID).data(6,1)  = 0;
data.Aero.Trim.Param(trimID).data(7,1)  = {'URDD2'};
data.Aero.Trim.Value(trimID).data(7,1)  = 0;
data.Aero.Trim.Param(trimID).data(8,1)  = {'URDD3'};
data.Aero.Trim.Value(trimID).data(8,1)  = URDD3;
data.Aero.Trim.Param(trimID).data(9,1)  = {'URDD4'};
data.Aero.Trim.Value(trimID).data(9,1)  = 0;
data.Aero.Trim.Param(trimID).data(10,1) = {'URDD5'};
data.Aero.Trim.Value(trimID).data(10,1) = 0;
data.Aero.Trim.Param(trimID).data(11,1) = {'URDD6'};
data.Aero.Trim.Value(trimID).data(11,1) = 0;
data.Aero.Trim.Param(trimID).data(12,1) = {'HEAD'};
data.Aero.Trim.Value(trimID).data(12,1) = 0;
data.Aero.Trim.Param(trimID).data(13,1) = {'CLIMB'};
data.Aero.Trim.Value(trimID).data(13,1) = 0;
data.Aero.Trim.Param(trimID).data(14,1) = {'BANK'};
data.Aero.Trim.Value(trimID).data(14,1) = 0;

for i = 1:nc
    data.Aero.Trim.Param(trimID).data(14+i,1) = LoadStruc.CS.Label(i);
    data.Aero.Trim.Value(trimID).data(14+i,1) = LoadStruc.CS.Value(i);
end
data.Aero.Trim.NC(trimID) = 14 + nc;

end