function  Param = ParasiteDragFuse(body,state,ref,Param)

%% Calculate Re based on local chord
Re_chord = state.rho*state.AS^2*body.Length/state.mu;

%% Roughness coefficient
k = 0.634*10^-5; % Smooth paint
%k = 0.405*10^-5; % Production sheet metal
%k = 0.152*10^-5; % Polished sheet metal
%k = 0.052*10^-5; % Smooth molded composite

%% Interference factor
Q = 1.0;

%% Calculate cutoff Reynolds number
Re_cutoff = 38.21*(body.Length/k).^1.053; % Subsonic
%Re_cutoff = 44.62*(geo.c/k).^1.053*state.M^(1.16); % Transonic/Supersonic

Re_min = min(Re_cutoff,Re_chord);

radius = max(body.radius);

%% Component Form Factor
f = body.Length/(2*radius);
FF = (1 + 60/(f^3) + f/400);

Swet_c = body.Length*pi*(radius)^2;
Swet_c = 0.9*Swet_c; % Accounts for nose and tail cone

%% Laminar/Turbulent Drag
Cf_turb = 0.455./(log10(Re_min).^2.58*(1+0.144*state.M^2)^0.65);

%% Combined Drag

if state.SIMXZ == 1
    Swet_c = 0.5*Swet_c;
end

Cd0= (Cf_turb.*FF*Q.*Swet_c)/ref.S_ref;
d0 = Cd0*(0.5*state.rho*state.AS^2*ref.S_ref);

Param.Swet_c = Swet_c;

Param.Cf_turb = Cf_turb;
Param.Cd0 = Cd0;
Param.d0 = d0;
Param.Re_chord = Re_chord;

end