function Param = ParasiteDragWing(geo,state,ref,Param)

mean_t_c = Param.t_c;
mean_x_c = Param.x_c;
Swet_c   = Param.Swet_c;
Sref_c   = Param.Sref_c;

%% Calculate Re based on local chord
Re_chord = state.rho*state.AS*geo.c/state.mu;

%% Critical Reynolds number and transition point
Re_crit = 1*10^6;
x_cr = state.mu*Re_crit/(state.rho*state.AS);
tr = x_cr./geo.c;

tr(tr>1) = 1;

%% Roughness coefficient
k = 0.634*10^-5; % Smooth paint
%k = 0.405*10^-5; % Production sheet metal
%k = 0.152*10^-5; % Polished sheet metal
%k = 0.052*10^-5; % Smooth molded composite

%% Interference factor
Q = 1.0;

%% Calculate cutoff Reynolds number
Re_cutoff = 38.21*(geo.c/k).^1.053; % Subsonic
%Re_cutoff = 44.62*(geo.c/k).^1.053*state.M^(1.16); % Transonic/Supersonic

Re_min = min(Re_cutoff,Re_chord);

%% Component Form Factor
FF = (1+(0.6./mean_x_c).*(mean_t_c) + 100*mean_t_c.^4).*(1.34*state.M.^0.18.*(cos(geo.dihed)).^0.28);

%% Laminar/Turbulent Drag
Cf_lam = 1.328./(sqrt(Re_min));
Cf_turb = 0.455./(log10(Re_min).^2.58*(1+0.144*state.M^2)^0.65);

%% Combined Drag
Cf_com = tr.*Cf_lam + (1-tr).*Cf_turb; % This is per unit span

Cd0 = sum(Cf_com.*FF*Q.*Swet_c)/ref.S_ref;
d0 = Cd0*(0.5*state.rho*state.AS^2*ref.S_ref);

Cd0_local =(Cf_com.*FF*Q.*Swet_c)./Sref_c;
d0_local =  Cd0_local.*0.5*state.rho*state.AS^2.*Sref_c;

Param.Cf_lam  = Cf_lam;
Param.Cf_turb = Cf_turb;
Param.tr = tr;
Param.Cd0 = Cd0;
Param.d0 = d0;
Param.Re_chord = Re_chord;
Param.cd0 = Cd0_local;
Param.d0 = d0_local;

end