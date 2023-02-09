
%L0 - smallest gust length
%L1 - largest gust length
%Alt - Alititude in meters
%U - Flight speed in m/s
%t - time vector
function beam_model = GustFamilyGen(beam_model,GustLengths,Alt,U,t,PlotFlag)

% GustLengths = linspace(L0,L1,NumGusts);
NumGusts = length(GustLengths);

A1 = -12/15000;
B1 = 56 * 0.3048;
A2 = -18/35000;
B2 = 0.3048*(44 + 18*15000/35000);

if Alt < 0.3048*15000
    U_ref = A1*Alt + B1;
else
    U_ref = A2*Alt + B2;
end

Fg = 1;

n_steps = length(t);

if length(t) == 1
    t = [];
end
if nargin < 8
    PlotFlag = 0;
end

for gustloop = 1:NumGusts
    H = GustLengths(gustloop);
    Uds = U_ref*Fg*(H/(0.3048*350))^(1/6);
    beam_model.Gust.U_ref(gustloop)  = U_ref;
    beam_model.Gust.ID(gustloop)     = gustloop;
    beam_model.Gust.DIR(gustloop)    = 3;
    beam_model.Gust.Amp(gustloop)    = Uds;
    beam_model.Gust.Tmax(gustloop)   = 2*H/(U);
    beam_model.Gust.X0(gustloop)     = 0;
    beam_model.Gust.funs{gustloop,1} = '1';
    beam_model.Gust.fun{gustloop,1}  = ['0.5*(1-cos(pi*(' num2str(U) '*t/' num2str(H) ')))'];
    Gust_time = zeros(1,n_steps);
    InpGustx = beam_model.Gust.Amp(gustloop)*eval(beam_model.Gust.fun{gustloop,1});
    [~,t0] = min(abs(t-beam_model.Gust.X0(gustloop)/U));
    [~,t1] = min(abs(t-beam_model.Gust.X0(gustloop)/U-beam_model.Gust.Tmax(gustloop)));
    Gust_time(1,t0:t1) = InpGustx(1:t1-t0+1);
    if ~isempty(t) && PlotFlag
        figure(222);
        hold on;
        plot(t,Gust_time,'-','Color',[0.5,0.5,0.5],'LineWidth',2);
    end
end

end