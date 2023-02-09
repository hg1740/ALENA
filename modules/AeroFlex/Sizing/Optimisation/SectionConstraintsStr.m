%SECTIONCONSTRAINTSSTR: The following scripts takes a set of optimisation
% variables X with other material and box properties to work out the
% streses in the each of the sections. It returns a set of contraint
% variables (CIN) which should ideally be less than zero. This then
% indicates that the contraints have been satisfied
%
%
% Constraints:
%
%
%
function [CIN, CEQ] = SectionConstraintsStr(X, data)

% --------------------------------
% Set up optimisation constraint variable
% --------------------------------
CIN = [];
CEQ = [];

% --------------------------------
% Recall optimistion variables
% --------------------------------
tskin = X(1:data.Nsec);
tspar = X(data.Nsec+1:2*data.Nsec);
Astr  = X(2*data.Nsec+1:end);

% --------------------------------
% Recover some vital box properties
% --------------------------------
C      = data.c_box;        % box chord
H      = data.h_spar;       % box height
E      = data.E;            % Young's modulus
NS     = data.NS;           % Number of stringers
Rpitch = data.RP;           % Rib pitch
KAs    = data.KAs;          % Stringer factor
d      = sqrt(Astr./KAs);   % Stringer depth
ds     = data.ds;           % Stringer thick/depth ratio
Spitch = C./(NS-1);         % Stringer Pitch 

% --------------------------------
% Recover Material allowables
% --------------------------------
smax  = data.smax;
shmax = data.shmax;

% --------------------------------
% Area around the mean line
% --------------------------------
Omega = (C-tspar).*(H-tskin); 

SF = data.SAF_MARG;

% --------------------------------
% Recover the internal loads of each of the sections
% --------------------------------
Ty = SF*data.Fy_Env;    % Vertical Shear
Tx = SF*data.Fz_Env;    % Horizontal Shear
Mx = SF*data.Mz_Env;    % Out of plane bending moment
My = SF*data.My_Env;    % In plane bending moment
Mt = SF*data.Mx_Env;    % Torsional moment

% --------------------------------
% Out of Plane Second moment of area
% --------------------------------
% stringer moment of inertia
%I1st = 1/12*ds*d.^4 + 2*(1/12*ds*d.^4 + ds*d.^4*(1+2*ds+ds^2)) + Astr.*(H/2-tskin - (d/2 + ds*d)).^2;
I1st = d.^4*((7/12)*ds + ds^2 + (8/12)*ds^3) + Astr.*(H/2 - tskin - (d/2 + ds*d)).^2;

% Skin/spar moment of inertia
I1sk = (1/12)*(C.*H.^3 - (C - 2*tspar).*(H-2*tskin).^3);

% overall moment of inertia
I1 = 2*(I1st.*NS) + I1sk;

% --------------------------------
% For-aft Second moment of area
% --------------------------------
% stringer moment of inertia
%I2st = 2*d.^4*ds*(1/12 + 1/4*ds^2 + ds^2/24);
I2st = d.^4*((8/12)*ds - ds^2 + (7/12)*ds^3);

% stringer moment of inertia relative to centroid
I2stp = 2*(NS.*I2st + Astr.*C.^2.*NS.*(-0.25 + (NS + 0.5)./(3*(NS+1))));

% Skin/spar moment of inertia
I2sk = (1/12)*(H.*C.^3 - (H - 2*tskin).*(C- 2*tspar).^3);

% overall moment of inertia
I2 = I2stp + I2sk;

% --------------------------------
% Bending Stress (Calculation)
% --------------------------------
% Skin stress due to bending
sigma_sk = (H/2 - tskin/2).*abs(Mx)./I1;

% stringer stress due to bending
sigma_st = (H/2 - tskin - 0.62*d).*abs(Mx)./I1;

% Spar stress due to bending
sigma_sp = (C/2 - tspar/2).*abs(My)./I2;

%sigma_VM = sqrt(sigma_sk.^2 + 3*tau_skin.^2);

if (abs(sigma_sk)<eps) 
  sigma_sk = eps;
  sigma_st = eps;
end

% --------------------------------
% Bending Stress (Constraint)
% --------------------------------
% Skin bending
Skb = sigma_sk/smax - 1; % Allowable stress due to bending

% Spar bending
Spb = sigma_sp/smax - 1; % Allowable stress due to bending

% Stringer Bending
Stb = sigma_st/smax - 1;

%% Shear Flow
% --------------------------------
% Note here that I'm assuming 3 stringers ... this can be extrapolated
% easily to an arbitrary number
% 7    8    9    10   1
% o----o----o----o----o
% |       |_          |
% |                   |
% o----o----o----o----o
% 6    5    4    3    2
% 
% Calculate the Boom Areas
%
A1 = tspar.*H/6 + Spitch.*tskin./2;
A2 = Spitch.*tskin + Astr.*((H-d*(1+2*ds))./H).^2;
%
% Calculate the Moment of inertia
%
Ix = A1.*H.^2 + NS.*A2.*H.^2/2;
% Iyy ?
% Check:
% http://sydney.edu.au/engineering/aeromech/structures/acs1-p61.html
%  Boom Area  z-Z          Dq             q0 = sum(Dq)
% --------------------------------------------------------------------
% 1   A1      H/2     Ty*A1*H/2Ix          Ty*A1*H/2Ix
% 2   A1      -H/2    -Ty*A1*H/2Ix              0
% 3   A2      -H/2    -Ty*A2*H/2Ix        -Ty*A2*H/2Ix
% 4   A2      -H/2    -Ty*A2*H/2Ix        -2*Ty*A2*H/2Ix 
% 5   A2      -H/2    -Ty*A2*H/2Ix        -3*Ty*A2*H/2Ix 
% 6   A1      -H/2    -Ty*A1*H/2Ix  -3*Ty*A2*H/2Ix - Ty*A1*H/2Ix 
% 7   A1      H/2     Ty*A1*H/2Ix         -3*Ty*A2*H/2Ix
% 8   A2      H/2     Ty*A2*H/2Ix         -2*Ty*A2*H/2Ix
% 9   A2      H/2     Ty*A2*H/2Ix           -Ty*A2*H/2Ix
% 10  A2      H/2     Ty*A2*H/2Ix                0

% if number of stringers = NS;
% then:
% the stringer before the spar q0 = -NS*Ty*A2*H/2Ix;

% Web No. 2A(Internal)             q0*2A      
% --------------------------------------------------------------------
%   1-2      H*C/2               Ty*A1*C*H^2/4Ix
%   2-3    H*C/(2*(Ns+1))                0   
%   3-4    H*C/(2*(Ns+1))      -Ty*A2*C*H^2/4Ix(Ns+1)   
%   4-5    H*C/(2*(Ns+1))      -2*Ty*A2*C*H^2/4Ix(Ns+1)  
%   5-6    H*C/(2*(Ns+1))      -3*Ty*A2*C*H^2/4Ix(Ns+1)  
%   6-7      H*C/2        -3*Ty*A2*C*H^2/4Ix - Ty*A1*C*H^2/4Ix  
%   7-8    H*C/(2*(Ns+1))      -3*Ty*A2*C*H^2/4Ix(Ns+1) 
%   8-9    H*C/(2*(Ns+1))      -2*Ty*A2*C*H^2/4Ix(Ns+1) 
%   9-10   H*C/(2*(Ns+1))       -Ty*A2*C*H^2/4Ix(Ns+1)
%   10-1   H*C/(2*(Ns+1))                 0

% Moment due to cut = -2sum((1:NS))*Ty*A2*C*H^2/4Ix(Ns+1) - NS*Ty*A2*C*H^2/4Ix
Mq0 = Ty.*A2.*C.*H.^2./(4*Ix) .* (-2*sum((1:NS))./(NS+1) - NS); 
%
% Calculate the torque inbalance
Myt = -Mt - Mq0;
%
% qt = Myt./(2*Omega)
%
% Final Shear flow:
%
% Web No.                      q                        
% -----------------------------------------------------------
%   1-2           Ty.*A1.*H./(2*Ix) + Myt./(2*Omega)   
%   2-3                 Myt./(2*Omega)                   
%   3-4           -Ty.*A2.*H./(2*Ix) + Myt./(2*Omega)      
%   4-5           -Ty.*A2.*H./Ix + Myt./(2*Omega)     
%   5-6     -Ty.*A2.*H./Ix - Ty.*A1.*H./(2*Ix)  + Myt./(2*Omega) 
%   6-7           -Ty.*A2.*H./Ix + Myt./(2*Omega)
%   7-8          -Ty.*A2.*H./(2*Ix)  + Myt./(2*Omega)    
%   8-1                  Myt./(2*Omega)          
%
% 
% Shear flow: Spar (1-2) or (5-6)
% Shear flow: Skin (2-3) or (3-4) or (4-5) or (6-7) or (7-8) or (8-1)

% tau_skN = abs(-Ty.*A2.*H./(2*Ix) + Myt./(2*Omega))./tskin;
% tau_spN = abs(Ty.*A1.*H./(2*Ix) + Myt./(2*Omega))./tspar;

% % --------------------------------
% % Shear Stress (Calculation)
% % --------------------------------
% % % Skin shear stress (Neo)
tau_skN = abs(Mt)./(2*Omega.*tskin) + abs(Tx)./(2*C.*tskin) + 0.5 *abs(Ty)./(2*H.*tskin);

% Spar shear stress (Neo)
tau_spN = abs(Mt)./(2*Omega.*tspar) + abs(Ty)./(2*H.*tspar) + 0.25*abs(Tx)./(2*C.*tspar);

% --------------------------------
% Shear Stress (Constraint)
% --------------------------------
% Skin shear
Sks = tau_skN/shmax - 1;

% Spar shear
Sps = tau_spN/shmax - 1;

% Von Mises Criterion:
Fmat_t = sqrt(sigma_sk.^2 + 3*tau_skN.^2)./smax - 1;
%
% %-------------------------------------------------------------------------------
% % NeoCASS
% %-------------------------------------------------------------------------------
% % spar
% %
% % web compression
% %
% tb = tspar./min(Rpitch, H);
% ab = max(Rpitch,H)./min(Rpitch,H);
% ks = -0.138*(ab.^3)+1.4076*(ab.^2)-4.6814*ab+10.467;
% tau = abs(Mt)./(2*Omega.*tspar) + abs(Ty)./(2*H.*tspar) + 0.25*abs(Tx)./(2*C.*tspar);
% Rsw = (abs(tau) .* 12 * (1-0.33^2))./((pi^2)*ks*E.*(tb.^2));
% %
% % web bending
% %
% ab = Rpitch./Spitch;
% tb = tspar./H;
% kb = 0.0012*(ab.^4)-0.0646*(ab.^3)+1.2644*(ab.^2)-10.08*ab+61.983;
% sol = (sigma_sk * 12 * (1-0.33^2))./((pi^2)*kb*E.*(tb.^2));
% check = (tanh(ab-3)+1)/2;
% Rbw = sol.*check;
% %
% % combined 
% %
% Msw = 1-1./(sqrt(Rsw.^2+Rbw.^2) + 0.001);
%-------------------------------------------------------------------------------

%% Buckling Constraints - I need to add skin stringer buckling here!
%ab  = max(Rpitch./Spitch,Spitch./Rpitch);
b   = min(Rpitch,Spitch);
tb  = tskin./b;

% Skin Buckling due to bending
kc  = 4; % simply supported edges ... assumes that the a/b is >> 1
K_C = kc*(pi^2)/(12 * (1-0.33^2));
sigma_cr_sk_bend = K_C.*E.*(tb.^2);

Rbp = (sigma_sk ./sigma_cr_sk_bend); % sigma over sigma critical

% Skin Buckling due to shear
ks  = 5.6; % simply supported edges ... assumes that the a/b is >> 1
K_S = ks*(pi^2)/(12 * (1-0.33^2));
sigma_cr_sk_shear = K_S*E.*(tb.^2);
Rsp = (tau_skN ./sigma_cr_sk_shear); % sigma over sigma critical

% Use the principle stress formulation
Msp = 1-2./( Rbp + sqrt(Rbp.^2+4*(Rsp.^2)) + 0.001);

%% Stringer Buckling
% Stringer local buckling (crippling stress)
% Create a for the nondimensional crippling stress curves for extrusions
% (Pg. 444 Airframe stress analysis and sizing - Niu);

% Crippling Stringer stress (sigma_cc)
A_free = -0.7885215;
B_free = 0.6193987;

A_fixed = -0.804637;
B_fixed = 1.211687474;

crip_x = sqrt(smax./E)*(1/ds);

sigma_co_nd = 1.45; % cut-off

sigma_cc1_nd = B_free*crip_x.^A_free;
sigma_cc2_nd = B_fixed*crip_x.^A_fixed;

sigma_cc1_nd = min(sigma_cc1_nd,sigma_co_nd);
sigma_cc2_nd = min(sigma_cc2_nd,sigma_co_nd);

sigma_cc_nd = (ds*d.^2*(2*sigma_cc1_nd + sigma_cc2_nd))./Astr;

sigma_cc = sigma_cc_nd*smax;

Bcr = sigma_st./sigma_cc - 1;

%% Column buckling of the skin-stringer panel
% Johnson Euler Method
% Linearly interpolate buckling coefficient for skin-stringer panel
k_c_sktr = 0.0386*Spitch./tskin + 2.0771;

k_c_sktr(Spitch./tskin < 40) =  3.62;
k_c_sktr(Spitch./tskin > 110) =  6.32;

Eff_width = tskin.*sqrt(k_c_sktr.*E./sigma_cc);

% Update the effective width so that it is bounded by the stringer pitch
Eff_width(Eff_width>Spitch) = Spitch(Eff_width>Spitch);

Eff_area  = Eff_width.*tskin;

L = Rpitch;
fix_coeff = 1.5; % As recommended by Nui for a typical wing section

% Centroid w.r.t top of skin -> Y = Sum(Ay)/A
y_cent = (Eff_width.*tskin.^2/2 + 3*ds*d.^2.*(tskin + ds*d + d/2))./(Astr+Eff_area);

% Sum up the second moment of areas of each section
sumI0  = 2*(1/12)*(d.*(ds*d).^3) + (1/12)*(ds*d.*d.^3) + (1/12)*Eff_width.*tskin.^3;

% Sum up -> Sum(Ay^2)
sumAyy = Eff_area.*(tskin/2).^2 + ds*d.^2.*((tskin + 0.5*ds*d).^2 +(tskin + ds*d + d/2).^2 +(tskin + 3*ds*d/2 + d).^2);

% AY^2
sumycA = y_cent.^2.*(Astr+Eff_area);

% Second moment of area
I = sumI0 + sumAyy - sumycA;

% Radius of gyration
radius = sqrt(I./(Astr+Eff_area));

% Segment slenderness ratio: K = L/(rho*sqrt(c))
K = L./(radius.*sqrt(fix_coeff));
K_crit = pi*sqrt(2*E./sigma_cc);

% Johnson Euler Stress Niu p.125 
Euler_co = K > K_crit;
John_co  = K <= K_crit;

% Allowable column stress based on Johnson Euler Formula
sigma_John  = sigma_cc.*(1-sigma_cc.*K.^2./(4*pi^2*E));

sigma_Euler = pi^2*E./(K.^2);

sigma_all = sigma_John.*John_co + sigma_Euler.*Euler_co;

Bglob = sigma_st./sigma_all - 1;
%Bglob = 1 - sigma_all./sigma_st;

% ---------------------------------------------
% t_str = ds*d;
%FEAS1 = 0.15*tskin - t_str;
%FEAS2 = t_str - 0.3*tskin;

% Farrar's efficiency factor

% As_bt = Astr./(Spitch.*tskin);
% ts_t = ds.*d./tskin;

% t     = ds.* d;
% FEAS1 = tskin - 1.25* t;
% FEAS2 = 0.75*t - tskin;

% FEAS1 = (ds.*d./tskin - 0.20001)/100;
% FEAS2 = (0.2 - ds.*d./tskin)*10000;
% ---------------------------------------------

%% Store Constraints
CIN = [Fmat_t;Sps;Msp;Bglob];%;instglob];%;FEAS1;FEAS2];
% skin_min = tskin == data.tskin_min;
% 
% for i = 1:data.Iteration
%     legendstring{1,i} = ['Iteration # ' num2str(i)];
% end
% 
% if data.PlotStress == 1
%     lcol = data.Colours(data.Iteration,:);
%     Cons_fig = figure(201);hold on;
%     set(Cons_fig, 'Position', [300, 50, 1000, 400]);
%     
%     subplot(2,3,1);plot(sigma_sk/smax,'color',lcol);hold on;
%     legend(legendstring,'Location','Best');
%     ylabel('Skin stress due to bending / Allowable bending stress');
%     subplot(2,3,2);plot(sigma_sp/smax,'color',lcol);hold on;
%     ylabel('Spar stress due to bending / Allowable bending stress');
%     subplot(2,3,3);plot(tau_skN/shmax,'color',lcol);hold on;
%     ylabel('Skin stress due to shear / Allowable shear stress');
%     subplot(2,3,4);plot(tau_spN/shmax,'color',lcol);hold on;
%     ylabel('Spar stress due to shear / Allowable shear stress');
%     subplot(2,3,5);plot(Fmat_t + 1,'color',lcol);hold on;
%     ylabel('Von Mises Stress/Allowable bending stress');
%     subplot(2,3,6);plot(Msp + 1,'color',lcol);hold on;
%     ylabel('Principle stress criteria for skin buckling');
%     
% end
% 
% if data.CombinedPlot == 1
%     Cons_fig = figure(202);
%     set(Cons_fig, 'Position', [300, 300, 600, 400],'defaultAxesFontSize',16);
%     hold off;
%     plot(data.y_mbox,Fmat_t + 1,'k-');hold on;
%     plot(data.y_mbox,Msp + 1,'r-');hold on;
%     plot(data.y_mbox,skin_min,'g-');hold on;
%     plot(data.y_mbox,Sps + 1,'b-');hold on;
%     plot(data.y_mbox,Bcr + 1,'m-');hold on;
%     plot(data.y_mbox,Bglob + 1,'c-');hold on;
%     %plot(data.y_mbox,1.24*d./tskin,'r--');hold on;
%     %plot(1./(abs(FEAS1) + 1),'y-');hold on;
%     ylim([0,1.2]);
%     ylabel('Constraint Activation');
%     xlabel('Spanwise location [m]');
%     legend('Material strength','Panel Buckling','Minimum Skin Thickness','Spar Shear','Stringer Crippling','Skin-Stringer Buckling','Stringer depth/Skin thickness','Location','Best');
% end

end


% figure;
% for idx = 1:32;
%     shear  = [];
%     
%     shear(:,1) =  Ty(idx)*A1(idx)*H(idx)/(2*Ix(idx)) + Myt(idx)/(2*Omega(idx));
%     shear(:,2) = Myt(idx)/(2*Omega(idx));
%     for i = 1:NS(idx)
%         shear(:,2+i) = -i*Ty(idx)*A2(idx)*H(idx)/(2*Ix(idx)) + Myt(idx)/(2*Omega(idx));
%     end
%     shear(:,2+NS(idx)+1) = -NS(idx)*Ty(idx)*A2(idx)*H(idx)/(2*Ix(idx)) - Ty(idx)*A1(idx)*H(idx)/(2*Ix(idx)) + Myt(idx)/(2*Omega(idx));
%     shear(:,2+NS(idx)+2) = -NS(idx)*Ty(idx)*A2(idx)*H(idx)/(2*Ix(idx)) + Myt(idx)/(2*Omega(idx));
%     for i = 1:NS(idx)
%         shear(:,2+NS(idx)+2+i) = -(NS(idx)-i)*Ty(idx)*A2(idx)*H(idx)/(2*Ix(idx)) + Myt(idx)/(2*Omega(idx));
%     end
%     
%     hold on;plot(shear)
% end
