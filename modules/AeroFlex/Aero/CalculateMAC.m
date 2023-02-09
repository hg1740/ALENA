%%CALCULATEMAC
% 
% 
%   Author: Dario Calderon 
function [mac, mac_LE, mac_AC] = CalculateMAC(AERO,idx)

if nargin < 2
    np = length(AERO.ID);
    idx = 1:np;
else
    np = length(idx);
end

geo = AERO.geo;

MACP  = zeros(np,1);
AREAP = zeros(np,1);
XLE   = zeros(np,1);
XAC   = zeros(np,1);

for i = idx
    [MACP(i), AREAP(i), XLE(i), XAC(i)] = MAC_patch(geo.c(i), geo.T(i), geo.SW(i), ...
        geo.b(i));
    XAC(i) = XAC(i) + geo.startx(i);
    XLE(i) = XLE(i) + geo.startx(i);
end

% Area weighted results
AREATOT = sum(abs(AREAP));
mac     = dot(MACP, abs(AREAP))/ AREATOT;
mac_LE  = dot(XLE, abs(AREAP)) / AREATOT;
mac_AC  = dot(XAC, abs(AREAP)) / AREATOT;

end
%-------------------------------------------------------
function [MAC, A, xLE, xAC] = MAC_patch(CROOT, TP, SW, SPAN)

% patch MAC
MAC = (2/3)*(1+TP+TP^2)*CROOT/(1+TP);

% patch Area
A = CROOT*(1+TP)*SPAN/2;

if TP == 1
    eta = 0.5;
else
    eta = interp1([CROOT,TP*CROOT],[0,1],MAC);
end

% try
%     eta = interp1([CROOT,TP*CROOT],[0,1],MAC);
% catch
%     eta = 0.5;
% end

xAC = tan(SW)*eta*SPAN + 0.25*CROOT;
xLE = xAC - 0.25*MAC;
end