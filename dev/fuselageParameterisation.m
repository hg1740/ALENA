function fuselageParameterisation
%fuselageParameterisation Creates a parameterisation of a fuselage using
%the method described in the PhD thesis - "Geometric Parameterisation and 
%Aerodynamic Shape Optimisation", Feng Zhu.

close all

%Set up figure
hF  = figure('Name', 'Fuselage geometry');
hAx = axes('Parent', hF, 'NextPlot', 'add');
xlabel(hAx, 'X [m]');
ylabel(hAx, 'Y [m]');
zlabel(hAx, 'Z [m]');
axis(hAx  , 'equal');

%% Cylindrical Fuselage

%Parameters
length   = 25;
width    = 4; 
height   = 3;
nPointsX = 20;
nPointsY = 25;

%Calculate the 3D coordinates of the elliptical midsection
[X, Y, Z] = getMidSectionCoords(length, width, height, nPointsX, nPointsY);

%Plot
plot3(hAx, X, Y, Z, 'ro')
surf(hAx, X, Y, Z);

%% Nose section

%Parameters
noseLength = length / 5;

ksi = linspace(0, 1, nPointsX);

coords = [ ...
    0, 0 ; ...
    3, 4 ; ...
    6, 4 ; ...
    9, 0];


%Bernstein polynomial
% nPow = 30;
nPow = size(coords, 1);
pow = (1 : nPow)';

Ki = factorial(nPow) ./ (factorial(pow) .* factorial(nPow - pow));
% Si = ksi ^ (pow') .* (1 - ksi) .^ (nPow - pow);

% bb = bsxfun(@power, ksi, pow);

%'ksi' raised to powers of r
ksiR  = arrayfun(@(r) ksi .^ r, pow, 'Unif', false);
ksiR  = vertcat(ksiR{:});

%(1 - ksi) raise to powers of (n - r)
ksiNR = arrayfun(@(nr) (1 - ksi) .^ nr, nPow - pow, 'Unif', false);
ksiNR = vertcat(ksiNR{:});

Si = Ki .* ksiR .* ksiNR;


hF  = figure('Name', 'Bernstein Polynomials');
hAx = axes('Parent', hF, 'NextPlot', 'add');
hL  = plot(ksi, Si);
str = arrayfun(@(n) sprintf('n = %i', n), pow, 'Unif', false);
legend(hAx, hL, str, 'Location', 'NorthEastOutside');
xlabel('\ksi [-]');
ylabel('Si [-]');
title(hAx, sprintf('Bernestien polynomials in increasing power up to n = %i', nPow));

S = sum(Si);

end

function [X, Y, Z] = getMidSectionCoords(length, width, height, nPointsX, nPointsY)
%getMidSectionCoords Returns the 3D coordinates of the elliptical fuselage
%midsection based on the desired length, width and height of the
%midsection.

if nargin < 4
    nPointsX = 20;
end
if nargin < 5
    nPointsX = 25;
end

%Non-dimensionsal vectors
eta = awi.model.LiftingSurface.cosspace(0, 1, nPointsY, 'bRow', true);
ksi = linspace(0, 1, nPointsX);

%Calculate upper surface coordinates using Equations 4.37, 4.38, 4.39 of
%Ref. [1].
X  = ksi .* length;
Yu = eta .* width;
Zu = 2 .* height .* eta.^0.5 .* (1 - eta) .^ 0.5;

%Lower surface is just mirror of the upper surface
Yl = fliplr(Yu);
Zl = fliplr(-Zu);

%Combine & stack
Y = [Yu, Yl];
Z = [Zu, Zl];
X = repmat([X, fliplr(X)]', [1, nPointsY * 2]); 
Y = repmat(Y, [nPointsX * 2, 1]);
Z = repmat(Z, [nPointsX * 2, 1]);

end

