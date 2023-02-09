%% Beam Reduction
%
% Calculate the equivalent bar properties for a hollow rectangular wing box
%
% Bar_prop = BeamReduction (chord, height, skin thickness, spar thickness)
%
%%%%%%%%%%%%
% Outputs:
%%%%%%%%%%%%
% A = cross-sectional area
% I1
% I2
% J
%
%% Example 1:
% c_box = [4,4,4,4];
% h_spar = [1.5,1,0.75,0.5];
% t_skin = [0.05,0.05,0.05,0.05];
% t_spar = [0.05,0.05,0.05,0.05];
% y_mbox = [2,6,10,14];
% y_lbox = [4,4,4,4];
% Bar_prop = BeamReduction(c_box, h_spar, t_skin, t_spar, y_mbox, y_lbox)
%
% 
% 
%   Author: Dario Calderon 
%%
function data = BeamReduction(data)


c_box   = data.c_box;
h_spar  = data.h_spar; 
t_skin  = data.t_skin;
t_spar  = data.t_spar;
nodex   = data.nodex;
nodey   = data.nodey;
nodez   = data.nodez;
A_str   = data.A_str;
KAs     = data.KAs;
NS      = data.NS;
rho     = data.Rho;

fid = 1;
%% Plot 2D/3D structure
plot3DW = 0;
plot2DS = 0;

NSec = length(c_box);

x_lbox = diff(nodex);
y_lbox = diff(nodey);
z_lbox = diff(nodez);

l_box = sqrt(x_lbox.^2 + y_lbox.^2 + z_lbox.^2);

if plot3DW == 1 && nargin > 4
    fprintf(fid,'\n\tPlotting 3D box\n')
    tic
    plot3Dwingbox(c_box, h_spar, t_skin, t_spar, nodex, nodey, nodez, NS);
    toc
end

if plot2DS == 1
    plot2DSection(c_box, h_spar, t_skin, t_spar);
end

%% Cross-sectional area CHECK!!!!!
A = 2*c_box.*t_skin + 2*(h_spar-2*t_skin).*t_spar + 2*NS.*A_str;

%% Internal area
A_int = c_box.*h_spar - A;

%% Internal Volume
V_int = A_int.*l_box;

d  = sqrt(A_str./KAs);
ds = KAs/3;

%% Out-of-plane second moment of area
% stringer moment of inertia
%I2st = 1/12*ds*d.^4 + 2*(1/12*ds*d.^4 + ds*d.^4*(1+2*ds+ds^2)) + Astr.*(h_spar/2 - t_skin - (d/2 + ds*d)).^2;
I2st = d.^4*((7/12)*ds + ds^2 + (8/12)*ds^3) + A_str.*(h_spar/2 - t_skin - (d/2 + ds*d)).^2;

% Skin/spar moment of inertia
I2sk = (1/12)*(c_box.*h_spar.^3 - (c_box - 2*t_spar).*(h_spar - 2*t_skin).^3);

% overall moment of inertia
I2 = 2*(I2st.*NS) + I2sk;

%% In-plane second moment of area
% stringer moment of inertia
%I1st = 2*d.^4*ds*(1/12 + 1/4*ds^2 + ds^2/24);
I1st = d.^4*((8/12)*ds - ds^2 + (7/12)*ds^3);
% stringer moment of inertia relative to centroid
I1stp = 2*(NS.*I1st + A_str.*c_box.^2.*NS.*(-0.25 + (NS + 0.5)./(3*(NS+1))));

% Skin/spar moment of inertia
I1sk = (1/12)*(h_spar.*c_box.^3 - (h_spar - 2*t_skin).*(c_box - 2*t_spar).^3);

% overall moment of inertia
I1 = I1stp + I1sk;

% Product moment of inertia
I12 = 2*NS.*ds.*d.^2;

%% Torsion
% Torsion constant for skin/spar box
Jsk = 4*(h_spar - t_skin).^2.*(c_box - t_spar).^2./(2*((h_spar - t_skin)./t_spar + (c_box - t_spar)./t_skin));

% Torsion constant of each individual stringer
K1 = 0.3; % Mechanics of Materials P318
Jst = 3*K1*ds.^3.*d.^4;

% Combined Torsion constant
J = Jsk + 2*NS.*Jst;

%% Store properties
data.A  = A;
data.I1 = I2;
data.I2 = I1;
data.I12 = I12;
data.J  = J;
data.A_int = A_int;
data.V_int = V_int;
data.ShearOffsets = zeros(NSec,9);
data.SkinWeight = rho*2*c_box.*t_skin.*l_box;
data.SparWeight = rho*2*(h_spar-2*t_skin).*t_spar.*l_box;
data.StrWeight  = rho*2*NS.*A_str.*l_box;

x1 = [];
x2 = [];
x3 = [];
for i = 1:length(nodex)-1
    x1(i,:) = [1,0,0];
    x2(i,:) = [nodex(i+1)-nodex(i),nodey(i+1)-nodey(i),nodez(i+1)-nodez(i)];
    x2(i,:) = x2(i,:)/norm(x2(i,:));
    x3(i,:) = cross(x1(i,:),x2(i,:));
    x3(i,:) = x3(i,:)/norm(x3(i,:));
    x1(i,:) = cross(x2(i,:),x3(i,:));
    x1(i,:) = x1(i,:)/norm(x1(i,:));
end

data.x1 = x1;
data.x2 = x2;
data.x3 = x3;
end
