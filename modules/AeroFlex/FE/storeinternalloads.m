%% The following function simply caluclates the average internal load across
% two evaluation points, and stores it into an array.

function BarForces = storeinternalloads(CFORCES)

wz = size(CFORCES,3);

Fx = zeros(wz,1);
Fy = zeros(wz,1);
Fz = zeros(wz,1);
Mx = zeros(wz,1);
My = zeros(wz,1);
Mz = zeros(wz,1);

for j = 1:size(CFORCES,3)
    Fx(j,1)    = 0.5*(CFORCES(1,1,j)+CFORCES(2,1,j));
    Fy(j,1)    = 0.5*(CFORCES(1,2,j)+CFORCES(2,2,j));
    Fz(j,1)    = 0.5*(CFORCES(1,3,j)+CFORCES(2,3,j));
    Mx(j,1)    = 0.5*(CFORCES(1,4,j)+CFORCES(2,4,j));
    My(j,1)    = 0.5*(CFORCES(1,5,j)+CFORCES(2,5,j));
    Mz(j,1)    = 0.5*(CFORCES(1,6,j)+CFORCES(2,6,j));
end

BarForces.Fx(:,1) = Fx;
BarForces.Fy(:,1) = Fy;
BarForces.Fz(:,1) = Fz;
BarForces.Mx(:,1) = Mx;
BarForces.My(:,1) = My;
BarForces.Mz(:,1) = Mz;

end