function [DW2] = fastmega(r1,r2)
%   - Uses the equation 7.37 from Aerodynamics for Engineers by Bertin &
%   Smith which is taken directly from Eqn B-8 of Application of the Vortex
%   Lattice Method to Propellor Performance Analysis by Masquelier (1982).
%
%% First part
% Cross product
% Calculate cross product
perm = [3,1,2];
a = permute(r1,perm);
b = permute(r2,perm);

% F1 = [a(2,:).*b(3,:)-a(3,:).*b(2,:);...
%     a(3,:).*b(1,:)-a(1,:).*b(3,:);...
%     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
% 
% F1 = reshape(F1,size(a));
% 
F1 = [a(2,:,:).*b(3,:,:)-a(3,:,:).*b(2,:,:);...
    a(3,:,:).*b(1,:,:)-a(1,:,:).*b(3,:,:);...
    a(1,:,:).*b(2,:,:)-a(2,:,:).*b(1,:,:)];

% Post-process.
F1 = ipermute(F1,perm);

%%
LF1=(sum(F1.^2,3));
F2 = zeros(size(F1));
F2(:,:,1)=F1(:,:,1)./(LF1);
F2(:,:,2)=F1(:,:,2)./(LF1);
F2(:,:,3)=F1(:,:,3)./(LF1);
%clear('F1')

%% Next part

Lr1=sqrt(sum(r1.^2,3));
Lr2=sqrt(sum(r2.^2,3));

R1 = zeros(size(r1));
R1(:,:,1)=r1(:,:,1)./Lr1;
R1(:,:,2)=r1(:,:,2)./Lr1;
R1(:,:,3)=r1(:,:,3)./Lr1;

R2 = zeros(size(r2));
R2(:,:,1)=r2(:,:,1)./Lr2;
R2(:,:,2)=r2(:,:,2)./Lr2;
R2(:,:,3)=r2(:,:,3)./Lr2;

L1=(R2-R1); %Unit vector
%clear('R1','R2')

%% Third part
R0=(r2-r1);

radial_distance=sqrt((LF1./(sum(R0.^2,3))));

%% combinging 2 and 3 - to give the scalar part of the equation
L2=  R0(:,:,1).*L1(:,:,1)...
    +R0(:,:,2).*L1(:,:,2)...
    +R0(:,:,3).*L1(:,:,3);

%% Downwash
DW = zeros(size(F2));
DW(:,:,1)=F2(:,:,1).*L2;
DW(:,:,2)=F2(:,:,2).*L2;
DW(:,:,3)=F2(:,:,3).*L2;

%Set any downwash terms that are close to the singularity to zero
%   - TODO - Specify the tolerance as an option!
near=1e-2;
DW2 = zeros(size(DW));
DW2(:,:,1)=DW(:,:,1).*(1-(radial_distance<near));
DW2(:,:,2)=DW(:,:,2).*(1-(radial_distance<near));
DW2(:,:,3)=DW(:,:,3).*(1-(radial_distance<near));



end
