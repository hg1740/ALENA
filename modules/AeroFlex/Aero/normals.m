%% Normals

function [N, DN, Camber]=normals(colloc,vortex,C_Slope)
step = size(colloc,1); % This was changed to call on the num of rows DARIO
[d, e, ~] = size(vortex);
a = e/2;
b = a + 1;
N = zeros(d, 3);
DN = zeros(d, 3);
Camber = zeros(d, 1);

for t=1:step	%Looping through panels
    
    alpha=C_Slope(t);
    
    for s=1:3						%Looping Through Dimensions.
        ra(s)=vortex(t,a,s);
        rb(s)=vortex(t,b,s);
        rc(s)=colloc(t,s);
    end
    r0=rb-ra;
    r0(1)=0;                    %fix to get normals to not point the right way
    r1=rc-ra;
    r2=rc-rb;
    n=cross(r1,r2);				%Passus to determine normal
    nl=sqrt(sum((n.^2),2));    %of panel at collocationpoint.
    R = n/nl;							%Normalizing normal.
    R2 = trot3(r0,R,-alpha);		%rotating wha trot
    r0 = r0/norm(r0,2);
    Rotmat = RotVecHinge(r0(1),r0(2),r0(3),-alpha);
    R3 = Rotmat*R';
    N(t,:) = R;
    DN(t,:) = R3'-R; % This works fine for a undeformed mesh
    Camber(t,1) = C_Slope(t);
end
end