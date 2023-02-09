function[panel]=tmesh2(wx,wy,wz,nx,ny,meshtype,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMESH: Essential function for TORNADO						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
% 	tmesh generated vertex points for						%
%	wing division given input arguments						%
%	division corners, numbers of panels in 					%
%	x- and y-direction										%
%															%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Division of Aeronautics		%
%				2000										%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Essential function for TORNADO					%
% Called by:	geometry									%
% Calls:			MATLAB 5.2 std fcns						%
%															%
% Loads: None												%
% Saves: none												%
% Input: wing division corners, nuber of elements in 		%
%			x- n' y-direction								%
% Output:Panel corner coordinates (nx5x3) Matrix			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panel = [];
if nx == 0
    return;
end
a1=[wx(1) wy(1) wz(1)];
b1=[wx(2) wy(2) wz(2)];

b2=[wx(3) wy(3) wz(3)];
a2=[wx(4) wy(4) wz(4)];
percent_cy=(0:ny)./ny;
percent_cx=(0:nx)./nx;

switch meshtype
    case 1
        %Linear lattice, both in x and y
    case 2
        %Linear in x, half cosine in y
        percent_cy=cos(pi/2*(1-percent_cy));
    case 3
        %Cosine in x, half cosine in y
        percent_cx=(cos(pi*(1-percent_cx))+1)/2 ;
        percent_cy=cos(pi/2*(1-percent_cy));
    case 4
        %Cosine in x, cosine in y
        percent_cx=(cos(pi*(1-percent_cx))+1)/2;
        percent_cy=(cos(pi*(1-percent_cy))+1)/2;
    case 5
        %secret, hush-hush ground effect special mission mesh.
        %percent_cx=1.7*percent_cx.^3 - 2.6*percent_cx.^2 + 1.9*percent_cx + 0;
        %percent_cy=1.7*percent_cy.^3 - 2.6*percent_cy.^2 + 1.9*percent_cy + 0;
        
        percent_cx=2.2*percent_cx.^3 - 3.3*percent_cx.^2 + 2.1*percent_cx;
        percent_cy=2.2*percent_cy.^3 - 3.3*percent_cy.^2 + 2.1*percent_cy;
        
        
    otherwise
        disp('NOT IMPLEMENTED')
        %Put new functione here for panel distribution scheme.
end

for i=1:ny+1
    perc_y=percent_cy(i);
    
    c1=b1-a1;
    l1=norm(c1);
    c1_hat=c1./l1;
    d1=(perc_y)*l1*c1_hat;
    m=a1+d1;
    
    c2=b2-a2;
    %l2=norm(c2);
    %c2_hat=c2./l2;
    %d2=(perc_y)*l2*c2_hat;
    d2=(perc_y)*c2;
    
    n=a2+d2;
    
    
    for j=1:nx+1
        
        perc_x=percent_cx(j);
        
        c3=n-m;
        %l3=norm(c3);
        %c3_hat=c3./l3;
        %d3=(perc_x)*l3*c3_hat;
        d3=(perc_x)*c3;
        p=m+d3;
        
        A(i,j,:)=[p];
        
    end
end

t=0;
for i=1:ny
    for j=1:nx
        t=t+1;
        panel(t,1,:)=A(i,j,:);
        panel(t,2,:)=A(i+1,j,:);
        panel(t,3,:)=A(i+1,j+1,:);
        panel(t,4,:)=A(i,j+1,:);
        panel(t,5,:)=A(i,j,:);
    end
end

if (b<0)
    panel(1:end,:,:) = panel(1:end,[2 1 4 3 2],:);
end


end%function
