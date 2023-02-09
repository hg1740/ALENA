function[x,y,z]=drawhinge(wx,wy,wz,fc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAWHINGE: subsidary function to TORNADO	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that draws the hinge line on		%
% on a wing division. It also returns the 	%
% coordinates on the foremost flap corners	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by: Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: WX,WY,WZ = wing cornerpoint coor-
%   		dinates.
%			fc is the percentage of total chord
%			built up by the flap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: 	graph (in figure (2))
%				flap cornerpoint coor-
%				dinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(fc);
    
    x=[];
    y=[];
    z=[];
    
else
    
    for i=0:1
        
        a=[wx(1+i) wy(1+i) wz(1+i)];
        b=[wx(4-i) wy(4-i) wz(4-i)];
        
        c=b-a;
        l=norm(c);
        c_hat=c./l;
        d=(1-fc(i+1))*l*c_hat;
        
        r=a+d;
        
        R1(i+1,:)=[r];
        R2(i+1,:)=[r];
    end
    
    x=[R1(:,1)'];
    y=[R1(:,2)'];
    z=[R1(:,3)'];
    
end
end