function [lattice]=setrudder3(rudder,deflection,lattice,geo)

[I K]=find(geo.flapped');
try
    wing=K(rudder);
    division=I(rudder);
catch
    terror(2)
    return
end

[q1 q2 q3]=size(lattice.VORTEX);

if q2==8
    tempV1=lattice.VORTEX(:,1,:);
    tempV2=lattice.VORTEX(:,8,:);
    lattice.VORTEX=lattice.VORTEX(:,2:7,:);
end

fsym=geo.fsym(wing,division);

mp=3;

t=1;
r=0;
[q6 q7]=size(geo.nx);
nr=((geo.nx+geo.fnx).*geo.ny).*((ones(q6,q7)+(geo.symetric'*ones(1,q7))));
[q4 q5]=size(nr);
for i=1:q4
    for j=1:q5
        if geo.flapped(i,j)==1
            r=r+1;
        end
        if r<rudder
            t=t+nr(i,j);
        end
    end
end

nx=geo.nx(wing,division);
ny=geo.ny(wing,division);
fnx=geo.fnx(wing,division);

a1=[lattice.XYZ(t+nx,1,1) lattice.XYZ(t+nx,1,2) lattice.XYZ(t+nx,1,3)];
b1=[lattice.XYZ(t+nx,2,1) lattice.XYZ(t+nx,2,2) lattice.XYZ(t+nx,2,3)];

a2=[lattice.XYZ(t+nx,2,1) -lattice.XYZ(t+nx,2,2) lattice.XYZ(t+nx,2,3)];
b2=[lattice.XYZ(t+nx,1,1) -lattice.XYZ(t+nx,1,2) lattice.XYZ(t+nx,1,3)];

h=b1-a1;				%defining hingeline SB-side
h1_hat=h./norm(h); %normalizing hingeline

h2=b2-a2;				%defining hingeline P-side
h2_hat=h2./norm(h2); %normalizing hingeline


s=nx+fnx;

for i=1:(nx+fnx)*ny*(1+geo.symetric(wing));
    %loop for trailing edge points
    rad2=t+i-1;
    
    if rad2 < t+(nx+fnx)*ny; %if wing is symmetric and
        %panel is on the SB-side
        a=a1;
        b=b1;
        h_hat=h1_hat;
        def=deflection;
    else							% if wing is on the P-side
        h_hat=h2_hat;
        a=a2;
        b=b2;
        if fsym==0;				%if flap deflection is anti-
            %symmetric
            def=-deflection;
        else
            def=deflection;
        end
    end
    
    for col=1:5:6
        
        p1(1)=lattice.VORTEX(rad2,col,1);
        p1(2)=lattice.VORTEX(rad2,col,2);
        p1(3)=lattice.VORTEX(rad2,col,3);
        if col<=mp
            r=p1-a;
            p2=trot3(h_hat,r,def);
            
            lattice.VORTEX(rad2,col,1)=p2(1)+a(1);
            lattice.VORTEX(rad2,col,2)=p2(2)+a(2);
            lattice.VORTEX(rad2,col,3)=p2(3)+a(3);
        else
            r=p1-b;
            p2=trot3(h_hat,r,def);
            
            lattice.VORTEX(rad2,col,1)=p2(1)+b(1);
            lattice.VORTEX(rad2,col,2)=p2(2)+b(2);
            lattice.VORTEX(rad2,col,3)=p2(3)+b(3);
        end
    end
end

for i=s:s:s*ny*(1+geo.symetric(wing))
    %stepping through number of strips
    for j=0:fnx-1 %stepping through number of flappanels
        ii=i-fnx;
        rad1=(t+ii+j);
        
        if rad1 < t+(nx+fnx)*ny; %if wing is symmetric and
            %panel is on the SB-side
            a=a1;
            b=b1;
            h_hat=h1_hat;
            def=deflection;
        else							% if wing is on the P-side
            h_hat=h2_hat;
            a=a2;
            b=b2;
            if fsym==0;				%if flap deflection is anti-
                %symmetric
                def=-deflection;
            else
                def=deflection;
            end
            
        end
        
        
        
        for k=0:3 %Vortex loop
            col=(k+mp-1);
            p1(1)=lattice.VORTEX(rad1,col,1);
            p1(2)=lattice.VORTEX(rad1,col,2);
            p1(3)=lattice.VORTEX(rad1,col,3);
            if col<=mp
                r=p1-a;
                p2=trot3(h_hat,r,def);
                
                lattice.VORTEX(rad1,col,1)=p2(1)+a(1);
                lattice.VORTEX(rad1,col,2)=p2(2)+a(2);
                lattice.VORTEX(rad1,col,3)=p2(3)+a(3);
            else
                r=p1-b;
                p2=trot3(h_hat,r,def);
                
                lattice.VORTEX(rad1,col,1)=p2(1)+b(1);
                lattice.VORTEX(rad1,col,2)=p2(2)+b(2);
                lattice.VORTEX(rad1,col,3)=p2(3)+b(3);
            end
            
        end
        
        
        %collocarion point rotation
        p1(1)=lattice.COLLOC(rad1,1);
        p1(2)=lattice.COLLOC(rad1,2);
        p1(3)=lattice.COLLOC(rad1,3);
        
        c=(a+b)./2;
        r=p1-c;
        p2=trot3(h_hat,r,def);
        
        lattice.COLLOC(rad1,1)=p2(1)+c(1);
        lattice.COLLOC(rad1,2)=p2(2)+c(2);
        lattice.COLLOC(rad1,3)=p2(3)+c(3);
        
        %Normals rotation
        
        p1(1)=lattice.N(rad1,1);
        p1(2)=lattice.N(rad1,2);
        p1(3)=lattice.N(rad1,3);
        
        c=(a+b)./2;
        r=p1;
        p2=trot3(h_hat,r,def);
        
        lattice.N(rad1,1)=p2(1);
        lattice.N(rad1,2)=p2(2);
        lattice.N(rad1,3)=p2(3);
        
        for k=0:4 %panelcoords
            col=(k+1);
            p1(1)=lattice.XYZ(rad1,col,1);
            p1(2)=lattice.XYZ(rad1,col,2);
            p1(3)=lattice.XYZ(rad1,col,3);
            %disp('************')
            if col<=1;
                r=p1-a;
                p2=trot3(h_hat,r,def);
                
                lattice.XYZ(rad1,col,1)=p2(1)+a(1);
                lattice.XYZ(rad1,col,2)=p2(2)+a(2);
                lattice.XYZ(rad1,col,3)=p2(3)+a(3);
                
            elseif col<=3
                r=p1-b;
                p2=trot3(h_hat,r,def);
                
                lattice.XYZ(rad1,col,1)=p2(1)+b(1);
                lattice.XYZ(rad1,col,2)=p2(2)+b(2);
                lattice.XYZ(rad1,col,3)=p2(3)+b(3);
            else
                r=p1-a;
                p2=trot3(h_hat,r,def);
                
                lattice.XYZ(rad1,col,1)=p2(1)+a(1);
                lattice.XYZ(rad1,col,2)=p2(2)+a(2);
                lattice.XYZ(rad1,col,3)=p2(3)+a(3);
            end
            
        end
    end
end

if q2==8
    lattice.VORTEX=[tempV1 lattice.VORTEX tempV2];
end
end
