function [C, Vor, N, DN, camber, twist,P, ndof, cdof, hinge] = VortCollNorm(fnx, ny, nx, fsym, fc, flapped, TW, foil, T, SW, c, dihed, b, sym, sx, sy, sz, meshtype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GEOMETRY: Essential function for TORNADO				 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the position of vortex-collocation-normals	 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Tomas Melin, KTH, Division of Aeronautics	 %
%				copyright 2000							 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTEXT:	Subsidary function for TORNADO				 %
% Called by:	setup									 %
% Calls:			MATLAB 5.2 std fcns, tmesh, drawhinge%
%					slope, normals						 %
% Loads:	none										 %
% Saves: none											 %
% Input: wing and division number						 %
% Output:coordinades for collocationpoints, vorticies and%
% 			Normals										 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TEP=[];
TEP1=[];
TEP2=[];
INF=[];
INF1=[];
INF2=[];

ox = sx;
oy = sy;
oz = sz;
neqns = (nx+fnx) * ny;

if sym==1
    ndof = 2 * neqns;
else
    ndof = neqns;
end

dx = (c*(1-fc(1))/nx); % panel chord at root

if flapped==1
    
    fdx = (c*fc(1)/fnx);
    
else
    
    fdx = 0;
    
end

a1 = ones(nx,1)*dx; %fixed panel chords
a2 = ones(fnx,1)*fdx; % control panel chords
dr=[a1' a2'];

%%%%%%%%%%%%%%%%%%%%%%%
%Calculates geometry, collocationpoints, panels and vortecies for a flat quad
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Plotting planform
%%%%%%%%%%%%%%%%%%

lem(1)=0.25*c;
lem(2)=0.25*T*c;
lem(3)=-0.75*T*c;
lem(4)=-0.75*c;

% Calculate the change due to preodminantly due to twist
% DX =[(1-cos(TW(1,1,1)))*cos(SW) (1-cos(TW(1,1,2)))*cos(SW) (1-cos(TW(1,1,2)))*cos(SW) (1-cos(TW(1,1,1)))*cos(SW)].*lem;
% DY =-[sin(TW(1,1,1))*sin(dihed)*cos(SW) sin(TW(1,1,2))*sin(dihed)*cos(SW) sin(TW(1,1,2))*sin(dihed)*cos(SW) sin(TW(1,1,1))*sin(dihed)*cos(SW)].*lem;
% DZ =[sin(TW(1,1,1))*cos(dihed) sin(TW(1,1,2))*cos(dihed) sin(TW(1,1,2))*cos(dihed) sin(TW(1,1,1))*cos(dihed)].*lem;
%

%
% wingx =[0 0.25*c+b*tan(SW)-0.25*T*c 0.25*c+b*tan(SW)+0.75*T*c c] + ox + DX;
% wingy =[0 b*cos(dihed) b*cos(dihed) 0] + oy + DY;
% wingz =[0 b*sin(dihed) b*sin(dihed) 0] + oz + DZ;

DX =[(1-cos(TW(1,1,1))) (1-cos(TW(1,1,2))) (1-cos(TW(1,1,2))) (1-cos(TW(1,1,1)))].*lem;
DY =-[0 0 0 0].*lem;
DZ =[sin(TW(1,1,1)) sin(TW(1,1,2)) sin(TW(1,1,2)) sin(TW(1,1,1))].*lem;

wingx = [0 0.25*c+b*tan(SW)-0.25*T*c 0.25*c+b*tan(SW)+0.75*T*c c] + ox;

if dihed == pi/2
    wingy = [0 b*cos(dihed) b*cos(dihed) 0] + oy;
    wingz = [0 b*sin(dihed) b*sin(dihed) 0] + oz;
else
    wingy =[0 b b 0] + oy;
    wingz = [0 b*tan(dihed) b*tan(dihed) 0] + oz;
end

hinge = zeros(1,2,3);
cdof = [];
%%%%%%%%%%%%%%%%%
%Plotting hinge %
%%%%%%%%%%%%%%%%%
if flapped==1
    
    [flapx flapy flapz] = drawhinge(wingx, wingy, wingz, fc);
    
    hinge(1,1,1) = flapx(1);
    hinge(1,2,1) = flapx(2);
    
    hinge(1,1,2) = flapy(1);
    hinge(1,2,2) = flapy(2);
    
    hinge(1,1,3) = flapz(1);
    hinge(1,2,3) = flapz(2);
    
    offset = 0;
    
    for i = 1:ny
        
        cdof = [cdof, [(offset + nx +1) : (offset + nx + fnx)]];
        
        offset = offset + nx + fnx;
        
    end
    
    %  if (b<0)
    %    hinge = hinge(1,2:-1:1,:);
    %  end
    
end

if flapped==0
    
    [p] = tmesh2(wingx, wingy, wingz, nx, ny, meshtype, b);
    % panel vertex coordinates ; stripes along flow
    PX(:,:) = p(:,:,1);
    PY(:,:) = p(:,:,2);
    PZ(:,:) = p(:,:,3);
    
else
    
    tempx=wingx(3:4);
    tempy=wingy(3:4);
    tempz=wingz(3:4);
    
    wingx(3:4)=fliplr(flapx(1:2));
    wingy(3:4)=fliplr(flapy(1:2));
    wingz(3:4)=fliplr(flapz(1:2));
    
    flapx(3:4)=tempx;
    flapy(3:4)=tempy;
    flapz(3:4)=tempz;
    
    [p]=tmesh2(wingx,wingy,wingz,nx,ny,meshtype,b);
    [q]=tmesh2(flapx,flapy,flapz,fnx,ny,meshtype,b);
    r=[];
    for i=1:ny
        count1=((1:nx)+(nx*(i-1)));
        count2=(1:fnx)+(fnx*(i-1));
        r=[r;p(count1,:,:);q(count2,:,:)];
    end
    
    PX(:,:)=r(:,:,1);
    PY(:,:)=r(:,:,2);
    PZ(:,:)=r(:,:,3);
end

nx = nx + fnx;

%%%%%%%%%%%%%%%%%%%
%Panel plot.
%Collocation point tensor generation & plot.
%Vortex tensor generation & plot.
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over all panels on quad. Determines panel corners, %
% vortex coo-rds, and collocation coo-rds		             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[EP, Yu, Yl, X_1_S, lemma_1_S_tot] = airfoil_mean_line(foil(1,1,1),0); %element inboard camber slope
[EP, Yu, Yl, X_2_S, lemma_2_S_tot] = airfoil_mean_line(foil(1,1,2),0); %element outboard camber slope

t = 0;
for j=0:(ny-1);
    %
    % determine the following data if control surface fraction is different
    cpi1 = nx*j+1;
    cpi2 = cpi1+nx-1;
    % leading edge midpoint
    refp = 0.5.*([PX(cpi1,1,:), PY(cpi1,1,:), PZ(cpi1,1,:)] + [PX(cpi1,2,:), PY(cpi1,2,:), PZ(cpi1,2,:)]);
    % determine panel strip chord
    mschord = norm(refp - 0.5.*([PX(cpi2,3,:), PY(cpi2,3,:), PZ(cpi2,3,:)] + [PX(cpi2,4,:), PY(cpi2,4,:), PZ(cpi2,4,:)]));
    %
    for i=0:(nx-1);
        t=t+1; % panel counter
        
        px=PX(t,:);
        py=PY(t,:);
        pz=PZ(t,:);
        
        if i==(nx-fnx-1) %if the panel is the rearest chordwise on wing, forward of flap.
            for s=0:(nx-fnx-1);
                HP(t-s,1,:)=[px(4) py(4) pz(4)];
                HP(t-s,2,:)=[px(3) py(3) pz(3)];
                % TEP=Trailing edge points, Vortex points on the trailing edge
                if sym==1
                    %Port side points
                    HP(t-s+neqns,1,:)=[px(3) -py(3) pz(3)];
                    HP(t-s+neqns,2,:)=[px(4) -py(4) pz(4)];
                end
            end
        end
        
        if i==(nx-1);		%if the panel is the rearest chordwise on both wing and flap
            for s=0:(nx-1);
                TEP1(t-s,1,:)=[px(4) py(4) pz(4)];
                TEP1(t-s,2,:)=[px(3) py(3) pz(3)];
                % TEP=Trailing edge points, Vortex points on the trailing edge
                if sym==1
                    %Port side points
                    TEP1(t-s+neqns,1,:)=[px(3) -py(3) pz(3)];
                    TEP1(t-s+neqns,2,:)=[px(4) -py(4) pz(4)];
                end
                
                for u=0:(fnx-1)	%Hinge points for flap (equals trailing points)
                    
                    HP(t-u,1,:)=[px(4) py(4) pz(4)];
                    HP(t-u,2,:)=[px(3) py(3) pz(3)];
                    
                    if sym==1
                        %Port side points
                        HP(t-u+neqns,1,:)=[px(3) -py(3) pz(3)];
                        HP(t-u+neqns,2,:)=[px(4) -py(4) pz(4)];
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Collocation point   %
        % tensor generation   %
        %%%%%%%%%%%%%%%%%%%%%%%
        
        mx=sum(px(1:4))/4;		%panel midpoint x-coord
        my=sum(py(1:4))/4;
        mz=sum(pz(1:4))/4;
        
        bkx=(px(3)+px(4))/2;		%panel rear edge avarage x-coord
        
        C1(t,1)=(mx+bkx)/2;				%SB-Collocation point x-coord.
        C1(t,2)=(py(3)+py(4)+2*my)/4;		%SB-Collocation point y-coord.
        C1(t,3)=(pz(3)+pz(4)+2*mz)/4;		%SB-Collocation point z-coord.
        if sym==1
            C2(t,1)=C1(t,1);					%P-Collpoint x-coord.
            C2(t,2)=-C1(t,2);					%P-Collpoint y-coord.
            C2(t,3)=C1(t,3);					%P-Collpoint z-coord.
        else
            C2=[];
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Vortex tensor generation and plot %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ax=((px(1)+px(4))/2+px(1))*0.5;	%vortex first point
        ay=(3*py(1)+py(4))/4;
        az=(3*pz(1)+pz(4))/4;
        
        bx=((px(2)+px(3))/2+px(2))*0.5;	%vortex second point
        by=(3*py(2)+py(3))/4;
        bz=(3*pz(2)+pz(3))/4;
        
        V1(t,1,1)=ax;
        V1(t,1,2)=ay;
        V1(t,1,3)=az;
        V1(t,2,1)=bx;
        V1(t,2,2)=by;
        V1(t,2,3)=bz;
        
        if sym==1;
            
            V1(t+neqns,1,1)=bx;
            V1(t+neqns,1,2)=-by;
            V1(t+neqns,1,3)=bz;
            V1(t+neqns,2,1)=ax;
            V1(t+neqns,2,2)=-ay;
            V1(t+neqns,2,3)=az;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Passus to compute camber slope at %
        % Station							%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ( fc(1) == fc(2) )
            a3=(sum(dr(1:i+1))-0.25*dr(1+i))/c;
        else
            a3 = norm(refp - C1(t,:)) / mschord;
        end
        %
        lemma_1_S(t) = interp1(X_1_S,lemma_1_S_tot,a3,'pchip'); %element inboard camber slope
        lemma_2_S(t) = interp1(X_2_S,lemma_2_S_tot,a3,'pchip'); %element outboard camber slope
        S(t)=(lemma_1_S(t)*(ny-j)+lemma_2_S(t)*(j))/ny; %average slope for panels on
        
        if sym==1
            S(t+neqns)=S(t);
        end
        
        twist(t,1) = 0.5*(TW(1,1,1) + TW(1,1,2));
    end
end

C=[C1;C2];
V=V1;
Vor=[TEP1(:,1,:) HP(:,1,:) V(:,:,:) HP(:,2,:) TEP1(:,2,:)];

if (b<0)
    Vor = Vor(:,6:-1:1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating normals              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, DN, camber] = normals(C,V,S);
V=Vor;

if sym==1
    PX2(:,1)=PX(:,2);
    PX2(:,2)=PX(:,1);
    PX2(:,3)=PX(:,4);
    PX2(:,4)=PX(:,3);
    PX2(:,5)=PX2(:,1);
    
    PY2(:,1)=PY(:,2);
    PY2(:,2)=PY(:,1);
    PY2(:,3)=PY(:,4);
    PY2(:,4)=PY(:,3);
    PY2(:,5)=PY2(:,1);
    
    PZ2(:,1)=PZ(:,2);
    PZ2(:,2)=PZ(:,1);
    PZ2(:,3)=PZ(:,4);
    PZ2(:,4)=PZ(:,3);
    PZ2(:,5)=PZ2(:,1);
    
    
    PX=[PX;PX2];
    PY=[PY;-PY2];
    PZ=[PZ;PZ2];
end

P(:,:,1)=PX;
P(:,:,2)=PY;
P(:,:,3)=PZ;
end