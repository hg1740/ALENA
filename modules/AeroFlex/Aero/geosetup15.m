
function [lattice, ref] = geosetup15(fid, geo, ref)

void=0;
try
    geo.meshtype;
catch
    geo.meshtype=ones(size(geo.T));
end
% find surfaces total movable
nt = size(geo.T,2);
TOTAL = zeros(size(geo.T));
for i = 1 :nt
    TOTAL(:,i) = geo.meshtype(:,i)==6;
    if ~isempty(find(TOTAL(:,i)==1))
        geo.meshtype(TOTAL(:,i)==1,i) = ones(size(TOTAL(TOTAL(:,i)==1,i)));
    end
end
clear nt

lattice.COLLOC = [];
lattice.VORTEX = [];
lattice.N      = [];
lattice.DN     = [];
lattice.XYZ    = [];
lattice.DOF    = [];
lattice.Camber = [];
lattice.Twist  = [];

X = [];
Y = [];
Z = [];
S = 1;
Cmac = 0;
CHORDS = [];

loopsperwing = geo.nelem;
noofloops = loopsperwing;
temp = 0;
noofwings = length(loopsperwing');

for s=1:noofwings			% Intermediate variable setuploop
    
    CHORDS(s,1) = geo.c(s);		% calculating chords of first element
    SX(s,1) = geo.startx(s);	% Element apex calculation
    SY(s,1) = geo.starty(s);	% Same-o
    SZ(s,1) = geo.startz(s);    % Same-o
    
end

t=0;							%resetting ticker variable

% draw aerodynamic patches without twist contribution

for s=1:noofwings
    
    for t=1:(noofloops(s))
        %Chord loop, generating chords for wing sections.
        %And startingpoints for partition-quads
        
        CHORDS(s,t+1)=CHORDS(s,t)*geo.T(s,t);	%calculating
        %element root-chord
        
        SX(s,t+1) = 0.25*CHORDS(s,t)+geo.b(s,t)*(tan(geo.SW(s,t))) -0.25*CHORDS(s,t+1) + SX(s,t);
        
        % This was added by Dario to allow correct orientation of panels
        if geo.dihed(s,t) == pi/2
            SY(s,t+1) = geo.b(s,t)*cos(geo.dihed(s,t)) + SY(s,t);
            SZ(s,t+1) = geo.b(s,t)*sin(geo.dihed(s,t)) + SZ(s,t);
        else
            SY(s,t+1) = geo.b(s,t) + SY(s,t);
            SZ(s,t+1) = geo.b(s,t)*tan(geo.dihed(s,t)) + SZ(s,t);
        end
        
    end
    
end

%MAIN GEOMETRY SETUP LOOP, CREATES Partition QUAD PANELS, VORTICIES AND COLL-POINTS
flapc = 0;
dofc = 0;
maxt = 0;

for s=1:noofwings
    for t=1:noofloops(s) %setuploop
        
        %s
        %t
        
        %geo.ny(s,t)
        [C, V, N2, DN, camber, twist, P, ndof, cdof, hinge] = VortCollNorm(geo.fnx(s,t),geo.ny(s,t),geo.nx(s,t),geo.fsym(s,t),geo.fc(s,t,:),geo.flapped(s,t),geo.TW(s,t,:),geo.foil(s,t,:),...
            geo.T(s,t),geo.SW(s,t),CHORDS(s,t),geo.dihed(s,t),geo.b(s,t),geo.symetric(s),SX(s,t),SY(s,t),SZ(s,t),geo.meshtype(s,t));
        
        if (geo.b(s,t)<0)
            %fprintf(fid, '\n\t\t### Warning: patch %d - part %d will be reverted (negative span given).', s,t);
        end
        
        lattice.DOF(s,t,1) = 1 + dofc;
        lattice.DOF(s,t,2) = dofc + ndof;
        
        lattice.COLLOC = [lattice.COLLOC; C];
        lattice.VORTEX = [lattice.VORTEX; V];
        lattice.N      = [lattice.N; N2];
        lattice.DN     = [lattice.DN; DN];
        lattice.Camber = [lattice.Camber; camber];
        lattice.Twist  = [lattice.Twist; twist];
        %figure(1);hold on
        %cdiff = (CHORDS(s,t)-(C(end,1)-C(1,1)))/2;
        %plot3((C(:,1)-0.5*(SX(s,1)+SX(s,2)))/(0.5*(CHORDS(s,1)+CHORDS(s,2))),C(:,2),-camber,'ro-');
        
        if geo.flapped(s,t)
            
            flapc = flapc +1;
            
            lattice.Control.Patch(flapc) = s;
            lattice.Control.Part(flapc) = t;
            lattice.Control.Hinge(flapc,:,:) = hinge;
            if TOTAL(s,t) == 1
                lattice.Control.DOF(flapc).data =  dofc+1:dofc + ndof;
            else
                lattice.Control.DOF(flapc).data = cdof + dofc;
            end
        end
        
        dofc = dofc + ndof;
        
        Cmgc(s,t) = CHORDS(s,t) * ((1+geo.T(s,t)))/2;
        S(s,t) = geo.b(s,t) * Cmgc(s,t);  % patch area
        
        if geo.symetric(s)==1
            
            S(s,t) = S(s,t)*2;
            
        end
        
        lattice.XYZ=[lattice.XYZ;P];
        
        if noofloops(s) > maxt
            
            maxt = noofloops(s);
            
        end
    end
    
end

%if (ref.b_ref==0)

%   ref.b_ref = sum(sum(abs(geo.b) + abs(geo.b) .* repmat(geo.symetric', 1, maxt),1)) / geo.nwing;

%ref.b_ref=ref.b_ref*(geo.symetric(1)+1);

%end

%if (ref.S_ref == 0)

%   ref.S_ref = sum(sum(abs(S) + abs(S) .* repmat(geo.symetric', 1, maxt),1)) / geo.nwing;

%end

%if (ref.C_mgc == 0)

%  ref.C_mgc = sum(sum(abs(Cmgc) + abs(Cmgc) .* repmat(geo.symetric', 1, maxt),1)) / geo.nwing;

%end

%if (ref.C_mac == 0)

%   [ref.C_mac void] = fCmac(CHORDS(1,1:2),geo.b(1,:),geo.SW(1,:),SX(1,:),SY(1,:),SZ(1,:),geo.dihed(1,:),geo.symetric(1)); %Main (first) wing Mean aerodymaic chord calculation

%end

%if isempty(ref.mac_pos)
%   [void ref.mac_pos]=fCmac(CHORDS(1,1:2),geo.b(1,:),geo.SW(1,:),...
%       SX(1,:),SY(1,:),SZ(1,:),geo.dihed(1,:),geo.symetric(1)); %Main (first) wing Mean aerodymaic chord calculation
%mac_pos=-mac_pos
%end

end