%% PLOTNEOCASSMODELFUNC Plot the strcutural elements of a NEOCASS model

%   Copyright 2016 University of Bristol
%   Private function.
function plotNeocassModelFunc(obj)

% Wing

h = figure(10);
if obj.Opts.Geom.addTail
    Parts = {'StbdWing','StbdHTP','VTP','PortWing','PortHTP'};
else
    Parts = {'Wing'};
end

for i = 1:length(Parts)
    Aeronodes = findobj(obj.Mdl.Grid,'type','Aerodynamic','-and','part',Parts{i});
    
    numnodes = numel(Aeronodes);
    
    Xnle   = struct2mat(Aeronodes(1:numnodes/2),'coord',1);
    Ynle   = struct2mat(Aeronodes(1:numnodes/2),'coord',2);
    Znle   = struct2mat(Aeronodes(1:numnodes/2),'coord',3);
    
    Xnte   = struct2mat(Aeronodes(numnodes/2+1:numnodes),'coord',1);
    Ynte   = struct2mat(Aeronodes(numnodes/2+1:numnodes),'coord',2);
    Znte   = struct2mat(Aeronodes(numnodes/2+1:numnodes),'coord',3);
    
    Xa = [Xnle;flipud(Xnte);Xnle(1)];
    Ya = [Ynle;flipud(Ynte);Ynle(1)];
    Za = [Znle;flipud(Znte);Znle(1)];
    
    plot3(Xa,Ya,Za,'m-','LineWidth',1);
    hold on;
    if obj.Opts.Geom.reflect == 1
        plot3(Xa,-Ya,Za,'m-','LineWidth',1);
        axis equal;
    end
    
end
% Plot the Caero Outlines
[h_w,h_c] = obj.Mdl.Caero.plot(10);

% Plot the beam connections
h_be = obj.Mdl.plotbar(10);

% Plot the beam nodes
h_b = obj.Mdl.Grid.plotbeam(10);

% Plot the aeronodes
h_a = obj.Mdl.Grid.plotaero(10);

% Plot the Aero connections
h_r0 = obj.Mdl.plotrbe0(10);

% Plot the Aero connections
h_r2 = obj.Mdl.plotrbe2(10);

axis equal;

figure(10);
hold on;
Conm2 = obj.Mdl.Conm2;
% Check for the size of the mass

Masses = [obj.Mdl.Conm2.m]';
minmass = min(Masses); % Size 1
maxmass = max(Masses); % Size 5

Sizemax = 6;
Sizemin = 6;
h_conm    = hggroup;
h_conmoff = hggroup;

for i = 1:numel(Conm2)
    
    Conmidx = find([obj.Mdl.Grid.id]' == Conm2(i).grid);
    coords  = obj.Mdl.Grid(Conmidx).coord + Conm2(i).offset;
    NodeC = obj.Mdl.Grid(Conmidx).coord;
    
    a = (Sizemax-Sizemin)/(maxmass-minmass);
    b = Sizemin - a*minmass;
    MS = a*obj.Mdl.Conm2(i).m + b;
    
    plot3(coords(1),coords(2),coords(3),'MarkerFaceColor','g','Marker','o','Markersize',MS,'MarkerEdgeColor','k','Parent',h_conm);
    plot3([coords(1) NodeC(1)],[coords(2) NodeC(2)],[coords(3) NodeC(3)],'g-','Parent',h_conmoff);
end

set(gca,'Xgrid','on','Ygrid','on','Zgrid','on')
axis equal

legend([h_w,h_c,h_be,h_b,h_a,h_r0,h_r2,h_conm,h_conmoff],'Wing','Controls','Beam Elements','Beam nodes','Aero nodes','RBE0','RBE2','Conm2','Conm2 Offset','Location','Best');

saveas(h,[obj.Dirs.writeFolder filesep 'ModelFigures' filesep 'AircraftModel.fig']);


end
