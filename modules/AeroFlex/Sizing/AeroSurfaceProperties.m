function [Param] = AeroSurfaceProperties(geo,ref)

%% Get Profile coordinates
Nsec = sum(geo.ny);
x_3 = [];
y_3 = [];
Zu_3 = [];
Zl_3 = [];
for i = 1:numel(geo.ny)
    
    [EP, Zu, Zl, ~, ~] = airfoil_mean_line(geo.foil(i,1,1),0);
    x_prof = (EP * geo.c(i) + geo.startx(i))';
    y_prof = (geo.starty(i)*ones(size(EP)))';
    z_prof = (geo.startz(i)*ones(size(EP)))';
    Zu = Zu * geo.c(i) + z_prof;
    Zl = Zl * geo.c(i) + z_prof;
    x_3 = [x_3;x_prof];
    y_3 = [y_3;y_prof];
    Zu_3 = [Zu_3;Zu];
    Zl_3 = [Zl_3;Zl];
end

% Plot the tip foil profile
[EP, Zu, Zl, ~, ~] = airfoil_mean_line(geo.foil(i,1,2),0);
x_prof = (EP * geo.c(i)*geo.T(i) + geo.startx(i) + 0.25*geo.c(i) + geo.b(i).*tan(geo.SW(i)) - 0.25*geo.T(i).*geo.c(i))'; 
y_prof = ((geo.starty(i)+geo.b(i))*ones(size(EP)))';
z_prof = ((geo.startz(i) + geo.b(i)*tan(geo.dihed(i)))*ones(size(EP)))';
Zu = Zu * geo.c(i)*geo.T(i) + z_prof;
Zl = Zl * geo.c(i)*geo.T(i) + z_prof;
x_3 = [x_3;x_prof];
y_3 = [y_3;y_prof];
Zu_3 = [Zu_3;Zu];
Zl_3 = [Zl_3;Zl];

%% Calculate Wetted surface area
x_diff  = diff(x_3,1,2);
y_diff  = diff(y_3,1,1);
Zu_diff = diff(Zu_3,1,2); 
Zl_diff = diff(Zl_3,1,2); 

% FIND THE X VALUES THAT CORRESPOND TO THE SPAR LOCATIONS
% ISOLATE THE UPPER AND LOWER CONTOURS AND ADD THOSE UP SEPARATELY

distu = sqrt(x_diff.^2+Zu_diff.^2);
distl = sqrt(x_diff.^2+Zl_diff.^2);

panelUpContLength = (distu(1:Nsec,:) + distu(2:Nsec+1,:));
panelLowerContLength = (distl(1:Nsec,:) + distl(2:Nsec+1,:));

numSurf = size(y_diff,2)-1;
panelSpanLength = (y_diff(:,1:numSurf)+y_diff(:,2:numSurf+1))/2;

panelareau = panelUpContLength.*panelSpanLength/2;
panelareal = panelLowerContLength.*panelSpanLength/2;

v1 = [geo.startx,geo.starty,geo.startz];
v2 = v1 + [0.25*geo.c + geo.b.*tan(geo.SW) - 0.25*geo.T.*geo.c, geo.b, geo.b.*tan(geo.dihed)];
v3 = v1 + [0.25*geo.c + geo.b.*tan(geo.SW) + 0.75*geo.T.*geo.c,geo.b,geo.b.*tan(geo.dihed)];
v4 = v1 + [geo.c,0.*geo.b,0.*geo.b];

Swet_c = sum(panelareau+panelareal,2);
Swet = sum(sum(panelareau+panelareal));

for i = 1:length(geo.startx)
    Sref_c(i,1) = polyarea([v1(i,1);v2(i,1);v3(i,1);v4(i,1)],[v1(i,2);v2(i,2);v3(i,2);v4(i,2)]); % neglects dihedral
end

Swet_Sref = Swet/ref.S_ref;

%% Find maximum thickness and location
[maxt,t_index] = max(Zu_3-Zl_3,[],2);
t_c = maxt./[geo.c;geo.c(end)*geo.T(end)];
maxtloc = [];
for i = 1:size(x_3,1)
    maxtloc = [maxtloc;x_3(i,t_index(i))];
end

x_c = (maxtloc - x_3(:,1))./[geo.c;geo.c(end)*geo.T(end)];
mean_x_c = 0.5*(x_c(1:end-1)+x_c(2:end));
mean_t_c = 0.5*(t_c(1:end-1)+t_c(2:end));

Param.Swet_c = Swet_c;
Param.Swet = Swet;
Param.Sref_c = Sref_c;
Param.Swet_Sref = Swet_Sref;
Param.x_c = mean_x_c;
Param.t_c = mean_t_c;
end
