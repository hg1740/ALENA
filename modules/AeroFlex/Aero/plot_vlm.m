%% plot_vlm
% The following script takes the lattice structure array as an input and
% plots the collocation points, panel vertices and vortex midpoints

function plot_vlm(lattice,XYZ_PLOT,NORMAL_PLOT,COLLOC_PLOT,VORTEX_PLOT,WAKE_PLOT,CAMBER_PLOT,TWIST_PLOT)

switch nargin
    case 1
        COLLOC_PLOT = 0;
        VORTEX_PLOT = 0;
        WAKE_PLOT = 0;
        NORMAL_PLOT = 0;
        XYZ_PLOT = 0;
        CAMBER_PLOT = 0;
        TWIST_PLOT = 0;
    case 2
        VORTEX_PLOT = 0;
        WAKE_PLOT = 0;
        COLLOC_PLOT = 0;
        CAMBER_PLOT = 0;
        NORMAL_PLOT = 0;
        TWIST_PLOT = 0;
    case 3
        CAMBER_PLOT = 0;
        WAKE_PLOT = 0;
        VORTEX_PLOT = 0;
        COLLOC_PLOT = 0;
        TWIST_PLOT = 0;
    case 4
        WAKE_PLOT = 0;
        CAMBER_PLOT = 0;
        VORTEX_PLOT = 0;
        TWIST_PLOT = 0;
    case 5
        WAKE_PLOT= 0;
        CAMBER_PLOT = 0;
        TWIST_PLOT=0;
    case 6
        CAMBER_PLOT = 0;
        TWIST_PLOT = 0;
    case 7
        TWIST_PLOT = 0;
    case 8
        
    otherwise
        COLLOC_PLOT = 0;
        VORTEX_PLOT = 0;
        WAKE_PLOT = 0;
        NORMAL_PLOT = 0;
        CAMBER_PLOT = 0;
        XYZ_PLOT = 1;
        CAMBER_PLOT = 0;
        TWIST_PLOT = 0;
        
end



[npwing, ~, ~] = size(lattice.XYZ);


figure(444);

if XYZ_PLOT
    plot3(lattice.XYZ(1:npwing,:,1)',lattice.XYZ(1:npwing,:,2)',lattice.XYZ(1:npwing,:,3)','-bo',...
        'MarkerSize', 1, 'MarkerFaceColor','b','LineWidth',2);
    hold on;
end

if COLLOC_PLOT
plot3(lattice.COLLOC(1:npwing,1),lattice.COLLOC(1:npwing,2),lattice.COLLOC(1:npwing,3),'bo',...
    'MarkerSize', 5, 'MarkerFaceColor','b');
hold on;
end

if VORTEX_PLOT 
    [npvort, nv, ~] = size(lattice.VORTEX);
switch nv
    case 5
        c = [1 2];
    case 6
        c = [3 4];
    case 8
        c = [4 5];  
    otherwise
        error('Wrong VORTEX dimension.');  
end

x1 = zeros(npwing,2);
x2 = zeros(npwing,2);
x3 = zeros(npwing,2);
v1 = zeros(npwing,1);
v2 = zeros(npwing,1);
v3 = zeros(npwing,1);

x1 = lattice.VORTEX(1:npwing,c,1);
x2 = lattice.VORTEX(1:npwing,c,2);
x3 = lattice.VORTEX(1:npwing,c,3);

%This takes the midpoint between the two spanwise edges of the panel
v1 = mean(x1,2);
v2 = mean(x2,2);
v3 = mean(x3,2);

% This plots the Vortex Location POINTS
plot3(v1, v2, v3,'ro',  'MarkerSize', 5, 'MarkerFaceColor','r');
hold on;
%legend('XYZ','COLLOC','Mean Vortex Segments','Location','Best');

end

if WAKE_PLOT
    for s = 1:npvort 
        w=0;
        for u = 1:nv-1
            w = w + 1;
            VX(w,:) = [lattice.VORTEX(s, u, 1) lattice.VORTEX(s, u+1, 1)];
            VY(w,:) = [lattice.VORTEX(s, u, 2) lattice.VORTEX(s, u+1, 2)];
            VZ(w,:) = [lattice.VORTEX(s, u, 3) lattice.VORTEX(s, u+1, 3)];
        end
        
        plot3(VX, VY, VZ,'k--','LineWidth',1);
        hold on; 
    end 
end

if NORMAL_PLOT
    
    quiver3(lattice.COLLOC(:,1),lattice.COLLOC(:,2),lattice.COLLOC(:,3),...
        -lattice.N(:,1),-lattice.N(:,2),-lattice.N(:,3),0.5,'Color','r');hold on;
    axis equal;
   
end

if CAMBER_PLOT && ~TWIST_PLOT
    hold on;
    colormap(jet);
    numpanels=size(lattice.XYZ,1);
    s=fill3(lattice.XYZ(1:numpanels,:,1)', lattice.XYZ(1:numpanels,:,2)', ...
        lattice.XYZ(1:numpanels,:,3)', -lattice.Camber','FaceColor','flat','EdgeColor','none');hold on;
end

if TWIST_PLOT && ~CAMBER_PLOT
    hold on;
    colormap(jet);
    numpanels=size(lattice.XYZ,1);
    s=fill3(lattice.XYZ(1:numpanels,:,1)', lattice.XYZ(1:numpanels,:,2)', ...
        lattice.XYZ(1:numpanels,:,3)', lattice.Twist','FaceColor','flat','EdgeColor','none');hold on;
end

if CAMBER_PLOT && TWIST_PLOT
        hold on;
    colormap(jet);
    numpanels=size(lattice.XYZ,1);
    s=fill3(lattice.XYZ(1:numpanels,:,1)', lattice.XYZ(1:numpanels,:,2)', ...
        lattice.XYZ(1:numpanels,:,3)', lattice.Twist' - lattice.Camber','FaceColor','flat','EdgeColor','none');hold on;
end
% 
% if isfield(lattice,'wingvel')
%             hold on;
%     colormap(jet);
%     numpanels=size(lattice.XYZ,1);
%     wingvel(:,1) = beam_model.Aero.Interp.Ic*lattice.wingvel(:,1);
%     wingvel(:,2) = beam_model.Aero.Interp.Ic*lattice.wingvel(:,2);
%     wingvel(:,3) = beam_model.Aero.Interp.Ic*lattice.wingvel(:,3);
%     
%     s=fill3(lattice.XYZ(1:numpanels,:,1)', lattice.XYZ(1:numpanels,:,2)', ...
%         lattice.XYZ(1:numpanels,:,3)', wingvel(:,2)','FaceColor','flat','EdgeColor','none');hold on;
% end



axis tight
axis equal
view([201,34]);

end