time.limits = [0.00 5.00];
time.dt     = 0.01;

Aero.gust.H       = 60;
Aero.gust.xg      =  0;
Aero.gust.yg      =  0;
Aero.gust.zg      =  0;
Aero.gust.wgmax   = 00;
Aero.gust.theta   = 00*pi/180;
Aero.gust_switch  =  1;

x_range            = [  0: -1.0:-Aero.gust.H*2 ];
y_range            = [  1  ];
z_gust             = zeros(length(x_range),length(y_range));

%%
% load('C:\Users\rc15645\Google Drive\AEROGUST\Code\Beam Model\Lite Version v2\CranfieldWorkshopGustloads\NL\gust_H_60_theta_0.mat','x_out1')
load('C:\Users\rc15645\Google Drive\AEROGUST\Code\Beam Model\Lite Version v2\RollTest_0g.mat','x_out1')

%%
clear gust_tmp
for ii = 1:length(x_range)
	gust_tmp(:,ii) = getGust(x_range(ii),y_range,0,Aero.gust);
end

time_in   = time.limits(1):time.dt:time.limits(2);
[maxx,maxy,maxz] = deal(-1e99);
[minx,miny,minz] = deal( 1e99);
max_f            =      -1e99;
max_x_f          =      -1e99;

force_dirn = 3;

for jj = 1:length(Matrices.n_elem)
    maxx = max([max(max(x_out1.x_p{jj}(1:3:end,:))),maxx]);
    minx = min([min(min(x_out1.x_p{jj}(1:3:end,:))),minx]);
    maxy = max([max(max(x_out1.x_p{jj}(2:3:end,:))),maxy]);
    miny = min([min(min(x_out1.x_p{jj}(2:3:end,:))),miny]);
    maxz = max([max(max(x_out1.x_p{jj}(3:3:end,:))),maxz]);
    minz = min([min(min(x_out1.x_p{jj}(3:3:end,:))),minz]);
%     for ii = 1:3
%         max_f = max([max(max(abs(x_out1.aeroForce{jj}(ii:6:end,:))))],max_f);
%     end
    max_x_f = max([max(max(abs(x_out1.x_f{jj}(force_dirn:6:end,:)))),max_x_f]);
end

patchcolor0 = [1 1 1];
patchcolor1 = [1 0 0];

maxx  = maxx+1;
minx  = minx-1;
maxy  = maxy+1;
miny  = miny-1;
maxz  = maxz+1;
minz  = minz-1;
% max_f = max_f/3;

close all

figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
cg        = zeros(3,length(time_in));
writerObj = VideoWriter('tmp.avi');
open(writerObj)

%% SciTech
x_surr = 30;
y_surr = 30;
z_surr = 30;%15;
%% AeroGust UAV
% x_surr = 7;
% y_surr = 7;
% z_surr = 7;
%%
for ii = 1:length(time_in)
%     n_elem_tot = 0;
%     surf(x_range,y_range,z_gust','LineStyle','none','FaceAlpha',0.5)
    CGa = Quat2Rot(x_out1.x_qa(:,ii));
    if isfield(Matrices,'RenderPoints_RB')
        if isempty(Matrices.RenderPoints_RB)
            scatter3(x_out1.x_pa(1,ii),x_out1.x_pa(2,ii),x_out1.x_pa(3,ii))
        else
            for jj = 1:size(Matrices.RenderPoints_RB.x,1)
                coords_rot = CGa*[Matrices.RenderPoints_RB.x(jj,:);Matrices.RenderPoints_RB.y(jj,:);Matrices.RenderPoints_RB.z(jj,:)] + repmat(x_out1.x_pa(:,ii),1,size(Matrices.RenderPoints_RB.x,2));
                X_out(jj,:) = coords_rot(1,:);
                Y_out(jj,:) = coords_rot(2,:);
                Z_out(jj,:) = coords_rot(3,:);
            end
            surf(X_out,Y_out,Z_out,'FaceColor','w')
        end
        hold on
    end
%       for jj = 1:size(Matrices.RenderPoints_RB.x,2)
%         plot3(X_out(:,jj),Y_out(:,jj),Z_out(:,jj),'k')
%     end
    axis equal
    for jj = 1:length(Matrices.n_elem)
        x_f_tmp{jj} = x_out1.x_f{jj}(:,ii);
        plot3(x_out1.x_p{jj}(end-2,1:ii),x_out1.x_p{jj}(end-1,1:ii),x_out1.x_p{jj}(end,1:ii),'k--')
        for kk = 1:Matrices.n_elem(jj)
            ind1 = [1:6] + (kk-1)*6;
            ind3 = [1:3] + (kk-1)*3;
            ind4 = [1:4] + (kk-1)*4;
            CGB = Quat2Rot(x_out1.x_q{jj}(ind4,ii));
            coords_upper = CGB*[Matrices.RenderPoints{jj}.x_upper(kk,:);
                                Matrices.RenderPoints{jj}.y_upper(kk,:);
                                Matrices.RenderPoints{jj}.z_upper(kk,:)] + repmat(x_out1.x_p{jj}(ind3,ii),1,size(Matrices.RenderPoints{jj}.x_upper,2));
            coords_lower = CGB*[Matrices.RenderPoints{jj}.x_lower(kk,:);
                                Matrices.RenderPoints{jj}.y_lower(kk,:);
                                Matrices.RenderPoints{jj}.z_lower(kk,:)] + repmat(x_out1.x_p{jj}(ind3,ii),1,size(Matrices.RenderPoints{jj}.x_upper,2));
            elem_colour  = patchcolor0*(1 - abs(x_out1.x_f{jj}(ind1(force_dirn),ii))/max_x_f) + patchcolor1*abs(x_out1.x_f{jj}(ind1(force_dirn),ii))/max_x_f;
            patch('XData',coords_upper(1,:),'YData',coords_upper(2,:),'ZData',coords_upper(3,:),'FaceColor',elem_colour)%,'white')
            hold on
            patch('XData',coords_lower(1,:),'YData',coords_lower(2,:),'ZData',coords_lower(3,:),'FaceColor',elem_colour)%,'white')
%             aero = CGB*x_out1.aeroForce{jj}(ind1(1:3),ii)/max_f;
%             plot3([x_out1.x_p{jj}(ind3(1),ii) x_out1.x_p{jj}(ind3(1),ii)+aero(1)],...
%                   [x_out1.x_p{jj}(ind3(2),ii) x_out1.x_p{jj}(ind3(2),ii)+aero(2)],...
%                   [x_out1.x_p{jj}(ind3(3),ii) x_out1.x_p{jj}(ind3(3),ii)+aero(3)],'r')
        end
    end
%     aero = CGa*x_out1.total_Aero(1:3,ii)/max_f/10;
%     plot3([x_out1.x_pa(1,ii) x_out1.x_pa(1,ii)+aero(1)],...
%           [x_out1.x_pa(2,ii) x_out1.x_pa(2,ii)+aero(2)],...
%           [x_out1.x_pa(3,ii) x_out1.x_pa(3,ii)+aero(3)],'r')
    MassProps = calculateCoM(x_f_tmp,Matrices,0);
    cg_vec(:,ii) = CGa*MassProps.cg + x_out1.x_pa(:,ii);
    scatter3(cg_vec(1,ii),cg_vec(2,ii),cg_vec(3,ii))
    title(['Time = ',num2str(time_in(ii)),'s'],'FontSize',14)
%     set(gca,'XLim',[minx maxx],'YLim',[miny maxy],'ZLim',[minz maxz])
%     plot3(x_range,ones(size(x_range))*(cg_vec(2,ii)+25),gust_tmp(3,:),'k','LineWidth',2)
%     plot3(x_range,gust_tmp(2,:),ones(size(x_range))*(cg_vec(3,ii)-15),'k','LineWidth',2)
    fill3(x_range,ones(size(x_range))*(cg_vec(2,ii)+30),gust_tmp(3,:),[0.88 1 1]);
    fill3(x_range,gust_tmp(2,:),ones(size(x_range))*(cg_vec(3,ii)-15),[0.88 1 1])
    set(gca,'XLim',[cg_vec(1,ii)-x_surr cg_vec(1,ii)+x_surr],'YLim',[cg_vec(2,ii)-y_surr cg_vec(2,ii)+y_surr],'ZLim',[cg_vec(3,ii)-z_surr cg_vec(3,ii)+z_surr])
    set(gca,'LineWidth',2)
%     set(h1,'FaceAlpha',0.1)
%     set(gca,'YLim',[-25 25],'ZLim',[-10 10])
    grid on
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    set(gca,'FontSize',14)
    hold off 
%     view([ 0 0]) %Side View
%     view([-90 0]) % Front View
    pause(0.001)
%     set(gcf,'Renderer','zbuffer');
    set(gcf,'Renderer','OpenGL');
    ax = gcf;
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
%     ti = ax.TightInset;
%     rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
%     marg = 24;
%     rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    rect = [0, 0, pos(3), pos(4)];
    F(ii) = getframe(gcf);%,rect
    writeVideo(writerObj,F(ii));
end

close(writerObj)