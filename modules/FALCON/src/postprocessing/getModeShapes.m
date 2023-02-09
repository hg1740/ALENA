function getModeShapes(SYSA,Matrices,x_dyn,modes,Sim,render_flag,ghost_flag,vid_flag,saveaspath)

Sim.Ind = getSystemInd(Matrices,Sim,[]);

if ~exist('vid_flag','var')
    vid_flag = 0;
end
if ~exist('render_flag','var')
    render_flag = 0;
end
if ~exist('ghost_flag','var')
    ghost_flag = 0;
end
if ~exist('saveaspath','var')
    saveaspath = pwd;
elseif strcmp(saveaspath(end),filesep)
    saveaspath = saveaspath(1:end-1);
end

[vecs,vals] = eig(full(SYSA));
[~,ind] = sort(abs(imag(diag(vals))));
vals    = diag(vals);
vals    = vals(ind);
vecs    = vecs(:,ind);
counter = 0;
del_ind = [];
for ii = 1:length(vals)
    if isreal(vals(ii))
        counter = counter + 1;
        del_ind = [del_ind;counter];
    elseif abs(real(vals(ii)))<1e-3&&abs(imag(vals(ii)))<1e-2
        counter = counter + 1;
        del_ind = [del_ind;counter];
    end
end
vals(  del_ind) = [];
vecs(:,del_ind) = [];
vals    = vals(  1:2:end);
freq    = abs(vals);
damp    = real(vals)./freq;
vecs    = vecs(:,1:2:end);

mult_fact = 2;

if vid_flag <= 0
    %% Still Images
    figure
    scatter(real(vals(modes)),imag(vals(modes)),'k.'),hold on
    grid on, grid minor
    title('Root Locus Plot')
    xlabel('Imaginary Part')
    ylabel('Real Part')
    set(gca,'FontSize',24)
    if vid_flag == -1
        saveas(gcf,[saveaspath,filesep,'Root Loci.fig'])
        close(gcf)
    end
    figure
    scatter(real(vals(modes))./abs(vals(modes)),abs(vals(modes))/2/pi,'k.'),hold on
    grid on, grid minor
    title('Root Locus Plot')
    xlabel('Damping')
    ylabel('Frequency (Hz)')
    set(gca,'FontSize',24)
    if vid_flag == -1
        saveas(gcf,[saveaspath,filesep,'Root Loci FrqDamp.fig'])
        close(gcf)
    end
    vectors   = real(vecs(:,modes));
    if render_flag == 0
        %% Just draw the sticks
        for ii = 1:length(modes)
            figure
            max_disp  = 0;
            for jj = 1:length(Matrices.n_elem)
                max_disp = max(max_disp,max(max(abs(vectors(Sim.Ind.x_p_ind{jj},ii)))));
            end
            for jj = 1:length(Matrices.n_elem)
                vectors_p{jj} = vectors(Sim.Ind.x_p_ind{jj},ii)*mult_fact/max_disp+x_dyn.x_p{jj};
                plot3(vectors_p{jj}(1:3:end),vectors_p{jj}(2:3:end),vectors_p{jj}(3:3:end),'k');hold on
                if ghost_flag
                    plot3(x_dyn.x_p{jj}(1:3:end),   x_dyn.x_p{jj}(2:3:end),   x_dyn.x_p{jj}(3:3:end)   ,'r')
                end
            end
            axis equal
            title(['Mode ',num2str(modes(ii)),', Frequency = ',num2str(freq(modes(ii))/2/pi),'Hz (Damping = ',num2str(damp(modes(ii))),')'])
            hold off
            pause(0.001)
            if vid_flag == -1
                saveas(gcf,[saveaspath,filesep,'Mode ',num2str(modes(ii)),'.fig'])
            end
        end
    else
        %% Render the whole aircraft
        for ii = 1:length(modes)
            figure
            max_disp  = 0;
            for jj = 1:length(Matrices.n_elem)
                max_disp = max(max_disp,max(max(abs(vectors(Sim.Ind.x_p_ind{jj},ii)))));
            end
            vectors_pa = vectors(Sim.Ind.x_pa_ind,ii)*mult_fact/max_disp+x_dyn.x_pa;
            vectors_qa = vectors(Sim.Ind.x_qa_ind,ii)*mult_fact/max_disp+x_dyn.x_qa;
            CGa = Quat2Rot(vectors_qa);
            if isfield(Matrices,'RenderPoints_RB')
                if isempty(Matrices.RenderPoints_RB)
                    scatter3(vectors_pa(1),vectors_pa(2),vectors_pa(3))
                else
                    for jj = 1:size(Matrices.RenderPoints_RB.x,1)
                        coords_rot = CGa*[Matrices.RenderPoints_RB.x(jj,:);Matrices.RenderPoints_RB.y(jj,:);Matrices.RenderPoints_RB.z(jj,:)] + repmat(vectors_pa,1,size(Matrices.RenderPoints_RB.x,2));
                        X_out(jj,:) = coords_rot(1,:);
                        Y_out(jj,:) = coords_rot(2,:);
                        Z_out(jj,:) = coords_rot(3,:);
                    end
                    surf(X_out,Y_out,Z_out,'FaceColor','w')
                end
                hold on
            end
            clear X_out Y_out Z_out
            if ghost_flag
                CGa = Quat2Rot(x_dyn.x_qa);
                if isfield(Matrices,'RenderPoints_RB')
                    if isempty(Matrices.RenderPoints_RB)
                        scatter3(x_dyn.x_pa(1),x_dyn.x_pa(2),x_dyn.x_pa(3))
                    else
                        for jj = 1:size(Matrices.RenderPoints_RB.x,1)
                            coords_rot = CGa*[Matrices.RenderPoints_RB.x(jj,:);Matrices.RenderPoints_RB.y(jj,:);Matrices.RenderPoints_RB.z(jj,:)] + repmat(x_dyn.x_pa,1,size(Matrices.RenderPoints_RB.x,2));
                            X_out(jj,:) = coords_rot(1,:);
                            Y_out(jj,:) = coords_rot(2,:);
                            Z_out(jj,:) = coords_rot(3,:);
                        end
                        surf(X_out,Y_out,Z_out,'FaceColor','w','FaceAlpha',0.25,'EdgeAlpha',0.25)
                    end
                    hold on
                end
                clear X_out Y_out Z_out
            end
            for jj = 1:length(Matrices.n_elem)
                if isfield(Matrices.RenderPoints_RB,'ExtraInd')
                    extra_parts = find(Matrices.RenderPoints_RB.ExtraInd(:,1)==jj);
                else
                    extra_parts = [];
                end
                vectors_p{jj} = vectors(Sim.Ind.x_p_ind{jj},ii)*mult_fact/max_disp+x_dyn.x_p{jj};
                vectors_q{jj} = vectors(Sim.Ind.x_q_ind{jj},ii)*mult_fact/max_disp+x_dyn.x_q{jj};
                for kk = 1:Matrices.n_elem(jj)
                    if isfield(Matrices.RenderPoints_RB,'ExtraInd')
                        extra_parts_ind = find(Matrices.RenderPoints_RB.ExtraInd(:,2)==kk);
                    else
                        extra_parts_ind = [];
                    end
                    ind1 = [1:6] + (kk-1)*6;
                    ind3 = [1:3] + (kk-1)*3;
                    ind4 = [1:4] + (kk-1)*4;
                    CGB = Quat2Rot(vectors_q{jj}(ind4));
                    coords_upper = CGB*[Matrices.RenderPoints{jj}.x_upper(kk,:);
                                        Matrices.RenderPoints{jj}.y_upper(kk,:);
                                        Matrices.RenderPoints{jj}.z_upper(kk,:)] + repmat(vectors_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                    coords_lower = CGB*[Matrices.RenderPoints{jj}.x_lower(kk,:);
                                        Matrices.RenderPoints{jj}.y_lower(kk,:);
                                        Matrices.RenderPoints{jj}.z_lower(kk,:)] + repmat(vectors_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                    patch('XData',coords_upper(1,:),'YData',coords_upper(2,:),'ZData',coords_upper(3,:),'FaceColor','white')
                    hold on
                    patch('XData',coords_lower(1,:),'YData',coords_lower(2,:),'ZData',coords_lower(3,:),'FaceColor','white')
                    if ~isempty(extra_parts)&&~isempty(extra_parts_ind)
                        for ee = 1:size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,1)
                            coords_rot = CGB*[Matrices.RenderPoints_RB.Extra{extra_parts}.x(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.y(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.z(ee,:)] + repmat(vectors_p{extra_parts}(ind3),1,size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,2));
                            X_out(ee,:) = coords_rot(1,:);
                            Y_out(ee,:) = coords_rot(2,:);
                            Z_out(ee,:) = coords_rot(3,:);
                        end
                        surf(X_out,Y_out,Z_out,'FaceColor','w')
                    end
                    if ghost_flag
                        CGB = Quat2Rot(x_dyn.x_q{jj}(ind4));
                        coords_upper = CGB*[Matrices.RenderPoints{jj}.x_upper(kk,:);
                                            Matrices.RenderPoints{jj}.y_upper(kk,:);
                                            Matrices.RenderPoints{jj}.z_upper(kk,:)] + repmat(x_dyn.x_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                        coords_lower = CGB*[Matrices.RenderPoints{jj}.x_lower(kk,:);
                                            Matrices.RenderPoints{jj}.y_lower(kk,:);
                                            Matrices.RenderPoints{jj}.z_lower(kk,:)] + repmat(x_dyn.x_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                        patch('XData',coords_upper(1,:),'YData',coords_upper(2,:),'ZData',coords_upper(3,:),'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        hold on
                        patch('XData',coords_lower(1,:),'YData',coords_lower(2,:),'ZData',coords_lower(3,:),'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        if ~isempty(extra_parts)&&~isempty(extra_parts_ind)
                            for ee = 1:size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,1)
                                coords_rot = CGB*[Matrices.RenderPoints_RB.Extra{extra_parts}.x(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.y(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.z(ee,:)] + repmat(x_dyn.x_p{extra_parts}(ind3),1,size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,2));
                                X_out(ee,:) = coords_rot(1,:);
                                Y_out(ee,:) = coords_rot(2,:);
                                Z_out(ee,:) = coords_rot(3,:);
                            end
                            surf(X_out,Y_out,Z_out,'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        end
                    end
                end
                
            end
            axis equal
            title(['Mode ',num2str(modes(ii)),', Frequency = ',num2str(freq(modes(ii))/2/pi),'Hz (Damping = ',num2str(damp(modes(ii))),')'])
            xlabel('X (m)')
            ylabel('Y (m)')
            zlabel('Z (m)')
            set(gca,'FontSize',24)
            hold off
            pause(0.001)
            if vid_flag == -1
                saveas(gcf,[saveaspath,filesep,'Mode ',num2str(modes(ii)),'.fig'])
                close(gcf)
            end
        end
    end
else
    %% Video
    if length(modes)>1
        fprintf('\nRequested number of modes greater than one! Video will only be generated for the first mode of the array!\n')
    end
    Mode2Plot = modes(1);
    figure
    if vid_flag == 2
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        writerObj = VideoWriter([saveaspath,filesep,'Mode',num2str(Mode2Plot),'.avi']);
        open(writerObj)
    end
    if render_flag == 1
        cycle_vec = 0:pi/10:4*pi;
    else
        cycle_vec = 0:pi/10:4*pi;
    end
    vectors   = real(vecs(:,Mode2Plot)*exp(1i*cycle_vec));
    if render_flag == 0
        %% Just draw the sticks
        max_disp  = 0;
        for jj = 1:length(Matrices.n_elem)
            max_disp = max(max_disp,max(max(abs(vectors(Sim.Ind.x_p_ind{jj},:)))));
        end                
        for jj = 1:length(Matrices.n_elem)
            vectors_p{jj} = vectors(Sim.Ind.x_p_ind{jj},:)*mult_fact/max_disp+repmat(x_dyn.x_p{jj},1,length(cycle_vec));
        end
        for ii = 1:length(cycle_vec)
            for jj = 1:length(Matrices.n_elem)
                plot3(vectors_p{jj}(1:3:end,ii),vectors_p{jj}(2:3:end,ii),vectors_p{jj}(3:3:end,ii),'k')
                hold on
                plot3(x_dyn.x_p{jj}(1:3:end),   x_dyn.x_p{jj}(2:3:end),   x_dyn.x_p{jj}(3:3:end)   ,'r')
            end
            axis equal
            if ii == 1
                axes_props = get(gcf,'children');
                xlim_init  = axes_props.XLim;
                ylim_init  = axes_props.YLim;
                zlim_init  = axes_props.ZLim;
            end
            set(gca,'XLim',xlim_init+[-1 1],'YLim',ylim_init+[-1 1],'ZLim',zlim_init+[-2 2])
            title(['Mode ',num2str(Mode2Plot),', Frequency = ',num2str(freq(Mode2Plot)/2/pi),'Hz (Damping = ',num2str(damp(Mode2Plot)),')'])
            hold off
            pause(0.001)
            if vid_flag == 2
                F(ii) = getframe(gcf);%,rect
                writeVideo(writerObj,F(ii));
            end
        end
    else
        %% Render the whole aircraft
        max_disp  = 0;
        for jj = 1:length(Matrices.n_elem)
            max_disp = max(max_disp,max(max(abs(vectors(Sim.Ind.x_p_ind{jj},:)))));
        end
        vectors_pa = vectors(Sim.Ind.x_pa_ind,:)*mult_fact/max_disp+repmat(x_dyn.x_pa,1,length(cycle_vec));
        vectors_qa = vectors(Sim.Ind.x_qa_ind,:)*mult_fact/max_disp+repmat(x_dyn.x_qa,1,length(cycle_vec));
        for jj = 1:length(Matrices.n_elem)
            vectors_p{jj} = vectors(Sim.Ind.x_p_ind{jj},:)*mult_fact/max_disp+repmat(x_dyn.x_p{jj},1,length(cycle_vec));
            vectors_q{jj} = vectors(Sim.Ind.x_q_ind{jj},:)*mult_fact/max_disp+repmat(x_dyn.x_q{jj},1,length(cycle_vec));
        end
        for ii = 1:length(cycle_vec)
            CGa = Quat2Rot(vectors_qa(:,ii));
            if isfield(Matrices,'RenderPoints_RB')
                if isempty(Matrices.RenderPoints_RB)
                    scatter3(vectors_pa(1,ii),vectors_pa(2,ii),vectors_pa(3,ii))
                else
                    for jj = 1:size(Matrices.RenderPoints_RB.x,1)
                        coords_rot = CGa*[Matrices.RenderPoints_RB.x(jj,:);Matrices.RenderPoints_RB.y(jj,:);Matrices.RenderPoints_RB.z(jj,:)] + repmat(vectors_pa(:,ii),1,size(Matrices.RenderPoints_RB.x,2));
                        X_out(jj,:) = coords_rot(1,:);
                        Y_out(jj,:) = coords_rot(2,:);
                        Z_out(jj,:) = coords_rot(3,:);
                    end
                    surf(X_out,Y_out,Z_out,'FaceColor','w')
                end
                hold on
                clear X_out Y_out Z_out
                if ghost_flag
                    CGa = Quat2Rot(x_dyn.x_qa);
                    if isfield(Matrices,'RenderPoints_RB')
                        if isempty(Matrices.RenderPoints_RB)
                            scatter3(x_dyn.x_pa(1),x_dyn.x_pa(2),x_dyn.x_pa(3))
                        else
                            for jj = 1:size(Matrices.RenderPoints_RB.x,1)
                                coords_rot = CGa*[Matrices.RenderPoints_RB.x(jj,:);Matrices.RenderPoints_RB.y(jj,:);Matrices.RenderPoints_RB.z(jj,:)] + repmat(x_dyn.x_pa,1,size(Matrices.RenderPoints_RB.x,2));
                                X_out(jj,:) = coords_rot(1,:);
                                Y_out(jj,:) = coords_rot(2,:);
                                Z_out(jj,:) = coords_rot(3,:);
                            end
                            surf(X_out,Y_out,Z_out,'FaceColor','w','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        end
                        hold on
                    end
                    clear X_out Y_out Z_out
                end
            end
            for jj = 1:length(Matrices.n_elem)
                if isfield(Matrices.RenderPoints_RB,'ExtraInd')
                    extra_parts = find(Matrices.RenderPoints_RB.ExtraInd(:,1)==jj);
                else
                    extra_parts = [];
                end
                for kk = 1:Matrices.n_elem(jj)
                    if isfield(Matrices.RenderPoints_RB,'ExtraInd')
                        extra_parts_ind = find(Matrices.RenderPoints_RB.ExtraInd(:,2)==kk);
                    else
                        extra_parts_ind = [];
                    end
                    ind1 = [1:6] + (kk-1)*6;
                    ind3 = [1:3] + (kk-1)*3;
                    ind4 = [1:4] + (kk-1)*4;
                    CGB = Quat2Rot(vectors_q{jj}(ind4,ii));
                    coords_upper = CGB*[Matrices.RenderPoints{jj}.x_upper(kk,:);
                                        Matrices.RenderPoints{jj}.y_upper(kk,:);
                                        Matrices.RenderPoints{jj}.z_upper(kk,:)] + repmat(vectors_p{jj}(ind3,ii),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                    coords_lower = CGB*[Matrices.RenderPoints{jj}.x_lower(kk,:);
                                        Matrices.RenderPoints{jj}.y_lower(kk,:);
                                        Matrices.RenderPoints{jj}.z_lower(kk,:)] + repmat(vectors_p{jj}(ind3,ii),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                    patch('XData',coords_upper(1,:),'YData',coords_upper(2,:),'ZData',coords_upper(3,:),'FaceColor','white')
                    hold on
                    patch('XData',coords_lower(1,:),'YData',coords_lower(2,:),'ZData',coords_lower(3,:),'FaceColor','white')
                    if ~isempty(extra_parts)&&~isempty(extra_parts_ind)
                        for ee = 1:size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,1)
                            coords_rot = CGB*[Matrices.RenderPoints_RB.Extra{extra_parts}.x(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.y(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.z(ee,:)] + repmat(vectors_p{extra_parts}(ind3,ii),1,size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,2));
                            X_out(ee,:) = coords_rot(1,:);
                            Y_out(ee,:) = coords_rot(2,:);
                            Z_out(ee,:) = coords_rot(3,:);
                        end
                        surf(X_out,Y_out,Z_out,'FaceColor','w')
                    end
                    if ghost_flag
                        CGB = Quat2Rot(x_dyn.x_q{jj}(ind4));
                        coords_upper = CGB*[Matrices.RenderPoints{jj}.x_upper(kk,:);
                                            Matrices.RenderPoints{jj}.y_upper(kk,:);
                                            Matrices.RenderPoints{jj}.z_upper(kk,:)] + repmat(x_dyn.x_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                        coords_lower = CGB*[Matrices.RenderPoints{jj}.x_lower(kk,:);
                                            Matrices.RenderPoints{jj}.y_lower(kk,:);
                                            Matrices.RenderPoints{jj}.z_lower(kk,:)] + repmat(x_dyn.x_p{jj}(ind3),1,size(Matrices.RenderPoints{jj}.x_upper,2));
                        patch('XData',coords_upper(1,:),'YData',coords_upper(2,:),'ZData',coords_upper(3,:),'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        hold on
                        patch('XData',coords_lower(1,:),'YData',coords_lower(2,:),'ZData',coords_lower(3,:),'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        if ~isempty(extra_parts)&&~isempty(extra_parts_ind)
                            for ee = 1:size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,1)
                                coords_rot = CGB*[Matrices.RenderPoints_RB.Extra{extra_parts}.x(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.y(ee,:);Matrices.RenderPoints_RB.Extra{extra_parts}.z(ee,:)] + repmat(x_dyn.x_p{extra_parts}(ind3),1,size(Matrices.RenderPoints_RB.Extra{extra_parts}.x,2));
                                X_out(ee,:) = coords_rot(1,:);
                                Y_out(ee,:) = coords_rot(2,:);
                                Z_out(ee,:) = coords_rot(3,:);
                            end
                            surf(X_out,Y_out,Z_out,'FaceColor','white','FaceAlpha',0.25,'EdgeAlpha',0.25)
                        end
                    end
                end
            end
            axis equal
            if ii == 1
                axes_props = get(gcf,'children');
                xlim_init  = axes_props.XLim;
                ylim_init  = axes_props.YLim;
                zlim_init  = axes_props.ZLim;
            end
            set(gca,'XLim',xlim_init+[-1 1],'YLim',ylim_init+[-1 1],'ZLim',zlim_init+[-2 2])
%             view([-90 0]) % Front View
            title(['Mode ',num2str(Mode2Plot),', Frequency = ',num2str(freq(Mode2Plot)/2/pi),'Hz (Damping = ',num2str(damp(Mode2Plot)),')'])
            hold off
            pause(0.001)
            if vid_flag == 2
                F(ii) = getframe(gcf);%,rect
                writeVideo(writerObj,F(ii));
            end
        end
    end
    pause(0.001)
    if vid_flag == 2
        close(writerObj)
    end
end