function PlotAircraftv2(x_in,Matrices,aeroFlag,cgFlag,tStep, ha)

if nargin < 6
    hf = figure;
    ha = axes('Parent', hf);
end

ha.NextPlot = 'add';
ha.Box      = 'on';

x_in = initDynSim(x_in,Matrices,Sim);
if isempty(tStep)
    tStep = 1;
end
if isempty(aeroFlag)
    aeroFlag = 0;
end
if isempty(cgFlag)
    cgFlag = 0;
end

if cgFlag
    for jj = 1:length(Matrices.n_elem)
        x{jj} = x_in{jj}(:,tStep);
    end
    MassProps = calculateCoM(x,Matrices,1);
end

struct_scale_fact = 3;
aero_scale_fact   = 3;

fprintf('\nPlotting Aircraft and Reference Frames...\n')
tic
CGa = Matrices.CGa;
for jj = 1:length(Matrices.n_elem)
    [coords ,CaB ] = strains2coords_all(x_in.x_f,Matrices,jj);
    for ii = 1:length(coords)
        coords_out{jj}(:,ii) = coords(:,:,ii);
    end
    if aeroFlag
        for ii = 1:size(Matrices.s_mp_ind{jj},2)
            ind  = [1:4] + (ii-1)*4;
            ind2 = [1:3] + (ii-1)*3;
            ind3 = [1:6] + (ii-1)*6;
            le = CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii))*Matrices.le{jj}(:,ii);
            te = CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii))*Matrices.te{jj}(:,ii);
            cp = CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii))*Matrices.cp{jj}(:,ii);
            scatter3(x_in.x_p{jj}(ind2(1),:)+le(1),x_in.x_p{jj}(ind2(2),:)+le(2),x_in.x_p{jj}(ind2(3),:)+le(3),'r.')
            hold on
            scatter3(x_in.x_p{jj}(ind2(1),:)+te(1),x_in.x_p{jj}(ind2(2),:)+te(2),x_in.x_p{jj}(ind2(3),:)+te(3),'r.')
            ex = CGa*CaB(:,1,Matrices.s_mp_ind{jj}(ii))/struct_scale_fact;
            ey = CGa*CaB(:,2,Matrices.s_mp_ind{jj}(ii))/struct_scale_fact;
            ez = CGa*CaB(:,3,Matrices.s_mp_ind{jj}(ii))/struct_scale_fact;
            plot3([x_in.x_p{jj}(ind2(1),:) x_in.x_p{jj}(ind2(1),:)+ex(1)],...
                  [x_in.x_p{jj}(ind2(2),:) x_in.x_p{jj}(ind2(2),:)+ex(2)],...
                  [x_in.x_p{jj}(ind2(3),:) x_in.x_p{jj}(ind2(3),:)+ex(3)],'g')
            plot3([x_in.x_p{jj}(ind2(1),:) x_in.x_p{jj}(ind2(1),:)+ey(1)],...
                  [x_in.x_p{jj}(ind2(2),:) x_in.x_p{jj}(ind2(2),:)+ey(2)],...
                  [x_in.x_p{jj}(ind2(3),:) x_in.x_p{jj}(ind2(3),:)+ey(3)],'g')
            plot3([x_in.x_p{jj}(ind2(1),:) x_in.x_p{jj}(ind2(1),:)+ez(1)],...
                  [x_in.x_p{jj}(ind2(2),:) x_in.x_p{jj}(ind2(2),:)+ez(2)],...
                  [x_in.x_p{jj}(ind2(3),:) x_in.x_p{jj}(ind2(3),:)+ez(3)],'g')
            CaA0 = CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii))*Matrices.CBA0{jj}(:,:,ii);
            ex = CaA0(:,1)/aero_scale_fact;
            ey = CaA0(:,2)/aero_scale_fact;
            ez = CaA0(:,3)/aero_scale_fact;
            plot3([x_in.x_p{jj}(ind2(1),:)+cp(1) x_in.x_p{jj}(ind2(1),:)+ex(1)+cp(1)],...
                  [x_in.x_p{jj}(ind2(2),:)+cp(2) x_in.x_p{jj}(ind2(2),:)+ex(2)+cp(2)],...
                  [x_in.x_p{jj}(ind2(3),:)+cp(3) x_in.x_p{jj}(ind2(3),:)+ex(3)+cp(3)],'r')
            plot3([x_in.x_p{jj}(ind2(1),:)+cp(1) x_in.x_p{jj}(ind2(1),:)+ey(1)+cp(1)],...
                  [x_in.x_p{jj}(ind2(2),:)+cp(2) x_in.x_p{jj}(ind2(2),:)+ey(2)+cp(2)],...
                  [x_in.x_p{jj}(ind2(3),:)+cp(3) x_in.x_p{jj}(ind2(3),:)+ey(3)+cp(3)],'r')
            plot3([x_in.x_p{jj}(ind2(1),:)+cp(1) x_in.x_p{jj}(ind2(1),:)+ez(1)+cp(1)],...
                  [x_in.x_p{jj}(ind2(2),:)+cp(2) x_in.x_p{jj}(ind2(2),:)+ez(2)+cp(2)],...
                  [x_in.x_p{jj}(ind2(3),:)+cp(3) x_in.x_p{jj}(ind2(3),:)+ez(3)+cp(3)],'r')
        end  
    end
    plot3(x_in.x_p{jj}(1:3:end,:),x_in.x_p{jj}(2:3:end,:),x_in.x_p{jj}(3:3:end,:),'k')
    hold on
    grid on
    axis equal
end
M_cg_rb_ = CGa*Matrices.M_cg_rb;
M_cg_tot = CGa*MassProps.cg;
scatter3(M_cg_rb_(1),M_cg_rb_(2),M_cg_rb_(3),'bx')
scatter3(M_cg_tot(1),M_cg_tot(2),M_cg_tot(3),'bo')
%     close all
toc