function [MassProps] = calculateCoM(x_f,Matrices,TextFlag)

M      = zeros(6);
M_memb = cell(1,length(Matrices.n_elem));

[coords ,CaB ] = strains2coords_all(x_f,Matrices);

for jj = 1:length(Matrices.n_elem)    
    M_memb{jj}     = zeros(6);
    for ii = 1:size(Matrices.s_mp_ind{jj},2)
        ds  = Matrices.s{jj}(ii+1) - Matrices.s{jj}(ii);
        M_e = Matrices.Aa2{jj}(:,:,ii)*ds;
        Ra_skew    = [eye(3) -skew(coords{jj}(:,:,Matrices.s_mp_ind{jj}(ii))+Matrices.p0{jj}); zeros(3) eye(3)];
        CaB_tot    = [CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) zeros(3); zeros(3) CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii))];
        M_memb{jj} = M_memb{jj} + Ra_skew'*CaB_tot*M_e*CaB_tot'*Ra_skew;
    end
    for ii = 1:size(Matrices.s_node_ind{jj},2)
        M_node     = Matrices.M_pt{jj}(:,:,ii);
        Ra_skew    = [eye(3) -skew(coords{jj}(:,:,Matrices.s_node_ind{jj}(ii))+Matrices.p0{jj}); zeros(3) eye(3)];
        CaB_tot    = [CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii)) zeros(3); zeros(3) CaB{jj}(:,:,Matrices.s_node_ind{jj}(ii))];
        M_memb{jj} = M_memb{jj} + Ra_skew'*CaB_tot*M_node*CaB_tot'*Ra_skew;
    end
    M = M + M_memb{jj};
end

if isfield(Matrices,'M_rb')
    M = M + Matrices.M_rb;
end

MassProps.M      = M;
MassProps.M_memb = M_memb;
MassProps.m      = M(1,1);
MassProps.cg     = [M(2,6);-M(1,6);M(1,5)]/MassProps.m;
MassProps.Mcg    = M - [zeros(3) -skew(MassProps.cg)*MassProps.m; skew(MassProps.cg)*MassProps.m -skew(MassProps.cg)*skew(MassProps.cg)*MassProps.m];

if TextFlag
    fprintf('\n************************************')
    fprintf('\n*  TOTAL AIRCRAFT MASS PROPERTIES  *')
    fprintf('\n************************************')
    fprintf(['\nTotal Mass of Aircraft: ',num2str(MassProps.m),'kg'])
    fprintf('\n************************************')
    if length(Matrices.n_elem)>1
        for jj = 1:length(Matrices.n_elem)
            fprintf(['\nMass of Member ',num2str(jj),': ',num2str(MassProps.M_memb{jj}(1,1)),'kg'])
        end
    end
    fprintf('\n************************************')
    fprintf( '\nCentre of Gravity Location:'                       )
    fprintf(['\nX: ',num2str(MassProps.cg(1)),'m'])
    fprintf(['\nY: ',num2str(MassProps.cg(2)),'m'])
    fprintf(['\nZ: ',num2str(MassProps.cg(3)),'m'])
    fprintf('\n************************************\n\n')
end
