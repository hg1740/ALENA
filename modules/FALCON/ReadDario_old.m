% function someoutputstructure = ReadDario(aircraft_data,hp)

%% Define the nodes of the members

% %Map AEROFLEX part names to ALENA part names
% %   - - TODO : Dario to update the part names to match the names in the
% %   Framework
% %   - TODO : Eventually we want to get rid of this and replace it with
% %   the handle to the 'LiftingSurface'/'Bluff-Body' object.
% map = {'Wing_R', 'StbdWing' ; 'Wing_L', 'PortWing' ; 'StbdHTP_R', 'StbdHTP' ; 'StbdHTP_L', 'PortHTP'};
% partID = {aircraft_data.PartId.Part};
% for i = 1 : size(map, 1)
%     partID(ismember(partID, map{i, 1})) = map(i, 2);
% end
% [aircraft_data.PartId.Part] = deal(partID{:});

%Grab unique part names
allPartIDs = unique({aircraft_data.PartId(:).Part},'stable');

%Find the wing
WingIdx = find(ismember(allPartIDs,'PortWing'));

for ii = 1:length(allPartIDs)
    part_idx     = ismember({aircraft_data.PartId.Part},allPartIDs{ii});
    part         = aircraft_data.PartId(part_idx);
    part_nodeidx = ismember({part.Type},'GRID');
    part_node    = part(part_nodeidx);
    part_aeroidx = ismember({part.Type},'CAERO');
    part_aero    = part(part_aeroidx);
    node_IDs{ii} = part_node.data;
    if isempty(part_aero)
        aero_IDs{ii} = 0;
    else
        aero_IDs{ii} = part_aero.data;
    end
end

fuse_nodes    = [];

flap_length_pc_chord = [0.25;0.25;0.25;0.25;0.25];
% %% Plot Aircraft
% if ~exist('hp','var')
%     hf = figure;
%     ha = axes('Parent',hf);
% else
%     ha = axes('Parent',uicontainer('Parent',hp));
% end
% scatter3(ha,aircraft_data.Node.Coord(:,1),aircraft_data.Node.Coord(:,2),aircraft_data.Node.Coord(:,3),'k.')
% axis(ha,'equal')
% hold(ha,'on')
% for jj = 1:length(allPartIDs)
%     for ii = 1:length(node_IDs{jj})
%         scatter3(ha,aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),1),...
%                     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2),...
%                     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3),'r')
%     end
% end

%% Iterate through the members
for jj = 1:numel(node_IDs)
    fprintf(['Setting up properties for Member ',num2str(jj),'...\n'])
    start_node(jj) = node_IDs{jj}(1);
    %% Number of Elements
    n_elem(jj) = length(node_IDs{jj})-1;
    %% Beam Lengths
    for ii = 1:n_elem(jj)
        s_coord(:,ii) = aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),:)-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),:);
    end
    l(jj) = sum(sqrt(sum(s_coord.^2, 1)));
    %% Surface Type
    type{jj}   = 'LS';
    %% Dihedral, Twist, Sweep
    local_dihedral{jj} =  zeros(1,n_elem(jj)+1);
    local_twist{jj}    =  zeros(1,n_elem(jj)+1);
    local_sweep{jj}    =  zeros(1,n_elem(jj)+1);
    
    %% Beam s-coordinates
    if 1%jj == 1 || jj == 2
        s{jj}      = cumsum([0 sqrt(sum(s_coord.^2, 1))]);%[0 : l(jj)/n_elem(jj)   : l(jj)];
    else
        s{jj}      = [0 : l(jj)/n_elem(jj)   : l(jj)];
    end
    clear s_coord
    s2{jj}     = sort([s{jj},s{jj}(2:end-1)-1e-5]);
    s_mp{jj}   = (s{jj}(1:end-1)+s{jj}(2:end))/2;%/4
    s_aero{jj} = [0 : l(jj)/n_elem(jj)*2  : l(jj)];%*2*4/2/4
    s_aero_mp{jj} = (s_aero{jj}(1:end-1)+s_aero{jj}(2:end))/2;%(s{1}(1:end-1)+s{1}(2:end))/2;
    s_cust{jj} = [];%14.3;
    
    s_out{jj}  = unique([s{jj},s_mp{jj},s_cust{jj},s_aero{jj},s_aero_mp{jj}]);
    s_out2{jj} = sort([s_out{jj},s{jj}(2:end-1)-1e-5]);
    
    s_node_ind{jj}  = zeros(size(s{jj}));
    for ii = 1:length(s{jj})
        s_node_ind{jj}(ii) = find(s{jj}(ii)==s_out{jj});
    end
    s_node_ind2{jj} = zeros(size(s2{jj}));
    for ii = 1:length(s2{jj})
        s_node_ind2{jj}(ii) = find(s2{jj}(ii)==s_out2{jj});
    end
    s_mp_ind{jj}    = zeros(size(s_mp{jj}));
    for ii = 1:length(s_mp{jj})
        s_mp_ind{jj}(ii) = find(s_mp{jj}(ii)==s_out{jj});
    end
    s_aero_ind{jj}  = zeros(size(s_aero{jj}));
    for ii = 1:length(s_aero{jj})
        s_aero_ind{jj}(ii) = find(s_aero{jj}(ii)==s_out{jj});
    end
    s_aero_mp_ind{jj} = zeros(size(s_aero_mp{jj}));
    for ii = 1:length(s_aero_mp{jj})
        s_aero_mp_ind{jj}(ii) = find(s_aero_mp{jj}(ii)==s_out{jj});
    end
    %% Beam Properties
    p{jj}          = zeros((n_elem(jj)+1)*3,1);
    p{jj}(1:3:end) = s{jj};

    p0{jj}   = aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(1),:)' - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{WingIdx}(1),:)';
    for ii = 1:n_elem(jj)
        %% Structural Setup        
        find1 = find(aircraft_data.Beam.Conn(:,1)==node_IDs{jj}(ii  ));
        find2 = find(aircraft_data.Beam.Conn(:,2)==node_IDs{jj}(ii  ));
        find3 = find(aircraft_data.Beam.Conn(:,1)==node_IDs{jj}(ii+1));
        find4 = find(aircraft_data.Beam.Conn(:,2)==node_IDs{jj}(ii+1));
        
        ID    = aircraft_data.Beam.ID( mode([find1',find2',find3',find4']));
        PID   = aircraft_data.Beam.PID(aircraft_data.Beam.ID==ID);
        
        if ii == 1
            ex = aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(2),:)'-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(1),:)';
            ey_tmp = aircraft_data.Beam.Orient(aircraft_data.Beam.ID==ID,:)';
            ez = cross(ex,ey_tmp);
            ey = cross(ez,ex);
            CaB0{jj} = [ex/norm(ex),ey/norm(ey),ez/norm(ez)];
        end
        
        ex_beam = aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),:)'-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),:)';
        ey_tmp = aircraft_data.Beam.Orient(aircraft_data.Beam.ID==ID,:)';
        ez = cross(ex_beam,ey_tmp);
        ey = cross(ez,ex_beam);
        CaBi{jj}(:,:,ii) = [ex_beam/norm(ex_beam),ey/norm(ey),ez/norm(ez)];
        if ii == 1
            CBB{jj}(:,:,ii) = (CaB0{jj}'          *CaBi{jj}(:,:,ii))';
        else
            CBB{jj}(:,:,ii) = (CaBi{jj}(:,:,ii-1)'*CaBi{jj}(:,:,ii))';
        end

        E     = aircraft_data.Mat.E(aircraft_data.PBeam.Mat(aircraft_data.PBeam.ID==PID));
        G     = aircraft_data.Mat.G(aircraft_data.PBeam.Mat(aircraft_data.PBeam.ID==PID));
        
        A        = mean(aircraft_data.PBeam.A(aircraft_data.PBeam.ID==PID).data);
        I1       = mean(aircraft_data.PBeam.I(aircraft_data.PBeam.ID==PID).data(2,:));
        I2       = mean(aircraft_data.PBeam.I(aircraft_data.PBeam.ID==PID).data(1,:));
        I12      = mean(aircraft_data.PBeam.I(aircraft_data.PBeam.ID==PID).data(3,:));
        J_       = mean(aircraft_data.PBeam.J(aircraft_data.PBeam.ID==PID).data);
        offset_a = aircraft_data.Beam.Offset(aircraft_data.PBeam.ID==PID,1:3);
        offset_b = aircraft_data.Beam.Offset(aircraft_data.PBeam.ID==PID,7:9);
        R_shear  = (offset_a'+offset_b')/2;
        
        ds       = s{jj}(ii+1)-s{jj}(ii);
        ex       = offset_b'+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),:)'-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),:)'-offset_a';
        ey_tmp   = aircraft_data.Beam.Orient(aircraft_data.Beam.ID==ID,:)';
        ez       = cross(ex,ey_tmp);
        ey       = cross(ez,ex);
        CBBbar   = CaBi{jj}(:,:,ii)'*[ex/norm(ex),ey/norm(ey),ez/norm(ez)];
        R_tilde  = [eye(3) -skew(CaBi{jj}(:,:,ii)'*R_shear);zeros(3) eye(3)];
        dsprime  = norm(ex);
        length_fact = ds/dsprime;
        
        K{jj}(:,:,ii)  = eye(6);
        K{jj}(1,1,ii)  = E*A;
        K{jj}(2,2,ii)  = G*A;
        K{jj}(3,3,ii)  = G*A;
        K{jj}(4,4,ii)  = G*J_;
        K{jj}(5,5,ii)  = E*I1;
        K{jj}(6,6,ii)  = E*I2;
        K{jj}(5,6,ii)  = E*I12; %TODO - When working with the HARTEN model we need to multiply this by 0.
        K{jj}(6,5,ii)  = E*I12; %TODO - When working with the HARTEN model we need to multiply this by 0.
        K{jj}(:,:,ii)  = R_tilde'*blkdiag(CBBbar,CBBbar)*K{jj}(:,:,ii)*blkdiag(CBBbar',CBBbar')*R_tilde*length_fact;
        
        CFFbar{jj}(:,:,ii) = inv(R_tilde'*blkdiag(CBBbar,CBBbar));
        
        C{jj}(:,:,ii)  = inv(K{jj}(:,:,ii));
        
        h = sqrt(sqrt(I1/I2)*A);
        b = A/h;
        
        rho = aircraft_data.Mat.Rho(aircraft_data.PBeam.Mat(aircraft_data.PBeam.ID==PID));
        
        M_pt       = zeros(6);
        
        % Search for lumped masses associated with the end element
        CONM2_list = find(aircraft_data.Conm2.Node==node_IDs{jj}(ii+1));
        for mm = 1:length(CONM2_list)
            M_pt = M_pt + blkdiag(CaBi{jj}(:,:,ii)',CaBi{jj}(:,:,ii)')*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(CONM2_list(mm),:));zeros(3)  eye(3)]*blkdiag(CaBi{jj}(:,:,ii),CaBi{jj}(:,:,ii));
        end
        
        M_pt_store{jj}(:,:,ii+1) = zeros(6);%M_pt;
        if ii == 1
            M_pt_store{jj}(:,:,1) = zeros(6);
        end
        
        m_pt0     =    M_pt(1,1);
        if m_pt0 == 0
            M_pt0 = zeros(6);
            cg_pt = [0;0;0];
        else
            cg_pt     = -[-M_pt(2,6);M_pt(1,6);-M_pt(1,5)]/m_pt0;
            M_pt0     =    M_pt - [zeros(3) -skew(cg_pt)*m_pt0; skew(cg_pt)*m_pt0 -skew(cg_pt)*skew(cg_pt)*m_pt0];
            cg_pt0    =    cg_pt;
        end
        
        J_pt0 = M_pt(4:6,4:6);
        
        m{jj}(ii)      = rho*A + m_pt0/ds;
        i1             = rho*A/12*(b^2 + h^2) ;
        i2             = rho*A*(b^2 + (s{jj}(ii+1)-s{jj}(ii))^2)/12;%1e-5;%
        i3             = rho*A*(h^2 + (s{jj}(ii+1)-s{jj}(ii))^2)/12;%1e-5;%
        J{jj}(:,:,ii)  = diag([i1,i2,i3]) + J_pt0/ds;
        cg{jj}(:,:,ii) = [0;0;0] + cg_pt0;
        
        % TODO - Remove hardcoded indices for the different members. We
        % should be looking for the 'VTP' or the 'WING' etc.
        if ~strcmpi(allPartIDs{jj},'vtp') % Main Wings
            %% Aero Setup
            clear chord_vec z_vec y_vec x_vec tw_vec
            for kk = 1:length(aero_IDs{jj})
                chord_vec(kk) = aircraft_data.Aero.geo.c(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                x_vec(kk)     = aircraft_data.Aero.geo.startx(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                y_vec(kk)     = aircraft_data.Aero.geo.starty(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                z_vec(kk)     = aircraft_data.Aero.geo.startz(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                tw_vec(kk)    = aircraft_data.Aero.geo.TW(aircraft_data.Aero.ID==aero_IDs{jj}(kk),:,1);
            end
            
            % TODO - Remove this from the code fudge needed to account fro
            % strip theory 
            tw_vec = smooth(tw_vec,20);
            
            aerotwist{jj}(ii) = interp1(y_vec,tw_vec, (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2,'linear','extrap');
            chord{jj}(ii)    = interp1(y_vec,chord_vec, (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2,'linear','extrap');
            le_a             = [interp1(y_vec,x_vec,    (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2,'linear','extrap');...
                                                        (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2;...
                                interp1(y_vec,z_vec,    (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2,'linear','extrap')];
            le{jj}(:,ii)     = CaBi{jj}(:,:,ii)'*[le_a(1) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),1)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),1))/2;...
                                                  le_a(2) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2;...
                                                  le_a(3) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2];
            ez = cross(le{jj}(:,ii),[1;0;0]);
            ex = cross(le{jj}(:,ii),ez);
            
            % Remove the x contribution 
            ex_beam(1) = 0;
            ex_beam = ex_beam/norm(ex_beam);
            
            %TODO - Robbie & Dario to fix the twist implementation.
            %Possibly an issue with sign convention between the port &
            %starboard wing.
            CBB_tw{jj}(:,:,ii) =  eye(3);%CaBi{jj}(:,:,ii)'*RotVecHinge(ex_beam(1),ex_beam(2),ex_beam(3),(pi/180)*aerotwist{jj}(ii))*CaBi{jj}(:,:,ii);
            CBA0{jj}(:,:,ii) = CBB_tw{jj}(:,:,ii)*[ex/norm(ex),le{jj}(:,ii)/norm(le{jj}(:,ii)),ez/norm(ez)];
            cp{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-0.25*chord{jj}(ii);0.0] + le{jj}(:,ii);
            te{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-1.00*chord{jj}(ii);0.0] + le{jj}(:,ii);
            tq{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-0.75*chord{jj}(ii);0.0] + le{jj}(:,ii);
        else  % VTP
            %% Aero Setup
            clear chord_vec z_vec y_vec x_vec tw_vec
            for kk = 1:length(aero_IDs{jj})
                chord_vec(kk) = aircraft_data.Aero.geo.c(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                x_vec(kk)     = aircraft_data.Aero.geo.startx(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                y_vec(kk)     = aircraft_data.Aero.geo.starty(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                z_vec(kk)     = aircraft_data.Aero.geo.startz(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
                tw_vec(kk)    = aircraft_data.Aero.geo.TW(aircraft_data.Aero.ID==aero_IDs{jj}(kk),:,1);

            end
            aerotwist{jj}(ii) = interp1(z_vec,tw_vec, (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2,'linear','extrap');

            chord{jj}(ii)    = interp1(z_vec,chord_vec, (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2,'linear','extrap');
            le_a             = [interp1(z_vec,x_vec,    (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2,'linear','extrap');...
                                interp1(z_vec,y_vec,    (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2,'linear','extrap');...
                                                        (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2];
            le{jj}(:,ii)     = CaBi{jj}(:,:,ii)'*[le_a(1) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),1)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),1))/2;...
                                                  le_a(2) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2;...
                                                  le_a(3) - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii+1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2];
            ez = cross(le{jj}(:,ii),[1;0;0]);
            ex = cross(le{jj}(:,ii),ez);
            CBA0{jj}(:,:,ii) = [ex/norm(ex),le{jj}(:,ii)/norm(le{jj}(:,ii)),ez/norm(ez)];
            cp{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-0.25*chord{jj}(ii);0.0] + le{jj}(:,ii);
            te{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-1.00*chord{jj}(ii);0.0] + le{jj}(:,ii);
            tq{jj}(:,ii)     = CBA0{jj}(:,:,ii)*[0.0;-0.75*chord{jj}(ii);0.0] + le{jj}(:,ii);
        end
        %% Member flaps
        flap_length      = flap_length_pc_chord(jj)*chord{jj}(ii);
        k                = (chord{jj}(ii)-flap_length)/chord{jj}(ii);
        theta_k          =  acos(1-2*k);
        dCLdd{jj}(:,ii)  = (2.00*(pi-theta_k)   + 2.00*sin(theta_k));
        dCMdd{jj}(:,ii)  = (0.25*sin(2*theta_k) - 0.50*sin(theta_k));
        
        xs{1}(ii)        = 0.0;
    end

    iaf.designation       = '0012';
    iaf.n                 = 10;
    iaf.HalfCosineSpacing = 1;
    iaf.wantFile          = 0;
    iaf.is_finiteTE       = 0;

    RenderPoints{jj}.y_upper = zeros(n_elem(jj),2*(iaf.n+1));
    RenderPoints{jj}.y_lower = zeros(n_elem(jj),2*(iaf.n+1));
    RenderPoints{jj}.z_upper = zeros(n_elem(jj),2*(iaf.n+1));
    RenderPoints{jj}.z_lower = zeros(n_elem(jj),2*(iaf.n+1));

    af                    = naca4gen(iaf);
    for ii = 1:n_elem(jj)+1
        %% Rendering Information       
        clear chord_vec z_vec y_vec x_vec
        for kk = 1:length(aero_IDs{jj})
            chord_vec(kk) = aircraft_data.Aero.geo.c(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
            x_vec(kk)     = aircraft_data.Aero.geo.startx(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
            y_vec(kk)     = aircraft_data.Aero.geo.starty(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
            z_vec(kk)     = aircraft_data.Aero.geo.startz(aircraft_data.Aero.ID==aero_IDs{jj}(kk));
        end
        
        if strcmpi(allPartIDs{jj},'vtp')
            renderchord =  interp1(z_vec,chord_vec, aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3),'linear','extrap');
            le_a        = [interp1(z_vec,x_vec,     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3),'linear','extrap');...
                           interp1(z_vec,y_vec,     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3),'linear','extrap');...
                                                    aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3)];
            z_upper(ii,:) = ones(size(af.xU'))*aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3);
            z_lower(ii,:) = ones(size(af.xU'))*aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3);
            y_upper(ii,:) = (1-af.xU')*renderchord + le_a(1);% - norm(le0);
            y_lower(ii,:) = (1-af.xL')*renderchord + le_a(1);% - norm(le0);
            x_upper(ii,:) = af.zU(end:-1:1)'*renderchord;% + le_a(3);
            x_lower(ii,:) = af.zL(end:-1:1)'*renderchord;% + le_a(3);
        else
            renderchord =  interp1(y_vec,chord_vec, aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2),'linear','extrap');
            le_a        = [interp1(y_vec,x_vec,     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2),'linear','extrap');...
                           aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2)                   ;...
                           interp1(y_vec,z_vec,     aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2),'linear','extrap')];
            x_upper(ii,:) = ones(size(af.xU'))*aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2);
            x_lower(ii,:) = ones(size(af.xU'))*aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2);
            y_upper(ii,:) = (1-af.xU')*renderchord + le_a(1);% - norm(le0);
            y_lower(ii,:) = (1-af.xL')*renderchord + le_a(1);% - norm(le0);
            z_upper(ii,:) = af.zU(end:-1:1)'*renderchord + le_a(3);
            z_lower(ii,:) = af.zL(end:-1:1)'*renderchord + le_a(3);
        end
        
        if ii > 1
            x_coord_upper = [y_upper(ii,:)';y_upper(ii-1,end:-1:1)'];
            y_coord_upper = [x_upper(ii,:)';x_upper(ii-1,end:-1:1)'];
            z_coord_upper = [z_upper(ii,:)';z_upper(ii-1,end:-1:1)'];
            %patch(ha,x_coord_upper,y_coord_upper,z_coord_upper,'white')
            x_coord_lower = [y_lower(ii,:)';y_lower(ii-1,end:-1:1)'];
            y_coord_lower = [x_lower(ii,:)';x_lower(ii-1,end:-1:1)'];
            z_coord_lower = [z_lower(ii,:)';z_lower(ii-1,end:-1:1)'];
            %patch(ha,x_coord_lower,y_coord_lower,z_coord_lower,'white')
            x_upper_B(ii-1,:) = x_coord_upper - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),1)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),1))/2;
            y_upper_B(ii-1,:) = y_coord_upper - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2;
            z_upper_B(ii-1,:) = z_coord_upper - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2;
            upper_coords = CaBi{jj}(:,:,ii-1)'*[x_upper_B(ii-1,:);y_upper_B(ii-1,:);z_upper_B(ii-1,:)];
            x_upper_B(ii-1,:) = upper_coords(1,:);
            y_upper_B(ii-1,:) = upper_coords(2,:);
            z_upper_B(ii-1,:) = upper_coords(3,:);
            x_lower_B(ii-1,:) = x_coord_lower - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),1)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),1))/2;
            y_lower_B(ii-1,:) = y_coord_lower - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),2)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),2))/2;
            z_lower_B(ii-1,:) = z_coord_lower - (aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii-1),3)+aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{jj}(ii),3))/2;
            lower_coords = CaBi{jj}(:,:,ii-1)'*[x_lower_B(ii-1,:);y_lower_B(ii-1,:);z_lower_B(ii-1,:)];
            x_lower_B(ii-1,:) = lower_coords(1,:);
            y_lower_B(ii-1,:) = lower_coords(2,:);
            z_lower_B(ii-1,:) = lower_coords(3,:);
        end
    end
    
    RenderPoints{jj}.x_upper = x_upper_B;
    RenderPoints{jj}.x_lower = x_lower_B;
    RenderPoints{jj}.y_upper = y_upper_B;% - norm(le0);
    RenderPoints{jj}.y_lower = y_lower_B;% - norm(le0);
    RenderPoints{jj}.z_upper = z_upper_B;
    RenderPoints{jj}.z_lower = z_lower_B;
    
    hF(1)  = figure;
    hAx = axes('Parent', hF(1), 'NextPlot', 'add');
    patch(hAx, x_upper_B, y_upper_B, z_upper_B, 'white');
    patch(hAx, x_lower_B, y_lower_B, z_lower_B, 'white');
    
    clear x_upper_B x_lower_B y_upper_B y_lower_B z_upper_B z_lower_B
    %% Members Aero Settings
    delta_flap{jj} = zeros(n_elem(jj),1);
    tau(jj)        = 1.0;
    lift_dist{jj}  = ones(size(s{jj}));
    %% Forces and initial conditions
    x{jj}     = zeros(n_elem(jj)*6,1);
    
    f               = zeros((n_elem(jj)+1)*6,1);
    f_a             = zeros((n_elem(jj)+1)*6,1);
    f_mp_pnt        = zeros( n_elem(jj)   *6,1);
    Matrices.f{jj}        =  f;
    Matrices.f_a{jj}      =  f_a;
    Matrices.f_mp_pnt{jj} =  f_mp_pnt;
    
    x2{1} = zeros(n_elem(jj)*12,1);
    f2    = zeros(n_elem(jj)*12,1);
    f_a2  = zeros(n_elem(jj)*12,1);
    Matrices.f2{jj}   =    f2;
    Matrices.f_a2{jj} =  f_a2;   
end
Matrices.grav_switch     = 1;

start_node = unique(start_node);
%% RB Data

% If there exists any nodes which are rigidly connected to any start nodes
% then lump them up onto the fuselage mass matrix. 

M_rb_fus = zeros(6);
for ii = 1:length(start_node)
    RB_CONM2_list    = find(aircraft_data.Conm2.Node==start_node(ii));
    M_cg_rb_fus      = aircraft_data.Node.Coord(aircraft_data.Node.ID==start_node(ii),:)-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{WingIdx}(1),:);
    for mm = 1:length(RB_CONM2_list)
        M_rb_fus         = M_rb_fus + [eye(3) zeros(3);skew(M_cg_rb_fus)  eye(3)]*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,RB_CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:)); zeros(3) eye(3)]*[eye(3) -skew(M_cg_rb_fus); zeros(3) eye(3)];
    end
end
AllRB_CONM2 = [];
% Iterate through the start nodes
if isfield(aircraft_data, 'RBE2')
    for ii = 1:length(start_node)
        % Find out if any of those start nodes are member of an Rbe2
        slave_nodes   = [aircraft_data.RBE2.IDS.data];
        Rebe2_list    = find(slave_nodes==start_node(ii));
        % Find out if there are any lumped masses connected to the master node
        RB_CONM2_list = find(aircraft_data.Conm2.Node==aircraft_data.RBE2.IDM(Rebe2_list));
        % if any RB_CONM2_list are contained in the the all list remove them
        [~,~,coin] = intersect(AllRB_CONM2,RB_CONM2_list);
        % Remove any nodes in the the AllRB_CONM2 from RB_CONM2_list
        RB_CONM2_list(coin) = [];
        % Store the conm2s so that they do not get repeated
        AllRB_CONM2  = [AllRB_CONM2;RB_CONM2_list];
        % Calculate the offset of the rigidly connected node to the start node
        M_cg_rb_fus   = aircraft_data.Node.Coord(aircraft_data.Node.ID==aircraft_data.RBE2.IDM(Rebe2_list),:)-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{WingIdx}(1),:);
        % Iterate through the conm2s and get the offsets
        for mm = 1:length(RB_CONM2_list)
            M_rb_fus         = M_rb_fus + [eye(3) zeros(3);skew(M_cg_rb_fus)  eye(3)]*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,RB_CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:)); zeros(3) eye(3)]*[eye(3) -skew(M_cg_rb_fus); zeros(3) eye(3)];
        end
    end
end

M_cg_full        = -[-M_rb_fus(2,6);M_rb_fus(1,6);-M_rb_fus(1,5)];%M_cg_rb_fus*M_rb_fus(1); 
Mf_tot           =    M_rb_fus(1);

Matrices.M_cg_rb =  M_cg_full/Mf_tot(1);
Matrices.M_rb    =  M_rb_fus;
% 
% theta_interp = [0:5:180];
% fuse_interp  = 0.5*(1-cos(theta_interp*pi/180));
% r_interp     = interp1(aircraft_data.BAero.geo.fs{2},aircraft_data.BAero.geo.Rs{2},fuse_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(fuse_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(2);
% 
% Z_shift      = interp1([0 5 35 44],[-0.1 0 0 1],X(:,1));
% Z            = Z + repmat(Z_shift,1,size(Z,2));
% 
% for ii = 1:size(X,1)
%     coords_rot = [X(ii,:);Y(ii,:);Z(ii,:)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X_out,Y_out,Z_out,'FaceColor','w')
% RenderPoints_RB.x = X_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),1);
% RenderPoints_RB.y = Y_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),2);
% RenderPoints_RB.z = Z_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),3);

%% Engines
% theta_interp = [0:30:180];
% eng_interp   = 0.5*(1-cos(theta_interp*pi/180));
% r_interp     = interp1(aircraft_data.BAero.geo.fs{3},aircraft_data.BAero.geo.Rs{3},eng_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(eng_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(3);
% clear X_out Y_out Z_out
% for ii = 1:size(X,1)
%     coords_rot = CaBi{1}(:,:,10)'*[X(ii,:) + aircraft_data.BAero.geo.ref_point(3,1) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),1);...
%                                    Y(ii,:) + aircraft_data.BAero.geo.ref_point(3,2) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),2);...
%                                    Z(ii,:) + aircraft_data.BAero.geo.ref_point(3,3) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),3)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X+aircraft_data.BAero.geo.ref_point(3,1),Y+aircraft_data.BAero.geo.ref_point(3,2),Z+aircraft_data.BAero.geo.ref_point(3,3),'FaceColor','w')
% RenderPoints_RB.Extra{1}.x = X_out;
% RenderPoints_RB.Extra{1}.y = Y_out;
% RenderPoints_RB.Extra{1}.z = Z_out;
% r_interp     = interp1(aircraft_data.BAero.geo.fs{1},aircraft_data.BAero.geo.Rs{1},eng_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(eng_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(1);
% clear X_out Y_out Z_out
% for ii = 1:size(X,1)
%     coords_rot = CaBi{2}(:,:,10)'*[X(ii,:) + aircraft_data.BAero.geo.ref_point(1,1) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),1);...
%                                    Y(ii,:) + aircraft_data.BAero.geo.ref_point(1,2) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),2);...
%                                    Z(ii,:) + aircraft_data.BAero.geo.ref_point(1,3) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),3)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X+aircraft_data.BAero.geo.ref_point(1,1),Y+aircraft_data.BAero.geo.ref_point(1,2),Z+aircraft_data.BAero.geo.ref_point(1,3),'FaceColor','w')
% RenderPoints_RB.Extra{2}.x = X_out;
% RenderPoints_RB.Extra{2}.y = Y_out;
% RenderPoints_RB.Extra{2}.z = Z_out;
% RenderPoints_RB.ExtraInd   = [1 10;
%                               2 10];
%% Aero Distributions
% load('\\ads\filestore\Engineering\Research\Projects\AWI\Models\Straight_Tapered\NoTwist\Cl_alpha.mat')
% load('\\ads\filestore\Engineering\Research\Projects\AWI\Models\Scitech2018Straight_Tapered\Cl_alphaAR10.mat')
% PanelArea = PanelArea';
% 
% wing1CL          = [Cl_alpha( 1:22);(Cl_alpha(23:2:33 ).*PanelArea(23:2:33 )+Cl_alpha(24:2:34 ).*PanelArea(24:2:34 ))./(PanelArea(23:2:33 )+PanelArea(24:2:34 ));Cl_alpha(35:38)];
% wing2CL          = [Cl_alpha(39:60);(Cl_alpha(61:2:71 ).*PanelArea(61:2:71 )+Cl_alpha(62:2:72 ).*PanelArea(62:2:72 ))./(PanelArea(61:2:71 )+PanelArea(62:2:72 ));Cl_alpha(73:76)];
% wing3CL          = [                (Cl_alpha(77:2:83 ).*PanelArea(77:2:83 )+Cl_alpha(78:2:84 ).*PanelArea(78:2:84 ))./(PanelArea(77:2:83 )+PanelArea(78:2:84 ));Cl_alpha(85)   ];
% wing4CL          = [                (Cl_alpha(86:2:92 ).*PanelArea(86:2:92 )+Cl_alpha(87:2:93 ).*PanelArea(87:2:93 ))./(PanelArea(86:2:92 )+PanelArea(87:2:93 ));Cl_alpha(94)   ];
% wing5CL          = [                (Cl_alpha(95:2:101).*PanelArea(95:2:101)+Cl_alpha(96:2:102).*PanelArea(96:2:102))./(PanelArea(95:2:101)+PanelArea(96:2:102));Cl_alpha(103)  ];
% 
% rootCL1          = interp1(s_mp{1},wing1CL/2/pi,s{1}(1),'linear','extrap');
% rootCL2          = interp1(s_mp{2},wing2CL/2/pi,s{2}(1),'linear','extrap');
% rootCL3          = interp1(s_mp{3},wing3CL/2/pi,s{3}(1),'linear','extrap');
% rootCL4          = interp1(s_mp{4},wing4CL/2/pi,s{4}(1),'linear','extrap');
% rootCL5          = interp1(s_mp{5},wing5CL/2/pi,s{5}(1),'linear','extrap');
% 
% n_wing           = length(wing1CL);
% Awing            = zeros(n_wing+2,n_wing+1);
% Bwing            = zeros(n_wing+1,n_wing);
% Cwing            = zeros(n_wing+1, 1);
% 
% Awing(1:n_wing+1,1:n_wing+1) =                                eye(n_wing+1);
% Awing(2:n_wing+2,1:n_wing+1) = Awing(2:n_wing+2,1:n_wing+1) + eye(n_wing+1);
% Awing                        = Awing(1:n_wing+1,1:n_wing+1)                ;
% Bwing(2:n_wing+1,1:n_wing)   =                              2*eye(n_wing)  ;
% Cwing(1)                     =                              1              ;
% 
% lift_dist{1}     = Awing\(Bwing*wing1CL/2/pi + Cwing*rootCL1);
% lift_dist{2}     = Awing\(Bwing*wing2CL/2/pi + Cwing*rootCL2);
% 
% n_htp            = length(wing3CL);
% Ahtp             = zeros(n_htp+2,n_htp+1);
% Bhtp             = zeros(n_htp+1,n_htp);
% Chtp             = zeros(n_htp+1,1);
% 
% Ahtp(1:n_htp+1,1:n_htp+1) =                             eye(n_htp+1);
% Ahtp(2:n_htp+2,1:n_htp+1) = Ahtp(2:n_htp+2,1:n_htp+1) + eye(n_htp+1);
% Ahtp                      = Ahtp(1:n_htp+1,1:n_htp+1)               ;
% Bhtp(2:n_htp+1,1:n_htp)   =                           2*eye(n_htp)  ;
% Chtp(1)                   =                           1             ;
% 
% lift_dist{3}     = Ahtp\(Bhtp*wing3CL/2/pi + Chtp*rootCL3);
% lift_dist{4}     = Ahtp\(Bhtp*wing4CL/2/pi + Chtp*rootCL4);
% lift_dist{5}     = Ahtp\(Bhtp*wing5CL/2/pi + Chtp*rootCL5);