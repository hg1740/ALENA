function beam_model = reflectmodel(beam_model,Aircraft_param,SizingInput)

% Reflect the wing
if ~isempty(Aircraft_param.Wing)
    [beam_model,beam_model.PartIDs.Wing] = reflectwing(beam_model,beam_model.PartIDs.Wing,beam_model.Aero.Wing);
end

% Reflect the HTP
if isfield(Aircraft_param,'StbdHTP')
    if ~isempty(Aircraft_param.StbdHTP)
        [beam_model,beam_model.PartIDs.StbdHTP] = reflectwing(beam_model,beam_model.PartIDs.StbdHTP,beam_model.Aero.StbdHTP);
    end
end

nodew_index = [];
for j = 1:length(beam_model.PartIDs.Wing.Nodes)
    index = find(beam_model.Node.ID == beam_model.PartIDs.Wing.Nodes(j)) ;
    nodew_index = [nodew_index,index];
end

if isfield(Aircraft_param,'StbdHTP')
    if ~isempty(Aircraft_param.StbdHTP)
        nodeh_index = [];
        for j = 1:length(beam_model.PartIDs.StbdHTP.Nodes)
            index = find(beam_model.Node.ID == beam_model.PartIDs.StbdHTP.Nodes(j)) ;
            nodeh_index = [nodeh_index,index];
        end
    end
end

if isfield(Aircraft_param,'VTP')
    if ~isempty(Aircraft_param.VTP)
        nodev_index = [];
        for j = 1:length(beam_model.PartIDs.VTP.Nodes)
            index = find(beam_model.Node.ID == beam_model.PartIDs.VTP.Nodes(j)) ;
            nodev_index = [nodev_index,index];
        end
    end
end

% Reflect the fuel
if SizingInput.AddFuel == 1
    ncom2 = beam_model.Info.ncom2;
    count = 0;
    
    % Fuel Indices
    conm2fuel_index = [];
    for j = 1:length(beam_model.PartIDs.Fuel.Conm2IDs)
        index = find(beam_model.Conm2.ID == beam_model.PartIDs.Fuel.Conm2IDs(j)) ;
        conm2fuel_index = [conm2fuel_index,index];
    end
    
    for iref = conm2fuel_index
        count=count+1;
        beam_model.Conm2.ID(count+ncom2)     = beam_model.Conm2.ID(iref) + 5000;
        if beam_model.Conm2.Node(iref) == beam_model.Node.ID(nodew_index(1))
            beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref);
        else
            beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref) + 5000;
        end
        beam_model.Conm2.CID(count+ncom2)    = beam_model.Conm2.CID(iref);
        beam_model.Conm2.M(:,:,count+ncom2)  = beam_model.Conm2.M(:,:,iref);
        beam_model.Conm2.Offset(count+ncom2,:) = beam_model.Conm2.Offset(iref,:)*diag([1,-1,1]);
        beam_model.PartIDs.Fuel.Conm2IDs = [beam_model.PartIDs.Fuel.Conm2IDs,beam_model.Conm2.ID(count+ncom2)];
    end
    beam_model.Info.ncom2 = ncom2+count;
end

% Reflect the Engine
if isfield(beam_model.PartIDs,'Engine')
    ncom2 = beam_model.Info.ncom2;
    count = 0;
    
    if strcmp(Aircraft_param.Engine.Parent,'Wing')
        nodeidx = nodew_index;
    elseif strcmp(Aircraft_param.Engine.Parent,'StbdHTP')
        nodeidx = nodeh_index;
    elseif strcmp(Aircraft_param.Engine.Parent,'VTP')
        nodeidx = nodev_index;
    end
    
    % Engine Indices
    conm2eng_index = [];
    for j = 1:length(beam_model.PartIDs.Engine.Conm2IDs)
        index = find(beam_model.Conm2.ID == beam_model.PartIDs.Engine.Conm2IDs(j)) ;
        conm2eng_index = [conm2eng_index,index];
    end
    
    for iref = conm2eng_index
        count=count+1;
        beam_model.Conm2.ID(count+ncom2)     = beam_model.Conm2.ID(iref) + 5000;
        if beam_model.Conm2.Node(iref) == beam_model.Node.ID(nodeidx(1))
            beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref);
        else
            beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref) + 5000;
        end
        beam_model.Conm2.CID(count+ncom2)    = beam_model.Conm2.CID(iref);
        beam_model.Conm2.M(:,:,count+ncom2)  = beam_model.Conm2.M(:,:,iref);
        beam_model.Conm2.Offset(count+ncom2,:) = beam_model.Conm2.Offset(iref,:)*diag([1,-1,1]);
        beam_model.PartIDs.Engine.Conm2IDs = [beam_model.PartIDs.Engine.Conm2IDs,beam_model.Conm2.ID(count+ncom2)];
    end
    beam_model.Info.ncom2 = ncom2+count;
    % Reflect the CAEROB entry
    
    nbody =  length(beam_model.Aero.body.ID);
    
    EngineBody = beam_model.PartIDs.Engine.BodyIdx;
    beam_model.Aero.body.ID(nbody+1) = beam_model.Aero.body.ID(EngineBody)+1;
    beam_model.Aero.body.CP(nbody+1) = beam_model.Aero.body.CP(EngineBody);
    beam_model.Aero.body.SET(nbody+1) = beam_model.Aero.body.SET(EngineBody)+1;
    beam_model.Aero.body.geo.fs{nbody+1} = beam_model.Aero.body.geo.fs{EngineBody};
    beam_model.Aero.body.geo.Rs{nbody+1} = beam_model.Aero.body.geo.Rs{EngineBody};
    beam_model.Aero.body.geo.Nelem(nbody+1) = beam_model.Aero.body.geo.Nelem(EngineBody);
    beam_model.Aero.body.geo.L(nbody+1) = beam_model.Aero.body.geo.L(EngineBody);
    beam_model.Aero.body.geo.ref_point(nbody+1,:) = [beam_model.Aero.body.geo.ref_point(EngineBody,1),...
        -beam_model.Aero.body.geo.ref_point(EngineBody,2),beam_model.Aero.body.geo.ref_point(EngineBody,3)];
    beam_model.Aero.body.Colour(nbody+1,:) = beam_model.Aero.body.Colour(EngineBody,:);
    beam_model.Info.nbaero = beam_model.Info.nbaero + 1;
    beam_model.Aero.body.lattice.Elem.Node{nbody + 1} = [beam_model.Aero.body.lattice.Elem.Node{EngineBody}(:,1),...
        -beam_model.Aero.body.lattice.Elem.Node{EngineBody}(:,2),beam_model.Aero.body.lattice.Elem.Node{EngineBody}(:,3)];
    beam_model.Aero.body.lattice.Elem.Conn{nbody + 1} = beam_model.Aero.body.lattice.Elem.Conn{EngineBody};
    %        beam_model.Info.nbody = nbody +1;
end

% Reflect the Fuselage
ncom2 = beam_model.Info.ncom2;
count = 0;

% Fuselage Indices
conm2fuse_index = [];
if isfield(beam_model.PartIDs,'Fuselage')
    for j = 1:length(beam_model.PartIDs.Fuselage.Conm2IDs)
        index = find(beam_model.Conm2.ID == beam_model.PartIDs.Fuselage.Conm2IDs(j)) ;
        conm2fuse_index = [conm2fuse_index,index];
    end
end

for iref = conm2fuse_index
    count=count+1;
    beam_model.Conm2.ID(count+ncom2)     = beam_model.Conm2.ID(iref) + 5000;
    
    if beam_model.Node.Coord(beam_model.Node.ID == beam_model.Conm2.Node(iref),2) == 0
        beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref);
    else
        beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref) + 5000;
    end
    beam_model.Conm2.CID(count+ncom2)    = beam_model.Conm2.CID(iref);
    beam_model.Conm2.M(:,:,count+ncom2)  = beam_model.Conm2.M(:,:,iref);
    beam_model.Conm2.Offset(count+ncom2,:) = beam_model.Conm2.Offset(iref,:)*diag([1,-1,1]);
    beam_model.PartIDs.Fuselage.Conm2IDs = [beam_model.PartIDs.Fuselage.Conm2IDs,beam_model.Conm2.ID(count+ncom2)];
end
beam_model.Info.ncom2 = ncom2+count;

% Reflect the Payload
ncom2 = beam_model.Info.ncom2;
count = 0;

% Fuselage Indices
conm2pay_index = [];
for j = 1:length(beam_model.PartIDs.Payload.Conm2IDs)
    index = find(beam_model.Conm2.ID == beam_model.PartIDs.Payload.Conm2IDs(j)) ;
    conm2pay_index = [conm2pay_index,index];
end

for iref = conm2pay_index
    count=count+1;
    beam_model.Conm2.ID(count+ncom2)     = beam_model.Conm2.ID(iref) + 5000;
    if beam_model.Node.Coord(beam_model.Node.ID == beam_model.Conm2.Node(iref),2) == 0
        beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref);
    else
        beam_model.Conm2.Node(count+ncom2)   = beam_model.Conm2.Node(iref) + 5000;
    end
    beam_model.Conm2.CID(count+ncom2)    = beam_model.Conm2.CID(iref);
    beam_model.Conm2.M(:,:,count+ncom2)  = beam_model.Conm2.M(:,:,iref);
    beam_model.Conm2.Offset(count+ncom2,:) = beam_model.Conm2.Offset(iref,:)*diag([1,-1,1]);
    beam_model.PartIDs.Payload.Conm2IDs = [beam_model.PartIDs.Payload.Conm2IDs,beam_model.Conm2.ID(count+ncom2)];
end
beam_model.Info.ncom2 = ncom2+count;


% Double the surface area
beam_model.Aero.ref.S_ref = 2*beam_model.Aero.ref.S_ref;

% Remove the symmetry condition
beam_model.Aero.state.SIMXZ = 0;

PlotRes = 0;
if PlotRes == 1
    
    tic
    [beam_model.Aero.lattice, beam_model.Aero.ref] = ...
        vlm_setup(1, beam_model.Aero.geo, beam_model.Aero.state, beam_model.Aero.ref);
    toc
    
    figure(444);hold on;
    
    tic
    plot_vlm(beam_model.Aero.lattice,1);
    toc
    
    %         tic
    %         fprintf(fid,'\n\tPlotting Airfoils\n');
    %         title(['Aspect Ratio = ' num2str(beam_model.Aero.Wing.ref.AR)]);
    %         plotairfoils(beam_model.Aero.geo,[beam_model.Optim.Wing.t_skin;beam_model.Optim.Wing.t_skin]);
    %         toc
    %
    %         fprintf('\n\tPlotting 3D box\n')
    %         tic
    %         plot3Dwingbox(beam_model.Optim.Wing.c_box, beam_model.Optim.Wing.h_spar, ...
    %             beam_model.Optim.Wing.t_skin,beam_model.Optim.Wing.t_spar, beam_model.Optim.Wing.nodex,...
    %             -beam_model.Optim.Wing.nodey,beam_model.Optim.Wing.nodez,[])
    %         toc
    %
    
    tic
    plot_conm(beam_model);
    toc
    
end

end