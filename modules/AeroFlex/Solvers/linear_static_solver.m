
function [NODEPOS,beam_model,BarForces] = linear_static_solver(F,beam_model,Kll,ldof,ENABLE_DISP)

%%%%%%%%%%%%%%%%%%%%%%%%%%% COUNTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ngrid = beam_model.Info.ngrid;
nbar  = beam_model.Info.nbar;
nbeam = beam_model.Info.nbeam;
ndof  = beam_model.Info.ndof;
ndof2 = beam_model.Info.ndof2;

%beam_model.WingNodes = SortWingNodes(beam_model);
WngNodes = beam_model.WingNodes;
ndof =  beam_model.Info.ndof;

if ~isempty(beam_model.RBE2)
    FRBE2 = RBE2Assembly2(beam_model.RBE2,F);
else
    FRBE2 = F;
end

dummyndof = size(FRBE2,1);
SOL = zeros(dummyndof, 1);
SOL(ldof) = Kll \ FRBE2(ldof,:);

if ~isempty(beam_model.RBE2)
    SOL = RBE2disp(beam_model.RBE2,SOL,ndof);
end

gdef = zeros(beam_model.Info.ngrid, 6);

for n = 1:beam_model.Info.ngrid
    
    dof = beam_model.Node.DOF(n, 1:6);
    index = find(dof);
    
    if ~isempty(index)
        gdef(n, index) = SOL(dof(index));
    end
end

% update interal database
beam_model.Res.SOL = 'Static linear';

% store nodal displacement
beam_model.Res.NDispl = gdef;
beam_model.Res.NRd = zeros(3, 3, beam_model.Info.ngrid);

% set delta Rot
for n = 1:beam_model.Info.ngrid
    beam_model.Res.NRd(:,:,n) = Rmat(gdef(n, 4:6)');
end

if (beam_model.Info.nrbe0 > 0)
    AERO_POS = update_aerobeam_node(beam_model.Info.ngrid, beam_model.Node, beam_model.Res.NDispl(:,1:3,1),...
        beam_model.Res.NRd-repmat(eye(3,3),[1,1,beam_model.Info.ngrid]));
    
    % update coord database with slave nodes position
    for n=1:beam_model.Info.ngrid
        ne = length(beam_model.Node.Aero.Index(n).data);
        if ne
            beam_model.Res.NDispl(beam_model.Node.Aero.Index(n).data, 1:3, 1) = AERO_POS(n).data';
        end
    end
    clear AERO_POS;
end

% store bar internal forces
beam_model.Res.Bar.CForces = [];
beam_model.Res.Beam.CForces = [];

% store bar strains and curvatures
beam_model.Res.Bar.CStrains = [];
beam_model.Res.Beam.CStrains = [];

% store bar stresses
beam_model.Res.Bar.CStresses = [];
beam_model.Res.Bar.CSM = [];
beam_model.Res.Beam.CStresses = [];
beam_model.Res.Beam.CSM = [];

% assembly BAR contributions directly in the undeformed position
[beam_model.Res.Bar.CForces, beam_model.Res.Bar.CStrains, beam_model.Res.Bar.CStresses, beam_model.Res.Bar.CSM] = ...
    get_bar_force_strain(beam_model.Info.nbar, beam_model.Bar, beam_model.PBar, beam_model.Mat, beam_model.Node, ...
    beam_model.Res.NDispl, beam_model.Param.FUSE_DP);

[beam_model.Res.Beam.CForces, beam_model.Res.Beam.CStrains, beam_model.Res.Beam.CStresses, beam_model.Res.Beam.CSM] = ...
    get_bar_force_strain(beam_model.Info.nbeam, beam_model.Beam, beam_model.PBeam, beam_model.Mat, beam_model.Node, ...
    beam_model.Res.NDispl, beam_model.Param.FUSE_DP);

NODEPOS = beam_model.Node.Coord + beam_model.Res.NDispl(:, 1:3);

BARR   = update_bar_rot(nbar, beam_model.Res.NRd, beam_model.Bar.Conn, beam_model.Bar.R, gdef);
BEAMR  = update_bar_rot(nbeam, beam_model.Res.NRd, beam_model.Beam.Conn, beam_model.Beam.R, gdef);

% Storing the internal loads at the midpoint of each of the beams
Fx = [];
Fy = [];
Fz = [];
Mx = [];
My = [];
Mz = [];

if ~isempty(beam_model.Res.Bar.CForces)
    for j = 1:size(beam_model.Bar.Colloc,3)
        Fx    = [Fx;mean([beam_model.Res.Bar.CForces(1,1,j),beam_model.Res.Bar.CForces(2,1,j)])];
        Fy    = [Fy;mean([beam_model.Res.Bar.CForces(1,2,j),beam_model.Res.Bar.CForces(2,2,j)])];
        Fz    = [Fz;mean([beam_model.Res.Bar.CForces(1,3,j),beam_model.Res.Bar.CForces(2,3,j)])];
        Mx    = [Mx;mean([beam_model.Res.Bar.CForces(1,4,j),beam_model.Res.Bar.CForces(2,4,j)])];
        My    = [My;mean([beam_model.Res.Bar.CForces(1,5,j),beam_model.Res.Bar.CForces(2,5,j)])];
        Mz    = [Mz;mean([beam_model.Res.Bar.CForces(1,6,j),beam_model.Res.Bar.CForces(2,6,j)])];
    end
end

if ~isempty(beam_model.Res.Beam.CForces)
    for j = 1:size(beam_model.Beam.Colloc,3)
        Fx    = [Fx;mean([beam_model.Res.Beam.CForces(1,1,j),beam_model.Res.Beam.CForces(2,1,j)])];
        Fy    = [Fy;mean([beam_model.Res.Beam.CForces(1,2,j),beam_model.Res.Beam.CForces(2,2,j)])];
        Fz    = [Fz;mean([beam_model.Res.Beam.CForces(1,3,j),beam_model.Res.Beam.CForces(2,3,j)])];
        Mx    = [Mx;mean([beam_model.Res.Beam.CForces(1,4,j),beam_model.Res.Beam.CForces(2,4,j)])];
        My    = [My;mean([beam_model.Res.Beam.CForces(1,5,j),beam_model.Res.Beam.CForces(2,5,j)])];
        Mz    = [Mz;mean([beam_model.Res.Beam.CForces(1,6,j),beam_model.Res.Beam.CForces(2,6,j)])];
    end
end

BarForces.Fx = Fx;
BarForces.Fy = Fy;
BarForces.Fz = Fz;
BarForces.Mx = Mx;
BarForces.My = My;
BarForces.Mz = Mz;

beam_model.Res.Bar.Colloc = bar_defo_colloc(beam_model.Info.nbar, beam_model.Bar, beam_model.Node.DOF, NODEPOS, beam_model.Res.NRd);
beam_model.Res.Beam.Colloc = bar_defo_colloc(beam_model.Info.nbeam, beam_model.Beam, beam_model.Node.DOF, NODEPOS, beam_model.Res.NRd);

if ENABLE_DISP==1
    figure;
    hold on;
    plot(beam_model.Node.Coord(WngNodes,2),beam_model.Res.NDispl(WngNodes,1),'b^--');
    hold on;
    plot(beam_model.Node.Coord(WngNodes,2),beam_model.Res.NDispl(WngNodes,2),'bo--');
    hold on;
    plot(beam_model.Node.Coord(WngNodes,2),beam_model.Res.NDispl(WngNodes,3),'b+--');    
    legend('x disp','y disp','z disp','Location','Best');
    
    for i = 1:beam_model.Info.nbar
        hold on;
        plot(beam_model.Res.Bar.Colloc(:,2,i),beam_model.Res.Bar.Colloc(:,3,i),'go');
        hold on;
    end
else
end

beam_model.Res.Bar.R  = BARR;
beam_model.Res.Beam.R = BEAMR;

end
