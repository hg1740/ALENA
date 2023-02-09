function x_dyn = initDynSim(x_in,Matrices,Sim)

getV_flag = 0;
if ~isfield(x_in,'x_f')
    for jj = 1:length(Matrices.n_elem)
        x_in.x_f{jj} = zeros(Matrices.n_elem(jj)*6,1);
    end
end
if ~isfield(x_in,'x_v')
    getV_flag = 1;
    for jj = 1:length(Matrices.n_elem)
        x_in.x_v{jj} = zeros(Matrices.n_elem(jj)*6,1);
    end
end
if isfield(x_in,'x_vg') && ~isfield(x_in,'x_va')
    x_in.x_va = [Matrices.CGa' zeros(3); zeros(3) Matrices.CGa']*x_in.x_vg;
end
if ~isfield(x_in,'x_va')
    if ~Sim.rb_flag%isfield(Sim,'rb_flag')
        x_in.x_va = zeros(6,length(Sim.time.limits(1):Sim.time.dt:Sim.time.limits(2)));
    else
        x_in.x_va = zeros(6,1);
    end
else
    if ~Sim.rb_flag && size(x_in.x_va,2) == 1
        x_in.x_va = repmat(x_in.x_va,1,length(Sim.time.limits(1):Sim.time.dt:Sim.time.limits(2)));
    end
end

for jj = 1:length(Matrices.n_elem)
    % Back-calculate orientation and positions from loads
    [coords{jj} ,CaB{jj} ] = strains2coords_all(x_in.x_f,Matrices,jj);
end
% Determine all initial inputs from initial conditions
for jj = 1:length(Matrices.n_elem)
    x_dyn.x_f{jj} = x_in.x_f{jj};
    x_dyn.x_v{jj} = x_in.x_v{jj};
    x_dyn.x_q{jj} = zeros(Matrices.n_elem(jj)*4,1);
    x_dyn.x_p{jj} = zeros(Matrices.n_elem(jj)*3,1);
    
    if Matrices.Parent{jj}(1)
        coords_parent = coords{Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)));
        CaB_parent    = CaB{   Matrices.Parent{jj}(1)}(:,:,Matrices.s_node_ind{Matrices.Parent{jj}(1)}(Matrices.Parent{jj}(2)))*Matrices.CaB0{Matrices.Parent{jj}(1)}';
        for ii = 1:size(Matrices.s_out{jj},2)
            coords{jj}(:,:,ii) = CaB_parent*coords{jj}(:,:,ii) + coords_parent;
            CaB{jj}(:,:,ii)    = CaB_parent*CaB{jj}(:,:,ii);
        end
    end
    
%     [coords ,CaB ] = strains2coords_all(x_dyn.x_f,Matrices,jj);
%     for ii = 1:length(coords)
%         coords_out{jj}(:,ii) = coords(:,:,ii);
%     end
    for ii = 1:size(Matrices.s_mp_ind{jj},2)
        ind  = [1:4] + (ii-1)*4;
        ind2 = [1:3] + (ii-1)*3;
        ind3 = [1:6] + (ii-1)*6;
        x_dyn.x_q{jj}(ind)       = Rot2Quat(Matrices.CGa*CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii)));
%         x_dyn.x_q{jj}(ind)       = x_dyn.x_q{jj}(ind)/norm(x_dyn.x_q{jj}(ind));
        x_dyn.x_p{jj}(ind2)      = Matrices.CGa*(coords{jj}(:,:,Matrices.s_mp_ind{jj}(ii)) + Matrices.p0{jj});
        if getV_flag
            CBa = [CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii))' zeros(3);zeros(3) CaB{jj}(:,:,Matrices.s_mp_ind{jj}(ii))']*[eye(3) -skew(coords{jj}(:,:,Matrices.s_mp_ind{jj}(ii)));zeros(3) eye(3)];
            x_dyn.x_v{jj}(ind3) = CBa*x_in.x_va(:,1);
        end
    end
end
x_dyn.x_va = zeros(6,1);
x_dyn.x_qa = zeros(4,1);
x_dyn.x_pa = zeros(3,1);

% Set initial velocity to the global velocity rotated into the local frame
x_dyn.x_va = x_in.x_va;%Matrices.CGa'*V;

% Set initial orientation in quaternions from the initial rotation matrix
x_dyn.x_qa   = Rot2Quat(Matrices.CGa);

if Sim.Soln == 2.2
    x_dyn.x = [x_dyn.x_v{1};x_dyn.x_f{1}];
end