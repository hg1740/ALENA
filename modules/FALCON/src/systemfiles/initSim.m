function x_out = initSim(x_in, Matrices, AnalysisParam)
%initSim Defines the initial conditions for the simulation. Different
%initial conditions are defined based on the solution type (i.e. static vs.
%dynamic).
%
% Detailed Description:
%   * Static simulation - In the static simulation the
%

%How many elements in each component
nElem = Matrices.n_elem;

%How many components?
nPart = numel(Matrices.n_elem);

switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %Local loads (forces and moments) for each element?
        if ~isfield(x_in, 'x_f')
            x_in.x_f = arrayfun(@(nE) zeros(nE * 6, 1), nElem, 'Unif', false);
            %         for jj = 1:nPart
            %             x_in.x_f{jj} = zeros(Matrices.n_elem(jj)*6,1);
            %         end
        end
        
        if isfield(x_in, 'x_vg') && ~isfield(x_in, 'x_va')
            %Transform global (freestream) velocity into aircraft system
            x_in.x_va = [Matrices.CGa' zeros(3); zeros(3) Matrices.CGa'] * x_in.x_vg;
        end
        
        if ~isfield(x_in, 'x_va')
            x_in.x_va = zeros(6,1);
        end        
        
        %Stash a copy of inputs
        x_out.x_f  = x_in.x_f;
        x_out.x_va = x_in.x_va;
        
    otherwise
        
        getV_flag = 0;
        if ~isfield(x_in,'x_f')
            for jj = 1:nPart
                x_in.x_f{jj} = zeros(Matrices.n_elem(jj)*6,1);
            end
        end
        if ~isfield(x_in,'x_v')
            getV_flag = 1;
            for jj = 1:nPart
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
        
        % Determine all initial inputs from initial conditions
        for jj = 1:nPart
            x_out.x_f{jj} = x_in.x_f{jj};
            x_out.x_v{jj} = x_in.x_v{jj};
            x_out.x_q{jj} = zeros(Matrices.n_elem(jj)*4,1);
            x_out.x_p{jj} = zeros(Matrices.n_elem(jj)*3,1);
            
            % Back-calculate orientation and positions from loads
            [coords ,CaB ] = strains2coords_all(x_out.x_f,Matrices,jj);
            for ii = 1:length(coords)
                coords_min{jj}(:,ii) = coords(:,:,ii);
            end
            for ii = 1:size(Matrices.s_mp_ind{jj},2)
                ind  = [1:4] + (ii-1)*4;
                ind2 = [1:3] + (ii-1)*3;
                ind3 = [1:6] + (ii-1)*6;
                x_out.x_q{jj}(ind)       = Rot2Quat(Matrices.CGa*CaB(:,:,Matrices.s_mp_ind{jj}(ii)));
                x_out.x_p{jj}(ind2)      = Matrices.CGa*(coords(:,:,Matrices.s_mp_ind{jj}(ii)) + Matrices.p0{jj});
                if getV_flag
                    CBa = [CaB(:,:,Matrices.s_mp_ind{jj}(ii))' zeros(3);zeros(3) CaB(:,:,Matrices.s_mp_ind{jj}(ii))']*[eye(3) -skew(coords(:,:,Matrices.s_mp_ind{jj}(ii)));zeros(3) eye(3)];
                    x_out.x_v{jj}(ind3) = CBa*x_in.x_va(:,1);
                end
            end
        end
        x_out.x_va = zeros(6,1);
        x_out.x_qa = zeros(4,1);
        x_out.x_pa = zeros(3,1);
        
        % Set initial velocity to the global velocity rotated into the local frame
        x_out.x_va = x_in.x_va;%Matrices.CGa'*V;
        
        % Set initial orientation in quaternions from the initial rotation matrix
        x_out.x_qa   = Rot2Quat(Matrices.CGa);
end

end