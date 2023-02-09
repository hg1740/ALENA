function Sim = initSystem(x0,Matrices,Aero,AnalysisParam)

Sim = struct();

% if ~isfield(Sim,'eps') %Done
%     Sim.eps = 1e-6;
% end
% if ~isfield(Sim,'damping') %Done
%     Sim.damping = 1.00;
% end
% if ~isfield(Sim,'tangent_update') %Done
%     Sim.tangent_update = 5;
% end
if ~isfield(Sim,'rb_flag')
    Sim.rb_flag = 0;
end
if ~isfield(Sim,'grav_vec')
    if ~isfield(Sim,'grav_fact')
        Sim.grav_fact = 0;
    end
    Sim.grav_vec = [0;0;-9.807;0;0;0]*Sim.grav_fact;
else
    if ~isfield(Sim,'grav_fact')
        Sim.grav_fact = 1;
    end
end
if ~isfield(Sim,'rbloads_flag')
    Sim.rbloads_flag = 0;
end

%Analysis dependent parameters
switch AnalysisParam.AnalysisType
    
    case 'static'
        
        %Define indices for system matrices
        Sim.Ind = getSystemInd(Matrices, Sim, [], AnalysisParam);
        
    otherwise
        
        if floor(Sim.Soln) == 2
            Sim.Ind      = getSystemInd(Matrices,Sim,[]);
            if Sim.time.dt == 0
                Sim.time_vec = 1;
            else
                Sim.time_vec = Sim.time.limits(1):Sim.time.dt:Sim.time.limits(end);
            end
            if ~isfield(Sim,'rho_inf')
                Sim.rho_inf = 0.99;
            end
            if ~isfield(Sim,'struct_damp')
                Sim.struct_damp = 0;
            end
            if Sim.Soln == 2.1
                Sim.speed_up = 1;
            elseif ~isfield(Sim,'speed_up')
                Sim.speed_up = 0;
            end
            
            alphm1      = 0.5*(3-Sim.rho_inf)/(1+Sim.rho_inf);
            alphf1      = 1/(1+Sim.rho_inf);
            Sim.gam1    = 0.5 + alphm1 - alphf1;
        end
        
end

end