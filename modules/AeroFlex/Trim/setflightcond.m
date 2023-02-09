function beam_model = setflightcond(beam_model,TRIM_INDEX)
%   - C.Szczyglowski (08/10/2019) This function appears to set the
%   different trim and state variables for the problem based on the load
%   cases that have been defined.

fid = beam_model.Param.FID;
index = find(beam_model.Aero.Trim.Select(TRIM_INDEX) == beam_model.Aero.Trim.ID);
dotpos = strfind(beam_model.Param.FILE,'.');
outf = [beam_model.Param.FILE(1:dotpos-1), '_vlm_', num2str(beam_model.Aero.Trim.ID(index)),'.txt'];
beam_model.Aero.lattice = beam_model.Aero.lattice_vlm;

if ~isempty(index)
    
    beam_model.Res = [];
    beam_model.Res.SOL = 'Static linear aerodynamic';
    beam_model.Res.Aero = [];
    
    % set flight mechanics trim params
    beam_model.Aero.Trim.FM = get_free_body_trim_params(beam_model.Aero.Trim.NC(index), ...
        beam_model.Aero.Trim.Param(index).data, beam_model.Aero.Trim.Value(index).data);
    
    [beam_model.Aero.state.alpha, beam_model.Aero.state.betha, beam_model.Aero.state.P, beam_model.Aero.state.Q, beam_model.Aero.state.R] = ...
        get_state_trim_vars(beam_model.Aero.Trim.FM);
    
    beam_model.Aero.state.alpha = D2R(beam_model.Aero.state.alpha);
    beam_model.Aero.state.betha = D2R(beam_model.Aero.state.betha);
    
    % make sure angles are converted
    beam_model.Aero.Trim.FM.Value(2) = D2R(beam_model.Aero.Trim.FM.Value(2));
    beam_model.Aero.Trim.FM.Value(3) = D2R(beam_model.Aero.Trim.FM.Value(3));
    
    % set Tornado state struct
    beam_model.Aero.state.ALT = beam_model.Aero.Trim.ALT(index);
    
    [beam_model.Aero.state.rho, beam_model.Aero.state.p, ...
        beam_model.Aero.state.T, beam_model.Aero.state.a, ...
        beam_model.Aero.state.mu] = ISA_h(beam_model.Aero.state.ALT);
    
    beam_model.Aero.state.AS = beam_model.Aero.Trim.Mach(index) * beam_model.Aero.state.a;
    beam_model.Aero.state.Mach(1) = beam_model.Aero.Trim.Mach(index);
    
    %    convert angular velocities to dimensional values
    beam_model.Aero.state.P = 2*beam_model.Aero.state.P*beam_model.Aero.state.AS/beam_model.Aero.ref.b_ref;
    beam_model.Aero.state.Q = 2*beam_model.Aero.state.Q*beam_model.Aero.state.AS/beam_model.Aero.ref.C_mgc;
    beam_model.Aero.state.R = 2*beam_model.Aero.state.R*beam_model.Aero.state.AS/beam_model.Aero.ref.b_ref;
    %
    if (beam_model.Param.GRDPNT==0)
        beam_model.Aero.geo.ref_point = zeros(1,3);
        beam_model.Aero.geo.CG = zeros(1,3);
        %         beam_model.Aero.geo.ref_point = beam_model.WB.CG;
        %         beam_model.Aero.geo.CG        = beam_model.WB.CG;
    else
        beam_model.Aero.geo.ref_point = beam_model.Node.Coord(beam_model.Param.GRDPNT,:);
        beam_model.Aero.geo.CG = beam_model.Node.Coord(beam_model.Param.GRDPNT,:);
    end
else
    error('Unable to find the required TRIM set %d.', TRIM_INDEX);
end

if beam_model.Aero.geo.nc
    beam_model.Aero.Trim.CS = get_control_surf_trim_params(beam_model.Aero.geo.nc, beam_model.Aero.Trim.NC(index), ...
        beam_model.Aero.Trim.Param(index).data, beam_model.Aero.Trim.Value(index).data, beam_model.Aero.lattice_vlm.Control);
    
    % check trim variables
    ncs = sum(beam_model.Aero.Trim.CS.Fixed);
    nfm = sum(beam_model.Aero.Trim.FM.Fixed);
    
    if (ncs + nfm) ~= beam_model.Aero.Trim.NC(index)
        
        error('Unable to fix trim variables. Wrong variable name given in TRIM card %d.', ...
            beam_model.Aero.Trim.ID(index));
    end
    
end

% check for duplicated AELINK labels
if (beam_model.Info.nlink)
    
    [labels, i] = unique(beam_model.Aero.Trim.Link.ID);
    if (length(labels) ~= beam_model.Info.nlink)
        
        n = [1 : beam_model.Info.nlink];
        dof = beam_model.Aero.Trim.Link.ID(setdiff(n, i));
        
        for k=1:length(dof)
            fprintf(fid, '\n\tWarning: duplicated labels for AELINK card: %d.', ...
                beam_model.Aero.Trim.Link.ID(dof(k)));
        end
        
        error('AELINK entries have duplicated labels.');
        
    end
end

beam_model.Aero.Trim.CS.MPC = [];   % store control surfaces contraint equations
beam_model.Aero.Trim.CS.Coeff = []; % store control surfaces contraint coefficients

% set output struct
beam_model.Res = [];
beam_model.Res.SOL = 'Static linear non-linear rigid trim';
beam_model.Res.FM.Value = repmat(beam_model.Aero.Trim.FM.Value,1,1); % current flight mechanics solution
beam_model.Res.CS.Value = beam_model.Aero.Trim.CS.Value; % current control surfaces solution
%%
beam_model.Res.CS.Fixed = zeros(1,beam_model.Aero.geo.nc); % set later

if beam_model.Aero.geo.nc
    [beam_model.Aero.Trim.CS.MPC, beam_model.Aero.Trim.CS.Coeff, beam_model.Res.CS.Fixed] = ...
        set_constr_eq(beam_model.Aero.geo.nc, beam_model.Aero.lattice.Control, beam_model.Aero.Trim.CS.Fixed, beam_model.Info.nlink, beam_model.Aero.Trim.Link);
end

%***********************************************************************************************************************
% set output struct
beam_model.Res.FM.Fixed = beam_model.Aero.Trim.FM.Fixed;
beam_model.Res.state = beam_model.Aero.state; % Tornado dummy state (changed during nonlinear solution search)
% set constraints for control surfaces
beam_model.Aero.lattice_defo=beam_model.Aero.lattice;

state = beam_model.Res.state;
% fprintf(fid, '\n\n - Aerodynamic state:');
% fprintf(fid, '\n\tFlight speed [m/s]: %g.', state.AS);
% fprintf(fid, '\n\tMach number: %g.', state.Mach(1));
% fprintf(fid, '\n\tAltitude [m]: %g.', state.ALT);
% fprintf(fid, '\n\tDensity [kg/m^3]: %g.', state.rho);

beam_model.Aero.lattice=beam_model.Aero.lattice_defo;

end
