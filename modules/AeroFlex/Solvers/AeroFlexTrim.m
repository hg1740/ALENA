%AEROFLEXTRIM:
%   Summary: The following script takes a beam model structure and runs a
%   linear/nonlinear trim analysis. The following code should not be
%   changed.

function Static = AeroFlexTrim(beam_model,linear)

% Initialise the updated lattice with the undeformed lattice
beam_model.Aero.lattice = beam_model.Aero.lattice_vlm;

% Generate the the FSI interpolation matrices
beam_model = aeroelastic_interface(beam_model);

% Setup some useful matrices
beam_model = set_state_matrices(beam_model);

% Initialise the Satic variable as empty (Results are stored here)
Static = [];

% Run through the various trim selections (This should be parallelised)
for iload = beam_model.SolParam.trimI
    
    % Reset the rigid body mass matrix and c.g. locations
    beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);
    
    % Setup the flight conditions for this trim selection
    beam_model                  = setflightcond(beam_model,iload);
    
    % Reset a few variables between trims
    beam_model.Res.WB.CG        = beam_model.WB.CG;
    beam_model.Aero.lattice     = beam_model.Aero.lattice_vlm;
    beam_model.Aero.AIC         = [];
       
    if linear == 0
        beam_model.Res.NRd  = beam_model.Node.R;
        
        AeroRot = [];
        for n = 1:size(beam_model.Aero.lattice_vlm.N,1)
            ez = -beam_model.Aero.lattice_vlm.N(n,:)';
            ex = [1,0,0]';
            ey = cross(ez,ex);
            AeroRot(:,:,n) = [ex,ey,ez];
        end
        beam_model.Aero.AeroRd          = AeroRot;
        beam_model.Res.AeroRd           = AeroRot;
        
        if beam_model.MeshType == 1 || beam_model.MeshType == 2
            beam_model                      = StripSlopes(beam_model);
            beam_model.Aero.geo.panelarea   = tarea(beam_model.Aero.lattice_vlm.XYZ);
        end
         
        %beam_model.Aero.Cl_alpha(:) = 2*pi;

        beam_model          = updateAeroRotation(beam_model);
        
        beam_model.SolParam.SetAngle = 0;
        beam_model.Res.state.alpha = beam_model.SolParam.SetAngle;
        
        beam_model.Res.FM.Value(2) = beam_model.Res.state.alpha;
        
        if beam_model.SolParam.Trim == 1
            beam_model.Res = rigid_trim(beam_model,beam_model.Res,beam_model.Aero);
        else
            beam_model.Res = rigid_aero(beam_model,beam_model.Res);
        end
        
    end
    
    % Call the linear/nonlinear aeroelastic trim solver
    if linear == 1
        [Static,beam_model] = solve_linear_trim_customisedV2(Static,beam_model,iload);
    else
        [Static,beam_model] = solve_nonlinear_trim_customised(Static,beam_model,iload);
    end
  
end

end