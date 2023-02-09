function beam_model = generateVTPAero(beam_model,Bar_prop,Aircraft_param,Parameters,n_loads)

SPLINE1Idx = Aircraft_param.VTP.Planform.SPLINE1Idx;
CAEROIdx   = Aircraft_param.VTP.Planform.CAEROIdx;

if strcmp(Parameters.SplineType,'Surface')
    
    [beam_model.Aero,beam_model.Info,beam_model.PartIDs.VTP.SplineIDs] = ...
        SPLINECreator(beam_model.Aero,beam_model.Aero.VTP,Bar_prop,...
        beam_model.Info,CAEROIdx,SPLINE1Idx);
    
elseif strcmp(Parameters.SplineType,'Beam')
    
    [beam_model.Aero,beam_model.Info,beam_model.PartIDs.VTP.SplineIDs] = ...
        SPLINECreatorLinear(beam_model.Aero,beam_model.Aero.VTP,Bar_prop,...
        beam_model.Info,CAEROIdx,SPLINE1Idx);
end

beam_model = Combine_Aero(beam_model,{'Wing','StbdHTP','VTP'});

beam_model.Info.ncaero = length(beam_model.Aero.ID);

% Add the VTP Control surface to the trim variables
NC = length(beam_model.Aero.Trim.Param(1).data);
for idx = 1:n_loads
    if ~isempty(beam_model.Aero.VTP.Trim.MasterSurf)
        beam_model.Aero.Trim.Param(idx).data(NC+1,1) = beam_model.Aero.VTP.Trim.MasterSurf(1);
        beam_model.Aero.Trim.Value(idx).data(NC+1,1) = 0.0;
    end
end

% Regenerate lattice to include the VTP
[beam_model.Aero.lattice, beam_model.Aero.ref] = ...
    vlm_setup(1, beam_model.Aero.geo, beam_model.Aero.state, beam_model.Aero.ref);

beam_model.Aero.lattice.Control.Name = beam_model.Aero.Control.Name;

beam_model.Aero = rmfield(beam_model.Aero,'Control');

beam_model.Aero.lattice_vlm = beam_model.Aero.lattice;

% Renumber the path and set fields
for n_vtp = 1:beam_model.Info.ninterp
    i_vtp = find(beam_model.Aero.ID == beam_model.Aero.Interp.Patch(n_vtp));
    m_vtp = find(beam_model.Aero.Set.ID == beam_model.Aero.Interp.Set(n_vtp));
    if ~isempty(i_vtp)
        beam_model.Aero.Interp.Patch(n_vtp) = i_vtp; % substitute with patch index for rapid accessing
    end
    if ~isempty(m_vtp)
        beam_model.Aero.Interp.Set(n_vtp) = m_vtp;
    end
end
end