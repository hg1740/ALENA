function beam_model = aeroelastic_interface(beam_model)

fid = 1;

if beam_model.Info.spline_type == 1
    SPLINE_TYPE = 'Linear beam spline';
    SPLINE_ID = 1;
elseif beam_model.Info.spline_type == 2
    SPLINE_TYPE = 'Surface spline';
    SPLINE_ID = 2;
else
    SPLINE_TYPE = 'UNKNOWN';
    SPLINE_ID = inf;
end

%fprintf(fid, '\n\tBuilding FS interpolation matrices...');
if beam_model.Info.ninterp == 0
    error('No interpolation set for aeroelastic interface given.');
end
%fprintf(fid, '\n\t - Method: %s', SPLINE_TYPE);

% VLM interface
switch(SPLINE_ID)
    case {1}
        [beam_model.Aero.Interp.Ic, beam_model.Aero.Interp.In, beam_model.Aero.Interp.Iv, beam_model.Aero.Interp.Imv] = ...
            interf_vlm_model1(fid, beam_model.Info.ncaero, beam_model.Node, beam_model.Aero);
    case {2,3}
        [beam_model.Aero.Interp.Ic, beam_model.Aero.Interp.In, beam_model.Aero.Interp.Iv, beam_model.Aero.Interp.Imv] = ...
            interf_vlm_model(fid, beam_model.Info.ncaero, beam_model.Node, beam_model.Aero,beam_model.MeshType);
        
end

end
