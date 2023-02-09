function beam_model = resetDOF(beam_model)

% Sort out the input decks ready for analysis
[beam_model.Node,beam_model.Info]     = ...
    set_NODE_DOF(beam_model.Param,beam_model.Node,beam_model.SPC,beam_model.Info);

if beam_model.Info.nrbe2 >0
    [beam_model.RBE2,beam_model.Node,beam_model.Info.nrbe2,beam_model.Info.ndof] = ...
        checkRBE2(1,beam_model.RBE2,beam_model.Node,beam_model.Param,beam_model.Bar,beam_model.Beam,beam_model.Celas,beam_model.Info.ngrid,beam_model.SPC,beam_model.Info.ndof);
    beam_model.Info.ndof2 = max(max(beam_model.Node.DOF2));
else
    beam_model.Node.DOF2 = beam_model.Node.DOF;
    beam_model.Node.DOF3 = beam_model.Node.DOF;
    beam_model.Info.ndof2 = max(max(beam_model.Node.DOF2));
end

end