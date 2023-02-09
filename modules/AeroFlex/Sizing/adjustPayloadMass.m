function beam_model = adjustPayloadMass(beam_model,LoadCase,Payload_Idx)

% Reassign the payload properties to the correct conm2 entry
PayloadMass = beam_model.Conm2.M(1:6,1:6,Payload_Idx);

beam_model.Conm2.M(1:6,1:6,Payload_Idx)  = zeros(6);

[beam_model.ConM,beam_model.Info] = ...
    sortCONMArray(beam_model.Conm1,beam_model.Conm2,beam_model.Node,beam_model.Coord,beam_model.Param,beam_model.Info);

beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);

% --------------------------
% Allocate payload properties
% --------------------------

% Calculate the
Aircraft_cg_x = LoadCase.PayloadCG*beam_model.Aero.ref.MAC + beam_model.Aero.ref.MAC_LE_x;
Payload_cg = [(Aircraft_cg_x*(beam_model.WB.MCG(1) + PayloadMass(1)) - beam_model.WB.CG(1)*beam_model.WB.MCG(1))/PayloadMass(1),0,0];

% Recalculate the offset of the lumped mass
NodeOffset      = Payload_cg - beam_model.Node.Coord((beam_model.Node.ID == beam_model.Conm2.Node(Payload_Idx)),:);

% Reassign the payload properties to the correct conm2 entry
beam_model.Conm2.M(1:3,1:3,Payload_Idx)  = eye(3) * PayloadMass(1,1);
beam_model.Conm2.M(4:6,4:6,Payload_Idx)  = PayloadMass([4:6],[4:6]);

beam_model.Conm2.Offset(Payload_Idx,1:3) = NodeOffset;

% --------------------------
% Recalculate model mass properties
% --------------------------
[beam_model.ConM,beam_model.Info] = ...
    sortCONMArray(beam_model.Conm1,beam_model.Conm2,beam_model.Node,beam_model.Coord,beam_model.Param,beam_model.Info);

beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);

beam_model.Res.WB.CG = beam_model.WB.CG;

end