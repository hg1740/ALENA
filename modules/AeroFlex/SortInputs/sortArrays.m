%% sortArrays

% Sort and check for duplicate entries

function deck = sortArrays(deck)

[deck.Param] = ...
    sortPARAMArray(deck.Param,deck.Node);

[deck.Node] = ...
    sortNODEArray(deck.Node,deck.Info);

%[deck.COORD1] = ...
%    sortCOORD1Array(deck.COORD1,deck.Node,deck.Info);

[deck.COORD2] = ...
    sortCOORD2Array(deck.COORD2,deck.Info);

[deck.Coord,deck.Info] = ...
    combineCOORDArrays(deck.Coord,deck.COORD1,deck.COORD2,deck.Info);

[deck.Node] = ...
    RotateNODECOORD(deck.Node,deck.Coord);

[deck.PBar,deck.Info]     = ...
    sortPBarArray(deck.PBar,deck.Mat,deck.Param,deck.Info);

[deck.PBeam,deck.Info]     = ...
    sortPBeamArray(deck.PBeam,deck.Mat,deck.Param,deck.Info);

[deck.Bar,deck.Node,deck.Info] = ...
    sortBarArray(deck.Bar,deck.PBar,deck.Node,deck.Mat,deck.Coord,deck.Info);

[deck.Beam,deck.Node,deck.Info] = ...
    sortBeamArray(deck.Beam,deck.PBeam,deck.Node,deck.Mat,deck.Coord,deck.Info);

[deck.Node, deck.RBE0]    = ...
    sortRBE0Array(deck.Node,deck.RBE0,deck.Info);

 [deck.ConM,deck.Info] = ...
     sortCONMArray(deck.Conm1,deck.Conm2,deck.Node,deck.Coord,deck.Param,deck.Info);

%[deck.ConM,deck.Info] = ...
%    groupCONMArray(deck.Conm2,deck.Node,deck.Coord,deck.Param,deck.Info);

[deck.Node] = ...
    sortNodeIndex(deck.Node,deck.RBE2,deck.Bar,deck.ConM,deck.Beam,deck.RBE0,deck.Celas,deck.Info);

[deck.Aero] = ...
    sortINTERPdata(deck.Aero,deck.Node,deck.Info);

[deck.SPC] = ...
    sortSPCArray(deck.Param,deck.SPC,deck.Node,deck.Info);

[deck.Node,deck.Info]     = ...
    set_NODE_DOF(deck.Param,deck.Node,deck.SPC,deck.Info);

[deck.Celas] = ...
    sortCELASArray(deck.Celas,deck.Node,deck.Info);

%[deck.Mat,deck.Info] = ...
%    sortMatArray(deck.Mat,deck.Info);

[deck.F] = ...
    sortFORCEArray(deck.F,deck.Node,deck.Coord,deck.Info);

[deck.M] = ...
    sortMOMENTArray(deck.M,deck.Node,deck.Coord,deck.Info);

[deck.F_FLW] = ...
    sortFLWArray(deck.F_FLW,deck.Node,deck.Coord,deck.Info);

if deck.Info.nrbe2 >0
    [deck.RBE2,deck.Node,deck.Info.nrbe2,deck.Info.ndof] = ...
        checkRBE2(1,deck.RBE2,deck.Node,deck.Param,deck.Bar,deck.Beam,deck.Celas,deck.Info.ngrid,deck.SPC,deck.Info.ndof);
else
    deck.Node.DOF2 = deck.Node.DOF;
    deck.Node.DOF3 = deck.Node.DOF;
end

[deck.RJoint] = ...
    sortRJOINTArray(deck.RJoint,deck.Node,deck.Info);

[deck.CBush] = ...
    sortCBUSHArray(deck.CBush,deck.PBush,deck.Node,deck.Info);

[deck.RBar] = ...
    sortRBARArray(deck.RBar,deck.Node,deck.Info);

[deck.Aero] = ...
    sortCAEROArray(deck.Aero,deck.Param,deck.Node,deck.Coord,deck.Info);

if ~isempty(deck.BAero.ID)
[deck.Aero,deck.Info] = ...
    sortBAEROArray(deck.BAero,deck.Aero,deck.Coord,deck.Info);
end

if deck.Info.naeros
    [deck.Aero] = ...
        sortAEROSArray(deck.Aero,deck.Aeros);
end

end