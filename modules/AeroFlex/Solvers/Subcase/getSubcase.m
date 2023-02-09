function beam_model = getSubcase(beam_model,subcase)

beam_model.Param.SPC  = beam_model.Param.SUBCASE{subcase}.SPC;
beam_model.Param.LOAD = beam_model.Param.SUBCASE{subcase}.LOAD;

%beam_model.Param.StructDamp = 0;
for idamp = 1:length(beam_model.Damping.ID)
    if beam_model.Damping.ID(idamp) == beam_model.Param.LOAD
        % Convert to percentage
        beam_model.Param.StructDamp = beam_model.Damping.G(idamp)*100;
        fprintf('\nStructural Damping set to %3.2f%%\n',beam_model.Param.StructDamp);
    end
end
beam_model.Param.GRAV = [];
GRAV = zeros(1,3);
for igrav = 1:length(beam_model.Grav.ID)
    if beam_model.Grav.ID(igrav) == beam_model.Param.LOAD
        % Convert to percentage
        GRAV = beam_model.Grav.Scale(igrav)*beam_model.Grav.Orient(igrav,:);
        fprintf('\nGravitational Vector set to [%2.2f, %2.2f, %2.2f] \n',GRAV(1),GRAV(2),GRAV(3));
        beam_model.Param.GRAV = GRAV;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%% Reset DOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam_model = resetDOF(beam_model);

%%%%%%%%%%%%%%%%%%%%%% Boundary Condtions (SPCD) %%%%%%%%%%%%%%%%%%%%%%%%%%
beam_model.Param.SPCD = [];
count = 0;
for ispcd = 1:beam_model.Info.nspcd
    if beam_model.SPCD.ID(ispcd) == beam_model.Param.SPC
        count = count+1;
        beam_model.Param.SPCD(count) = ispcd;
        fprintf('\nSPCD # %i\n',beam_model.SPCD.ID(ispcd));
    end
end

%%%%%%%%%%%%%%%%%%%%%% Extract Gust Case %%%%%%%%%%%%%%%%%%%%%%%%%%
if beam_model.SolParam.NGusts == 0
    beam_model.SolParam.GustID = [];
else
    beam_model.SolParam.GustID = subcase;
end

%%%%%%%%%%%%%%%%%%%%%% Extract Wing Set Angle %%%%%%%%%%%%%%%%%%%%%%%%%%
if length(beam_model.SolParam.WingAngle) > 1
    beam_model.SolParam.SetAngle = beam_model.SolParam.WingAngle(subcase);  
else
    beam_model.SolParam.SetAngle = beam_model.SolParam.WingAngle(1);
end
%beam_model.SolParam.SetAngle = 5;
end