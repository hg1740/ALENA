function data = PrintMasses(data)

% Fuselage
if isfield(data.PartIDs,'Fuselage')
    Mass_Fuselage = 0;
    if isfield(data.PartIDs.Fuselage,'CBar')
        
        for i = 1:length(data.PartIDs.Fuselage.CBar)
            bar_index = find(data.Bar.ID == data.PartIDs.Fuselage.CBar(i));
            Mass_Fuselage = Mass_Fuselage + data.Bar.M(1,1,1,bar_index)+data.Bar.M(1,1,2,bar_index)+data.Bar.M(1,1,3,bar_index);
        end
    end
    data.WeightBreakdown.Fuselage = Mass_Fuselage;
end
% Wing
if isfield(data.PartIDs,'Wing')
    Mass_Wing = 0;
    if isfield(data.PartIDs.Wing,'CBar')
        for i = 1:length(data.PartIDs.Wing.CBar)
            bar_index = find(data.Bar.ID == data.PartIDs.Wing.CBar(i));
            Mass_Wing = Mass_Wing + data.Bar.M(1,1,1,bar_index)+data.Bar.M(1,1,2,bar_index)+data.Bar.M(1,1,3,bar_index);
        end
    end
    if isfield(data.PartIDs.Wing,'CBeam')
        for i = 1:length(data.PartIDs.Wing.CBeam)
            beam_index = find(data.Beam.ID == data.PartIDs.Wing.CBeam(i));
            Mass_Wing = Mass_Wing + data.Beam.M(1,1,1,beam_index)+data.Beam.M(1,1,2,beam_index)+data.Beam.M(1,1,3,beam_index);
        end
    end
    data.WeightBreakdown.Wing = Mass_Wing;
end
if isfield(data.PartIDs,'HTP')
    % HTP
    Mass_HTP = 0;
    if isfield(data.PartIDs.HTP,'CBar')
        for i = 1:length(data.PartIDs.HTP.CBar)
            bar_index = find(data.Bar.ID == data.PartIDs.HTP.CBar(i));
            Mass_HTP = Mass_HTP + data.Bar.M(1,1,1,bar_index)+data.Bar.M(1,1,2,bar_index)+data.Bar.M(1,1,3,bar_index);
        end
    end
    if isfield(data.PartIDs.HTP,'CBeam')
        for i = 1:length(data.PartIDs.HTP.CBeam)
            beam_index = find(data.Beam.ID == data.PartIDs.HTP.CBeam(i));
            Mass_HTP = Mass_HTP + data.Beam.M(1,1,1,beam_index)+data.Beam.M(1,1,2,beam_index)+data.Beam.M(1,1,3,beam_index);
        end
    end
    data.WeightBreakdown.HTP = Mass_HTP;
end

end