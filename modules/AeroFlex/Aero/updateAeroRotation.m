function beam_model = updateAeroRotation(beam_model)

if (beam_model.Info.nrbe0 > 0)
    for n = 1:beam_model.Info.ngrid
        ne = length(beam_model.Node.Aero.Index(n).data);
        if ne
            beam_model.Res.NRd(:,:,beam_model.Node.Aero.Index(n).data) = repmat(beam_model.Res.NRd(:,:,n),1,1,length(beam_model.Node.Aero.Index(n).data));
        end
    end
end

NR = permute(beam_model.Res.NRd,[1,3,2]);

AeroRd(:,1,:) = (beam_model.Aero.Interp.Ic*NR(:,:,1)')';
AeroRd(:,2,:) = (beam_model.Aero.Interp.Ic*NR(:,:,2)')';
AeroRd(:,3,:) = (beam_model.Aero.Interp.Ic*NR(:,:,3)')';

for n = 1:size(AeroRd,3)
    beam_model.Res.AeroRd(:,:,n) = AeroRd(:,:,n)*beam_model.Aero.AeroRd(:,:,n);
end

%beam_model.Res.AeroRd = AeroRd;

end