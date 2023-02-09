function Printbox(beam_model,Optim)

Nsec = Optim.Nsec;

fprintf('\n\tWing properties');
fprintf('\n\t- AR = %d',beam_model.Aero.ref.AR);
fprintf('\n\t- Root chord  = %-.2f m',beam_model.Aero.ref.C_mac);
fprintf('\n\t- Wing Length = %-.2f m',beam_model.Aero.ref.b_ref);
fprintf('\n\t- Wing mass   = %-.2f Kg',Optim.M);
fprintf('\n\t- Wing surface area = %-.2f m^2',beam_model.Aero.ref.S_ref);
fprintf('\n\tOptimisation inital condition');
if Optim.Type == 2
    fprintf('\n\t- Root skin thickness = %-.2f mm',Optim.X0(1)*1000);
    fprintf('\n\t- Tip skin thickness  = %-.2f mm',Optim.X0(Nsec)*1000);
    fprintf('\n\t- Root spar thickness = %-.2f mm',Optim.X0(Nsec + 1)*1000);
    fprintf('\n\t- Tip spar thickness  = %-.2f mm',Optim.X0(2*Nsec)*1000);
elseif Optim.Type == 2
    fprintf('\n\t- Root skin thickness = %-.2f mm',Optim.X0(1)*1000);
    fprintf('\n\t- Tip skin thickness  = %-.2f mm',Optim.X0(Nsec)*1000);
    fprintf('\n\t- Root spar thickness = %-.2f mm',Optim.X0(Nsec + 1)*1000);
    fprintf('\n\t- Tip spar thickness  = %-.2f mm',Optim.X0(2*Nsec)*1000);
    fprintf('\n\t- Root stringer thickness = %-.2f mm^2',Optim.X0(2*Nsec + 1)*1e6);
    fprintf('\n\t- Tip stringer thickness = %-.2f mm^2',Optim.X0(3*Nsec)*1e6);
end

fprintf('\n');

end