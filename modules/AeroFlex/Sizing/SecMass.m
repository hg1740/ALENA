
function [SecondaryMass,beam_model] = SecMass(beam_model,Aircraft_param,OPTIM)

% The following formulae are taken from Torenbeek (1992)
% ------------------------------------
%
% Inputs:   - beam_model (structure)
%               - beam_model.Breguet.MTOW
%               - beam_model.PartIDs.Wing
%               - beam_model.PartIDs.HTP
%               - OPTIM.Wing.V_int
%               - OPTIM.HTP.V_int
%           - Aircraft_param (structure)
%               - Aircraft_param.Wing.Planform.fspar
%               - Aircraft_param.Wing.Planform.aspar
%               - Aircraft_param.Wing.AeroSurface.Swet
%               - Aircraft_param.Wing.AeroSurface.Swet
%
% Outputs:  - SecondaryMass (structure)
%           - beam_model (updated)
%
% ------------------------------------

% ------------------------------------
% Leading Edge sections from Torenbeek
% ------------------------------------
% Kfle = 1.0 no le devices, 1.4 with le devices N/m^2
Kfle = 1.4;
Wfle_sfle = 75*Kfle * (1+ sqrt(Aircraft_param.Weights.MTOW*9.81*10^(-6)));
% Multiply 2 by as the wetted surface area is for the half model
sfle = 2*mean(OPTIM.Wing.f_spar)*Aircraft_param.Wing.AeroSurface.Swet;

% ------------------------------------
% Trailing Edge sections from Torenbeek
% ------------------------------------
%Df = 0,DF = 45,DF = 105, single/double/triple slotted flaps N/m^2
Df = 0;
Wfte_sfte = 0.65*60*(1+1.6*sqrt(Aircraft_param.Weights.MTOW*9.81*10^(-6))) + Df;
sfte = (1-mean(OPTIM.Wing.a_spar))*2*Aircraft_param.Wing.AeroSurface.Swet;

% ------------------------------------
% Leading Edge high lift devices from Torenbeek
% ------------------------------------
Wslat_sslat = 160*(1+0.7*sqrt(Aircraft_param.Weights.MTOW*9.81*10^(-6)));
sslat = mean(OPTIM.Wing.f_spar)*Aircraft_param.Wing.AeroSurface.Swet;

% ------------------------------------
% Trailing Edge high lift devices from Torenbeek
% ------------------------------------
Ktef = 1.0;
Wtef_stef = 0.65*100*Ktef*(1+sqrt(Aircraft_param.Weights.MTOW*9.81*10^(-6)));

% ------------------------------------
% Aileron Mass from Torenbeek
% ------------------------------------
Wa_sa = 125*(1+0.5*(Aircraft_param.Weights.MTOW*9.81*10^(-6))^(0.25));
sa = 2*0.19*mean(OPTIM.Wing.f_spar)*Aircraft_param.Wing.AeroSurface.Swet;

% ------------------------------------
% Estimate Paint weight
% ------------------------------------
Paint_thickness = 0.15*10^(-3); % (150 microns)
paint_density = 1; % 1.0 Kg/liter
Paint_mass = 2*Paint_thickness*Aircraft_param.Wing.AeroSurface.Swet*1000*paint_density;

% ------------------------------------
% Store Masses in the SecondaryMass structure
% ------------------------------------
SecondaryMass.Slat_mass        = Wslat_sslat*sslat/9.81;
SecondaryMass.FixedTE_mass     = Wfte_sfte*sfte/9.81;
SecondaryMass.FixedLE_mass     = Wfle_sfle*sfle/9.81;
SecondaryMass.TEHighLift_mass  = Wtef_stef*(sfte)/9.81;
SecondaryMass.Aileron_mass     = Wa_sa*(sa)/9.81;
SecondaryMass.Paint_mass       = Paint_mass;
SecondaryMass.Total = SecondaryMass.Slat_mass + SecondaryMass.FixedTE_mass + ...
    SecondaryMass.FixedLE_mass + SecondaryMass.TEHighLift_mass + SecondaryMass.Aileron_mass + SecondaryMass.Paint_mass;

% ------------------------------------
% Linearly Interpolate Volume to the nodes
% ------------------------------------

NodeRefVol = interp1(OPTIM.Wing.y_mbox,OPTIM.Wing.V_int,OPTIM.Wing.nodey,'linear','extrap');

% ------------------------------------
% Store Slat CG locations
% ------------------------------------
SecondaryMass.Slat_cg = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
Slat_refVol = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
NodeIdx = [];
slat_in   = 0.0;
slat_out  = 1.0;
slat_perc = 0.02;
for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    if (slat_in <= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end)) && ...
        (slat_out >= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end))
        
         SecondaryMass.Slat_cg(i) = OPTIM.Wing.nodex(i) - 0.5*(mean(OPTIM.Wing.a_spar) + ...
             mean(OPTIM.Wing.f_spar))*OPTIM.Wing.c_aero(i) + slat_perc*OPTIM.Wing.c_aero(i);
         NodeIdx = [NodeIdx,i];
    end
end
Slat_refVol(NodeIdx(1:end)) = NodeRefVol(NodeIdx(1:end));
SecondaryMass.Slat_node = NodeIdx;

% ------------------------------------
% Store Fixed TE CG locations
% ------------------------------------
SecondaryMass.FixedTE_cg = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
FixedTE_refVol = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
NodeIdx = [];
FixedTE_in   = 0.0;
FixedTE_out  = 0.65;
FixedTE_perc = 0.5*(1-mean(OPTIM.Wing.a_spar)) + mean(OPTIM.Wing.a_spar);
for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    if (FixedTE_in <= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end)) && ...
        (FixedTE_out >= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end))
        
         SecondaryMass.FixedTE_cg(i) = OPTIM.Wing.nodex(i) - 0.5*(mean(OPTIM.Wing.a_spar) + ...
             mean(OPTIM.Wing.f_spar))*OPTIM.Wing.c_aero(i) + FixedTE_perc*OPTIM.Wing.c_aero(i);
         NodeIdx = [NodeIdx,i];
    end
end
FixedTE_refVol(NodeIdx(1:end)) = NodeRefVol(NodeIdx(1:end));
SecondaryMass.FixedTE_node = NodeIdx;

% ------------------------------------
% Store Fixed LE CG locations
% ------------------------------------
SecondaryMass.FixedLE_cg = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
FixedLE_refVol = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
NodeIdx = [];
FixedLE_in   = 0.0;
FixedLE_out  = 1.0;
FixedLE_perc = (2/3)*mean(OPTIM.Wing.f_spar);
for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    if (FixedLE_in <= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end)) && ...
        (FixedLE_out >= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end))
        
         SecondaryMass.FixedLE_cg(i) = OPTIM.Wing.nodex(i) - 0.5*(mean(OPTIM.Wing.a_spar) + ...
             mean(OPTIM.Wing.f_spar))*OPTIM.Wing.c_aero(i) + FixedLE_perc*OPTIM.Wing.c_aero(i);
         NodeIdx = [NodeIdx,i];
    end
end
FixedLE_refVol(NodeIdx(1:end)) = NodeRefVol(NodeIdx(1:end));
SecondaryMass.FixedLE_node = NodeIdx;

% ------------------------------------
% Store TE High Lift CG locations
% ------------------------------------
SecondaryMass.TEHighLift_cg = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
TEHighLift_refVol = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
NodeIdx = [];
TEHighLift_in   = 0.0;
TEHighLift_out  = 0.65;
TEHighLift_perc = (3/4)*(1-mean(OPTIM.Wing.a_spar)) + mean(OPTIM.Wing.a_spar);
for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    if (TEHighLift_in <= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end)) && ...
        (TEHighLift_out >= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end))
        
         SecondaryMass.TEHighLift_cg(i) = OPTIM.Wing.nodex(i) - 0.5*(mean(OPTIM.Wing.a_spar) + ...
             mean(OPTIM.Wing.f_spar))*OPTIM.Wing.c_aero(i) + TEHighLift_perc*OPTIM.Wing.c_aero(i);
         NodeIdx = [NodeIdx,i];
    end
end
TEHighLift_refVol(NodeIdx(1:end)) = NodeRefVol(NodeIdx(1:end));
SecondaryMass.TEHighLift_node = NodeIdx;

% ------------------------------------
% Store Aileron CG locations
% ------------------------------------
SecondaryMass.Aileron_cg = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
Aileron_refVol = zeros(length(beam_model.PartIDs.Wing.BeamNodes),1);
NodeIdx = [];
Aileron_in  = 0.69;
Aileron_out = 0.9;
Aileron_perc = (0.5)*(1-mean(OPTIM.Wing.a_spar)) + mean(OPTIM.Wing.a_spar);

for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    if (Aileron_in <= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end)) && ...
        (Aileron_out >= OPTIM.Wing.nodey(i)/OPTIM.Wing.nodey(end))
        
         SecondaryMass.Aileron_cg(i) = OPTIM.Wing.nodex(i) - 0.5*(mean(OPTIM.Wing.a_spar) + ...
             mean(OPTIM.Wing.f_spar))*OPTIM.Wing.c_aero(i) + Aileron_perc*OPTIM.Wing.c_aero(i);
          NodeIdx = [NodeIdx,i];
    end
end
Aileron_refVol(NodeIdx(1:end)) = NodeRefVol(NodeIdx(1:end));
SecondaryMass.Aileron_node = NodeIdx;

SecondaryMass.Paint_cg = [];
% --------------------------
% Distribute the Wing Secondary Mass
% --------------------------
beam_model.PartIDs.Wing.SecConm2IDs = [];
for i = 1:length(beam_model.PartIDs.Wing.BeamNodes)
    
    lumpedmass_FixedLE    = SecondaryMass.FixedLE_mass*FixedLE_refVol(i)/sum(FixedLE_refVol)/2;
    lumpedmass_FixedTE    = SecondaryMass.FixedTE_mass*FixedTE_refVol(i)/sum(FixedTE_refVol)/2;
    lumpedmass_Slat       = SecondaryMass.Slat_mass*Slat_refVol(i)/sum(Slat_refVol)/2;
    lumpedmass_TEHighLift = SecondaryMass.TEHighLift_mass*TEHighLift_refVol(i)/sum(TEHighLift_refVol)/2;
    lumpedmass_Aileron    = SecondaryMass.Aileron_mass*Aileron_refVol(i)/sum(Aileron_refVol)/2;
    
    FixedLE_offset = SecondaryMass.FixedLE_cg(i)-OPTIM.Wing.nodex(i);
    FixedTE_offset = SecondaryMass.FixedTE_cg(i)-OPTIM.Wing.nodex(i);
    Slat_offset = SecondaryMass.Slat_cg(i)-OPTIM.Wing.nodex(i);
    TEHighLift_offset = SecondaryMass.TEHighLift_cg(i)-OPTIM.Wing.nodex(i);
    Aileron_offset = SecondaryMass.Aileron_cg(i)-OPTIM.Wing.nodex(i);
    
    lumpedmass = (lumpedmass_FixedLE + lumpedmass_FixedTE + lumpedmass_Slat + lumpedmass_TEHighLift + lumpedmass_Aileron);
    
    if lumpedmass > 0
        
        mass_cg = (lumpedmass_FixedLE*FixedLE_offset + lumpedmass_FixedTE*FixedTE_offset + ...
            lumpedmass_Slat*Slat_offset + lumpedmass_TEHighLift*TEHighLift_offset + ...
            lumpedmass_Aileron*Aileron_offset)/lumpedmass;
        
        mass_I = -(lumpedmass_FixedLE*crossm([FixedLE_offset-mass_cg,0,0])*crossm([FixedLE_offset-mass_cg,0,0])+ ...
            lumpedmass_FixedTE*crossm([FixedTE_offset-mass_cg,0,0])*crossm([FixedTE_offset-mass_cg,0,0]) + ...
            lumpedmass_Slat*crossm([Slat_offset-mass_cg,0,0])*crossm([Slat_offset-mass_cg,0,0]) + ...
            lumpedmass_TEHighLift*crossm([TEHighLift_offset-mass_cg,0,0])*crossm([TEHighLift_offset-mass_cg,0,0]) + ...
            lumpedmass_Aileron*crossm([Aileron_offset-mass_cg,0,0])*crossm([Aileron_offset-mass_cg,0,0]));
        
        [beam_model.Conm2,beam_model.Info] = ...
            CONM2Creator(beam_model.Conm2,beam_model.Info,410000 + i,beam_model.PartIDs.Wing.BeamNodes(i),...
            lumpedmass,[mass_cg,0,0],mass_I(1,1),mass_I(2,2),mass_I(3,3));
        
        beam_model.PartIDs.Wing.SecConm2IDs = [beam_model.PartIDs.Wing.SecConm2IDs,410000 + i];
    end
    
end

% --------------------------
% Distribute the HTP Secondary Mass based on 1% of the MTOW
% --------------------------
SecondaryMass.Tail_mass = 0;
if isfield(beam_model.PartIDs,'StbdHTP')
    beam_model.PartIDs.StbdHTP.SecConm2IDs = [];
    for i = 1:length(beam_model.PartIDs.StbdHTP.BeamNodes)
        SecondaryMass.Tail_mass = 0.0025*Aircraft_param.Weights.MTOW; % 1% of the MTOW!!!
        if i == 1
            lumpedmass = SecondaryMass.Tail_mass*OPTIM.StbdHTP.V_int(1)/sum(OPTIM.StbdHTP.V_int)/2;
        elseif i == length(beam_model.PartIDs.StbdHTP.BeamNodes)
            lumpedmass = SecondaryMass.Tail_mass*OPTIM.StbdHTP.V_int(i-1)/sum(OPTIM.StbdHTP.V_int)/2;
        else
            lumpedmass1 = SecondaryMass.Tail_mass*OPTIM.StbdHTP.V_int(i-1)/sum(OPTIM.StbdHTP.V_int)/2;
            lumpedmass2 = SecondaryMass.Tail_mass*OPTIM.StbdHTP.V_int(i)/sum(OPTIM.StbdHTP.V_int)/2;
            lumpedmass = lumpedmass1 + lumpedmass2;
        end
        
        [beam_model.Conm2,beam_model.Info] = ...
            CONM2Creator(beam_model.Conm2,beam_model.Info,420000 + i,beam_model.PartIDs.StbdHTP.BeamNodes(i),...
            lumpedmass,[0,0,0]);
        beam_model.PartIDs.StbdHTP.SecConm2IDs = [beam_model.PartIDs.StbdHTP.SecConm2IDs,420000 + i];
    end
end

SecondaryMass.Tail.Total = 2*SecondaryMass.Tail_mass;

end
