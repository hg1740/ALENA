function PART = InitialiseOptimVariables(PART)

if PART.Type == 2
    PART.SOL = zeros(2*PART.Nsec,1);
    PART.SOL(1:PART.Nsec) = PART.t_skin;
    PART.SOL(PART.Nsec+1:end) = PART.t_spar;
elseif PART.Type == 3
    PART.SOL = zeros(3*PART.Nsec,1);
    PART.SOL(1:PART.Nsec)  = PART.t_skin;
    PART.SOL(PART.Nsec+1:2*PART.Nsec)  = PART.t_spar;
    PART.SOL(2*PART.Nsec+1:end)  = PART.A_str;
elseif PART.Type == 4
    PART.SOL = zeros(3*PART.Nsec,1);
    PART.SOL(1:PART.Nsec) = PART.t_skin;
    PART.SOL(2*PART.Nsec+1:end) = PART.t_str;
    PART.SOL(PART.Nsec+1:2*PART.Nsec) = PART.t_spar;
end

end