function [PBEAM,INFO,IDs] = PBEAMCreator(PBEAM,Bar_prop,INFO,startID)

npbeam = INFO.npbeam;
for i = 1:length(Bar_prop.A)
    npbeam = npbeam + 1;
    PBEAM.ID(npbeam) = startID + i;
    % Assumes material ID = 1;
    PBEAM.PID(npbeam) = 1;
    
    %ENDA
    PBEAM.A(npbeam).data(1)   = Bar_prop.A(i);
    PBEAM.I(npbeam).data(1,1) = Bar_prop.I1(i);
    PBEAM.I(npbeam).data(2,1) = Bar_prop.I2(i);
    PBEAM.I(npbeam).data(3,1) = 0*Bar_prop.I1(i);
    PBEAM.J(npbeam).data(1)   = Bar_prop.J(i);
    PBEAM.RhoNS(npbeam).data(1) = 0.0;
    
    %ENDB
    PBEAM.A(npbeam).data(2)   = Bar_prop.A(i);
    PBEAM.I(npbeam).data(1,2) = Bar_prop.I1(i);
    PBEAM.I(npbeam).data(2,2) = Bar_prop.I2(i);
    PBEAM.I(npbeam).data(3,2) = 0*Bar_prop.I1(i);
    PBEAM.J(npbeam).data(2)   = Bar_prop.J(i);
    PBEAM.RhoNS(npbeam).data(2) = 0.0;
    
    PBEAM.Mat(npbeam) = Bar_prop.Mat;
    PBEAM.X_L(npbeam).data(1) = 0;
    PBEAM.X_L(npbeam).data(2) = 1;
    PBEAM.Kshear(npbeam,1:2) = zeros(1,2);
    PBEAM.NSI(npbeam,1:2)  = zeros(1,2);
    PBEAM.NSCG(npbeam,1:4) = zeros(1,4);
    PBEAM.NA(npbeam,1:4)   = zeros(1,4);
    PBEAM.Str_point(:,:,npbeam) = zeros(4,2);
    PBEAM.SI(npbeam) = 0;
    PBEAM.Type(npbeam) = 0;

    PBEAM.DA(:, :, npbeam) = zeros(1,2);
    PBEAM.DI(:, :, npbeam) = zeros(3,2);
    PBEAM.DJ(:, :, npbeam) =  zeros(1,2);
    PBEAM.DRhoNS(:, :, npbeam) =  zeros(1,2);
    
    PBEAM.DA(1, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.A(npbeam).data,	  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    PBEAM.DI(1, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.I(npbeam).data(1,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    PBEAM.DI(2, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.I(npbeam).data(2,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    PBEAM.DI(3, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.I(npbeam).data(3,:), [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    PBEAM.DJ(1, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.J(npbeam).data,	  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
    PBEAM.DRhoNS(1, 1:2, npbeam) = interp1(PBEAM.X_L(npbeam).data,...
        PBEAM.RhoNS(npbeam).data,  [0.5 - 1/(2*sqrt(3)), 0.5 + 1/(2*sqrt(3))]);
        
end

IDs = PBEAM.ID(INFO.npbeam + 1:INFO.npbeam + length(Bar_prop.A));
INFO.npbeam = npbeam;
end