function [PBEAM,INFO] = sortPBeamArray(PBEAM,MAT,PARAM,INFO)

npbeam = INFO.npbeam;

if npbeam >0
      
    [PBEAM.ID, index] = sort(PBEAM.ID);
    
    [labels, i] = unique(PBEAM.ID);
    
    if (length(labels) ~= npbeam) 
        n = [1 : npbeam];
        error('Beam property entries have duplicated labels: %d', PBEAM.ID(setdiff(n, i)));
    end
    
    PBEAM.Mat = PBEAM.Mat(index);
    
    for n=1:npbeam
        
        i = find(MAT.ID == PBEAM.Mat(n));
        if isempty(i)
            error('Unable to find Material %d for Beam property %d.', PBEAM.Mat(n), PBEAM.ID(n));
        else
            PBEAM.Mat(n) = i; % substitute MAT ID with MAT index for rapid accessing
        end
        
    end
    
    PBEAM.A = PBEAM.A(index);
    PBEAM.I = PBEAM.I(index);
    PBEAM.J = PBEAM.J(index);
    PBEAM.RhoNS = PBEAM.RhoNS(index);
    PBEAM.X_L = PBEAM.X_L(index);

    % node propoerties
    PBEAM.Kshear = PBEAM.Kshear(index, :);
    PBEAM.NSI = PBEAM.NSI(index,:);
    PBEAM.NSCG = PBEAM.NSCG(index,:);
    PBEAM.NA = PBEAM.NA(index,:);
    % section propoerties
    PBEAM.DA = PBEAM.DA(:, :, index);
    PBEAM.DI = PBEAM.DI(:, :, index);
    PBEAM.DJ =  PBEAM.DJ(:, :, index);
    PBEAM.DRhoNS =  PBEAM.DRhoNS(:, :, index);
end

end