function [PBAR,INFO] = sortPBarArray(PBAR,MAT,PARAM,INFO)

npbar = INFO.npbar;
    if npbar >0
              
        [PBAR.ID, index] = sort(PBAR.ID);
        [labels, i] = unique(PBAR.ID);
        
        if (length(labels) ~= npbar)
            
            n = [1 : npbar];
            error('Bar property entries have duplicated labels: %d', PBAR.ID(setdiff(n, i)));
            
        end
        
        %
        PBAR.Mat = PBAR.Mat(index);
        PBAR.Type = PBAR.Type(index);
        PBAR.SI = PBAR.SI(index);
        
        for n=1:npbar
            
            i = find(MAT.ID == PBAR.Mat(n));
            if isempty(i)
                error('Unable to find Material %d for Bar property %d.', PBAR.Mat(n), PBAR.ID(n));
            else
                PBAR.Mat(n) = i; % substitute MAT ID with MAT index for rapid accessing
            end
            
        end
        
        PBAR.A = PBAR.A(index);
        PBAR.I = PBAR.I(index,:);
        PBAR.J = PBAR.J(index);
        PBAR.RhoNS = PARAM.WTMASS .* PBAR.RhoNS(index);
        PBAR.Kshear = PBAR.Kshear(index,:);
        PBAR.Str_point = PBAR.Str_point(:,:,index);
        PBAR.RhoNSV = zeros(size(PBAR.RhoNS));
         
    end
end