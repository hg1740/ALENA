function [coords, CaB] = strains2coords_all(x, Matrices)
%strains2coords_all Returns the coordinates and orientations at each
%structural node based on the current value of the local beam loads (forces
%and moments) and the initial position/orientation of the undeformed beam.

nPart = numel(Matrices.n_elem);

%Preallocate
nS     = cellfun(@numel, Matrices.s_out);
coords = arrayfun(@(nS) zeros(3, 1, nS), nS, 'Unif', false);
CaB    = arrayfun(@(nS) zeros(3, 3, nS), nS, 'Unif', false);

%Direction of axial extension
e1 = [1 ; 0 ; 0];
%e1 = zeros(3,1);e1(1) = 1;

%For each part...
for jj = 1:nPart
        
    %Part data
    CaB{jj}(:,:,1) = Matrices.CaB0{jj};     %Initial orientation of part w.r.t aircraft
    nElem          = Matrices.n_elem(jj);    
            
    %Indexing    
    counter  = 1;
    counter2 = 1;
    ub = cumsum(repmat(6, [1, nElem]));
    lb = [1, ub(1 : end - 1) + 1];
    
    %Multiply forces & moments by compliance matrix to obtain strains and
    %curvatures 
    % -> straincurv = [eps_x ; eps_y ; eps_z ; kappa_x ; kappa_y ; kappa_z]
    straincurv = Matrices.Cfull{jj} * x{jj};
    
    for ii = 1 : nElem
        
        %Grab strains and curvatures for this element with pre-curve
        eps_i   = straincurv(lb(ii) : ub(ii) - 3)      + Matrices.eps0{jj}(:,:,ii)   + e1;
        kap_i   = skew(straincurv(ub(ii) - 2 : ub(ii)) + Matrices.kappa0{jj}(:,:,ii))    ;
        
        %Which nodes lie in the element
        idx = and( ...
            Matrices.s_out{jj} >  Matrices.s{jj}(ii), ...
            Matrices.s_out{jj} <= Matrices.s{jj}(ii + 1));
        
        %Calculate distance along each element for points along element
        ds_node = (Matrices.s_out{jj}(idx) - Matrices.s{jj}(ii));
        
        %Calculate new coordinates and orientations for each node
        %   - N.B. Strains and curvatures are constant across the element 
        for kk = 1 : nnz(idx)
            
            %Exponential term
            expon   = [zeros(1, 4); eps_i, kap_i] * ds_node(kk);
            
            %         expon_tot = expon_tot + expon;
            coord_CaB = expm(expon')*[coords{jj}(:,:,counter2)';Matrices.CBB{jj}(:,:,ii)*CaB{jj}(:,:,counter2)'];
            CaB{jj}(:,:,counter+1)    = coord_CaB(2:end,:)';
            coords{jj}(:,:,counter+1) = coord_CaB(1,:)';
            
            counter = counter + 1;
        end
        counter2 = counter;
    end
end

end
% coord_CaB_final = [coords(:,:,1),CaB(:,:,1)]*expm(expon_tot);