function [coords,CaB] = strains2coords(x,n_elem,s)

% error('Not the latest version!')

CaB    = zeros(3,3,n_elem+1);
coords = zeros(3,1,n_elem+1);

CaB(:,:,1) = eye(3);

e1 = zeros(3,1);e1(1) = 1;

expon_tot = zeros(4,4);

for ii = 1:n_elem
    ind  = [1:6] + (ii-1)*6;
    eps_i   = x(ind(1:3))  + e1;
    kap_i   = skew(x(ind(4:6)));
    
    expon   = [zeros(1,4);eps_i,kap_i]*(s(ii+1)-s(ii));
    
    expon_tot = expon_tot + expon;    
    coord_CaB = expm(expon')*[coords(:,:,ii)';CaB(:,:,ii)'];
    CaB(:,:,ii+1)    = coord_CaB(2:end,:)';
    coords(:,:,ii+1) = coord_CaB(1,:)';
end

coord_CaB_final = [coords(:,:,1),CaB(:,:,1)]*expm(expon_tot);