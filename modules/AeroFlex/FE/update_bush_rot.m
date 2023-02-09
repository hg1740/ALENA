

function UPDR = update_bush_rot(nbush, CBUSH, DELTAR, BUSHR)

	UPDR = zeros(3, 3, nbush);
    
    for k=1:nbush
        n1 = CBUSH.Node(k,1);
        n2 = CBUSH.Node(k,2);
        % update nodal sol
        UPDR(:,:,1,k) = DELTAR(:,:,n1) *  BUSHR(:,:,1,k);
        UPDR(:,:,2,k) = DELTAR(:,:,n2) *  BUSHR(:,:,2,k);
    end

end
