function [dw,DW] = fastdw(lattice,DW)

one_by_four_pi=1/(4*pi);

[psize, ~, ~] = size(lattice.VORTEX);
lemma         = ones(1,psize);

if isempty(DW)
    
    mCOLLOC        = zeros(psize, psize, 3);
    mCOLLOC(:,:,1) = lattice.COLLOC(:,1)*lemma;
    mCOLLOC(:,:,2) = lattice.COLLOC(:,2)*lemma;
    mCOLLOC(:,:,3) = lattice.COLLOC(:,3)*lemma;
    
    [csize, ~, ~]  = size(lattice.COLLOC);
    lemmac         = ones(1,csize);
    %LDW = cross_vseg_mex(lattice.VORTEX,mCOLLOC,lemmac);
    
    LDW = cross_vseg(lattice.VORTEX,mCOLLOC,lemmac);
    
    %Move into parent function
    LDW((isnan(LDW(:,:,:,:))))=0;
    %DDW(find((isnan(DDW(:,:,:,:)))))=0;
    
    DW = -squeeze(sum(LDW,3))*one_by_four_pi;
    %DW_TV=-squeeze(sum(DDW,3))*one_by_four_pi;
end

mN(:,:,1) = lattice.N(:,1)*lemma;
mN(:,:,2) = lattice.N(:,2)*lemma;
mN(:,:,3) = lattice.N(:,3)*lemma;

dw = sum(DW.*mN,3);
%dw_TV=sum(DW_TV.*mN,3);
end
