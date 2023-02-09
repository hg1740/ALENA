function LDW = cross_vseg(VORTEX,mCOLLOC,lemma)
%This function loops throught the corners of each of the vortex filaments
%and constructs the r1 and r2 vectors before passing them to 'fastmega' to
%calculate the downwash terms.

[m1size, m2size, ~] = size(mCOLLOC);
[~, v2size, ~] = size(VORTEX);

LDW = zeros(m1size,m2size,v2size-1,3);  %Edit C.Szczyglowski (9/10/2019) Used to be "LDW = zeros(m1size,m2size,4,3);"
lr1 = zeros(size(mCOLLOC));
lr2 = zeros(size(mCOLLOC));

for j = 1:(v2size-1)

    lr1(:,:,1)=(VORTEX(:,j,1)*lemma)';
    lr1(:,:,2)=(VORTEX(:,j,2)*lemma)';
    lr1(:,:,3)=(VORTEX(:,j,3)*lemma)';

    lr2(:,:,1)=(VORTEX(:,j+1,1)*lemma)';
    lr2(:,:,2)=(VORTEX(:,j+1,2)*lemma)';
    lr2(:,:,3)=(VORTEX(:,j+1,3)*lemma)';
    
    r1=lr1-mCOLLOC;
    r2=lr2-mCOLLOC;
    
    %LDW(:,:,j,:)= mexmega(r1,r2);
    
    LDW(:,:,j,:)= fastmega(r1,r2);
    %     LDW(:,:,j,:)= mega(r1,r2);
    %     F1 = cross_opt(r1,r2,m1size,m2size);
    %     LDW(:,:,j,:) = reshape(F1,m1size,m2size,3);
end

LDW(isnan(LDW(:,:,:,:)))=0;

end