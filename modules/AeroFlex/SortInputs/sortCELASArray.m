function [CELAS] = sortCELASArray(CELAS,NODE,INFO)

ncelas = INFO.ncelas;

if ncelas>0
    fprintf('\n\tSorting Celas database...');
    CelasID = sort(CELAS.ID);
    if ~isempty(find(diff(CelasID)==0,1))
        error('\n\t: duplicated Celas.');
    end
    for i = 1 : ncelas
        if CELAS.Node(i,1) == CELAS.Node(i,2)
            error('\n\tCelas grids must be different %d. ',CELAS.ID(i));
        end
        CelasDOF1 = zeros(1,6);
        CelasDOF2 = zeros(1,6);
        for j= 1: length(CELAS.DOF(i,1).data)
            CelasDOF1(j) = str2double(CELAS.DOF(i,1).data(j));
        end
        for j= 1: length(CELAS.DOF(i,2).data)
            CelasDOF2(j) = str2double(CELAS.DOF(i,2).data(j));
        end
        CELAS.DOF(i,1).data = NODE.DOF(NODE.ID == CELAS.Node(i,1),CelasDOF1(CelasDOF1~=0));
        CELAS.DOF(i,2).data = NODE.DOF(NODE.ID == CELAS.Node(i,2),CelasDOF2(CelasDOF2~=0));
        
    end
    
    fprintf('done.');
end
end
