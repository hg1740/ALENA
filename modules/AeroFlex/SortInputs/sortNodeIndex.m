function NODE = sortNodeIndex(NODE,RBE2,BAR,CONM,BEAM,RBE0,CELAS,INFO)

nbar   = INFO.nbar;
nbeam  = INFO.nbeam;
mtot   = INFO.nconm;
ncelas = INFO.ncelas;
nrbe2  = INFO.nrbe2;
ngrid  = INFO.ngrid;
nrbe0  = INFO.nrbe0;

% MODEL dofs. Discard nodes not having mass, stiffness or RBE2
NODE.Index = int32(zeros(1,ngrid));
index = [];
for n = 1:nbar
    index = [index, BAR.Conn(n,:)];
end
NODE.Index(index) = int32(1);

index = [];
for n = 1:nbeam
    index = [index, BEAM.Conn(n,:)];
end
NODE.Index(index) = int32(1);

index = [];
for n = 1:mtot
    index = [index, CONM.Node(n)];
end
NODE.Index(index) = int32(1);

index = [];
for n=1:ncelas
    index = [index, find(NODE.ID == CELAS.Node(n,1))];
    index = [index, find(NODE.ID == CELAS.Node(n,2))];
end
NODE.Index(index) = int32(1);

index = [];
for n=1:nrbe2
    index = [index, find(NODE.ID == RBE2.IDM(n))];
    for k=1:length(RBE2.IDS(n).data)
        index = [index, find(NODE.ID == RBE2.IDS(n).data(k))];
    end
end
NODE.Index(index) = int32(1);

index = sort(find(NODE.Index));
i = [1:length(index)];
NODE.Index(index) = int32(i);

% set reference frame for each node in the basic reference frame
NODE.R = zeros(3, 3, ngrid);
for j=1:ngrid
    NODE.R(:,:,j) = eye(3);
end

% check master aero nodes
for n=1:nrbe0
    
    if NODE.Index(RBE0.Master(n)) == 0
        error('Wrong master node %d for aero set %d: node is not master.', NODE.ID(RBE0.Master(n)), RBE0.ID(n));
    end
    
    index = NODE.Index(RBE0.Node(n).data);
    
    i = find(index);
    
    if ~isempty(i)
        error('Wrong slave node %d for aero set %d.', NODE.ID(RBE0.Node(n).data(i(1))), RBE0.ID(n));
    end
    
end

end