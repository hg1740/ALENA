function [CAERO,INFO,IDs] = SET1Creator(NODE,CAERO,INFO,RBE0,setID)

nset = INFO.nset;
nset = nset + 1;
CAERO.Set.ID(nset) = setID;

CAERO.Set.Node(1,nset).data = NODE(1:length(NODE))';

IDs = setID;
INFO.nset = nset;
end