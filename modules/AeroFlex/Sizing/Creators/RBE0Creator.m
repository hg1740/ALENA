function [RBE0,INFO,IDs] = RBE0Creator(RBE0,INFO,startID,MasterID,SlaveLE,SlaveTE,NumGrids)

nrbe0 = INFO.nrbe0;
for i = 1:NumGrids
    nrbe0 = nrbe0 + 1;
    RBE0.ID(nrbe0) = startID + i;
    RBE0.Master(nrbe0) = MasterID + i;
    RBE0.Node(nrbe0).data = [SlaveLE + i,SlaveTE + i];
end
IDs = RBE0.Master(INFO.nrbe0 + 1:INFO.nrbe0 + NumGrids);
INFO.nrbe0 = nrbe0;
end