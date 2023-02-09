function [RBE2,INFO,IDs] = RBE2Creator(RBE2,INFO,startID,MasterID,SlaveID,DOF,IDs)

nrbe2 = INFO.nrbe2;

nrbe2 = nrbe2 + 1;
RBE2.ID(nrbe2) = startID;
RBE2.IDM(nrbe2) = MasterID;
RBE2.GDL(nrbe2).data = DOF;
RBE2.IDS(nrbe2).data = SlaveID;

IDs = [IDs,RBE2.ID(nrbe2)];
INFO.nrbe2 = nrbe2;
end