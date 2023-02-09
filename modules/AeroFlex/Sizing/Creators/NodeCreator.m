function [NODE,INFO,IDs]  = NodeCreator(NODE,X,Y,Z,INFO,startID)

ngrid = INFO.ngrid;

for i = 1:length(Y)
    ngrid = ngrid +1;
    NODE.ID(i + INFO.ngrid) = startID + i;
    NODE.Coord(i + INFO.ngrid,1:3) = [X(i),Y(i),Z(i)];
    NODE.CS(i + INFO.ngrid) = 0;
    NODE.CD(i + INFO.ngrid) = 0;
end

IDs = NODE.ID(INFO.ngrid + 1:INFO.ngrid + length(Y));
INFO.ngrid = ngrid;
end

