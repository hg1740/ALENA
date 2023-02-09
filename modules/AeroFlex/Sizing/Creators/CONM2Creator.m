function [CONM2,INFO] = CONM2Creator(CONM2,INFO,ID,NODE,PAYLOAD,OFFSET,Ixx,Iyy,Izz)

ncom2 = INFO.ncom2;
ncom2 = ncom2 +1;

CONM2.ID(ncom2) = ID;
CONM2.Node(ncom2) = NODE;
CONM2.CID(ncom2) = 0;
CONM2.M(:,:,ncom2) = zeros(6,6);
CONM2.M(1,1,ncom2) = PAYLOAD;
CONM2.M(2,2,ncom2) = CONM2.M(1,1,ncom2);
CONM2.M(3,3,ncom2) = CONM2.M(1,1,ncom2);
switch nargin
    case 9
        CONM2.M(4,4,ncom2) = Ixx;
        CONM2.M(5,5,ncom2) = Iyy;
        CONM2.M(6,6,ncom2) = Izz;
end

CONM2.Offset(ncom2,1:3) = OFFSET;

INFO.ncom2 = ncom2;
end

