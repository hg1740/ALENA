function [BEAM,INFO,IDs] = CBEAMCreator(BEAM,Bar_prop,INFO,startID,NODE)

nbeam = INFO.nbeam;
for i = 1:length(Bar_prop.A)
    nbeam = nbeam + 1;
    BEAM.ID(nbeam) = startID + i;
    BEAM.PID(nbeam) = startID + i;
    BEAM.Conn(nbeam,1) = startID + i;
    BEAM.Conn(nbeam,2) = startID + i + 1;
    BEAM.Conn(nbeam,3) = 0;
    
    % Precalculate the orientation vector based on the direction of beam  
    x1 = NODE.Coord(NODE.ID == startID + i,:) - NODE.Coord(NODE.ID == startID + i + 1,:);
    x2 = [0,0,1];
    x3 = cross(x1,x2);
    if norm(x3) == 0
        error('Wrong Orientation Vector for Beam %i',BEAM.ID(nbeam));
    end
    x2 = cross(x3,x1);
    x2 = x2./norm(x2);
    
    BEAM.Orient(nbeam,1) = x2(1);
    BEAM.Orient(nbeam,2) = x2(2);
    BEAM.Orient(nbeam,3) = x2(3);

    %BEAM.Offset(nbeam,1:9) = zeros(1,9);
    BEAM.Offset(nbeam,1:9) = Bar_prop.ShearOffsets(i,:);
    BEAM.OffsetT(nbeam) = 1;
    BEAM.beamg0(nbeam) = false;
end
IDs = BEAM.ID(INFO.nbeam + 1:INFO.nbeam + length(Bar_prop.A));
INFO.nbeam = nbeam;
end