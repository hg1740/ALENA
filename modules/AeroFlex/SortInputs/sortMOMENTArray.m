function [MOMENT] = sortMOMENTArray(MOMENT,NODE,COORD,INFO)

nmoment = INFO.nm;
if nmoment >0
    
    MOMENT = get_load_dir(nmoment, MOMENT, NODE, COORD);

end

end