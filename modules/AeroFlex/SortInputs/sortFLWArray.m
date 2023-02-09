function [F_FLW] = sortFLWArray(F_FLW,NODE,COORD,INFO)

nflw = INFO.nflw;

if nflw >0  
    F_FLW = get_flw_load_dir(nflw, F_FLW, NODE, COORD);
end

end