function [FORCE] = sortFORCEArray(FORCE,NODE,COORD,INFO)

nforce = INFO.nf;
if nforce >0
    
    FORCE = get_load_dir(nforce, FORCE, NODE, COORD);

end

end