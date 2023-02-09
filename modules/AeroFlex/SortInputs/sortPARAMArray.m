function [PARAM] = sortPARAMArray(PARAM,NODE)

if PARAM.GRDPNT ~= 0
    idx = find(NODE.ID == PARAM.GRDPNT);
    PARAM.GRDPNT = idx;
end

end