function [SPC] = sortSPCArray(PARAM,SPC,NODE,INFO)

nspc = INFO.nspc;

    if nspc > 0
        if PARAM.SPC
            % check if all nodes given are defined
            spcindex = find( SPC.ID == PARAM.SPC);
            
            for j=1:length(spcindex)
                
                [~, i] = intersect(NODE.ID, SPC.Nodes(spcindex(j)).list);
                
                if length(i) ~= length(SPC.Nodes(spcindex(j)).list)
                    error('Unable to determine nodes given in SPC set %d.', SPC.ID(spcindex(j)));
                end
                
                SPC.Nodes(spcindex(j)).list = i;
                
            end
        end
    else
    end

end