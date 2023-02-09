classdef FEModelAircraft < awi.fe.FEModel
    
    methods
        function updateDmiEntry(obj)
            ch = flatlist(obj);
            res = {};
            for i = 1:length(ch)
               if ~isempty(ch(i).AeroPanels) 
                   res{end+1} = ch(i).AeroPanels.definePanels;
               end
            end
            AoAs = [];
            IDs = [];
            for i = 1:length(res)
                % deal with nan AoA
                aoa = res{i}.AoA(:);
                aoa(isnan(aoa)) = 0;
                AoAs = [AoAs;aoa];
                IDs = [IDs;res{i}.panelID(:)];
            end
            [~,I] = sort(IDs);
            if isempty(obj.DMI)
                DMI = awi.fe.DMI;
                DMI.NAME = 'W2GJ';
                DMI.DATA = [];

                %Set up parts
                addPart(obj, 'DMI', DMI);

                %Add the objects to the FE model
                addFEData(obj, DMI);
            end
            obj.DMI.DATA = AoAs(I);
        end
    end
end

