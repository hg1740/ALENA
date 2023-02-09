function [SPC,INFO] = SPCCreator(SPC,INFO,SPCID,SPCDOF,SPCNode)

nspc = INFO.nspc;

nspc = nspc + 1;

SPC.ID(nspc) = SPCID;
constr = num2str(SPCDOF);
SPC.DOF(nspc).list = zeros(1, length(constr));

if length(constr) > 6
    fclose(fp);error('Constrained too many dofs in SPC set %d.', SPC.ID(nspc));
end

for i = 1:length(constr)
    SPC.DOF(nspc).list(i) = int32(str2double(constr(i)));
end

SPC.Nodes(nspc).list(nspc) = SPCNode;
SPC.Nodes(nspc).list = sort(SPC.Nodes(nspc).list);
j = find(SPC.Nodes(nspc).list);
SPC.Nodes(nspc).list = SPC.Nodes(nspc).list(j);

INFO.nspc = nspc;

end