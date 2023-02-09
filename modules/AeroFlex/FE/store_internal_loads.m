%STORE_INTERNAL_LOADS: Store the internal loads into a structure 
%
% TODO - See if I can use this everywhere, as I have another similar
% function used in the main solvers
function Optim = store_internal_loads(model,Optim,Parts,trimidx) 

if ~isempty(model.Res.Bar.CForces)
    IntForces = model.Res.Bar.CForces;
    nbar = model.Info.nbar;
else
    IntForces = model.Res.Beam.CForces;
    nbar = model.Info.nbeam;
end

IntBarEnd = zeros(nbar,12);
for i = 1 : nbar
    IntBarEnd(i,1:6)  = IntForces(1,:,i);
    IntBarEnd(i,7:12) = IntForces(2,:,i);
end

IntBarMid = 0.5*(IntBarEnd(:,1:6)+IntBarEnd(:,7:12));

for j = 1:length(Parts)
    
    index = [];
    
    if isfield(model.PartIDs.(Parts{j}),'CBar')
        for i = 1:length(model.PartIDs.(Parts{j}).CBar)
            index = [index,find(model.Bar.ID == model.PartIDs.(Parts{j}).CBar(i))];
        end
    end
    
    if isfield(model.PartIDs.(Parts{j}),'CBeam')
        for i = 1:length(model.PartIDs.(Parts{j}).CBeam)
            index = [index,find(model.Beam.ID == model.PartIDs.(Parts{j}).CBeam(i))];
        end
    end
    if ~isempty(index)
        Optim.(Parts{j}).Fx(:,trimidx) = IntBarMid(index,1); % Axial
        Optim.(Parts{j}).Fy(:,trimidx) = IntBarMid(index,2); % Horizontal
        Optim.(Parts{j}).Fz(:,trimidx) = IntBarMid(index,3); % Vertical
        Optim.(Parts{j}).Mx(:,trimidx) = IntBarMid(index,4); % Torque
        Optim.(Parts{j}).My(:,trimidx) = IntBarMid(index,5); % Out plane Bending Moment
        Optim.(Parts{j}).Mz(:,trimidx) = IntBarMid(index,6); % In plane Bending Moment
    end
    
end

end