function plot_beam(beam_model,NODE)

figure(444);

% ---------------------- %
% Plot Bar nodes
% ---------------------- %
BarNodes = [];

for i = 1: beam_model.Info.nbar
    BarNodes = [BarNodes; beam_model.Bar.Conn(i,:)'];
    hold on;
    plot3(NODE(beam_model.Bar.Conn(i,:),1),NODE(beam_model.Bar.Conn(i,:),2),...
        NODE(beam_model.Bar.Conn(i,:),3),'k-','LineWidth',1.5);
    hold on;
    plot3(NODE(beam_model.Bar.Conn(i,:),1) + beam_model.Bar.Offset(i,[1,4,7])',NODE(beam_model.Bar.Conn(i,:),2) + beam_model.Bar.Offset(i,[2,5,8])',...
        NODE(beam_model.Bar.Conn(i,:),3) + beam_model.Bar.Offset(i,[3,6,9])','r-','LineWidth',1.5);
end

BarNodes = sort(unique(BarNodes));
hold on;
plot3(NODE(BarNodes,1),NODE(BarNodes,2),NODE(BarNodes,3),'ko');

% ---------------------- %
% Plot Beam nodes
% ---------------------- %
BeamNodes = [];
for i = 1: beam_model.Info.nbeam
    BeamNodes = [BeamNodes; beam_model.Beam.Conn(i,:)'];
    hold on;
    plot3(NODE(beam_model.Beam.Conn(i,:),1),NODE(beam_model.Beam.Conn(i,:),2),...
        NODE(beam_model.Beam.Conn(i,:),3),'k-','LineWidth',1.5);
     hold on;
    plot3(NODE(beam_model.Beam.Conn(i,:),1) + beam_model.Beam.Offset(i,[1,4,7])',NODE(beam_model.Beam.Conn(i,:),2) + beam_model.Beam.Offset(i,[2,5,8])',...
        NODE(beam_model.Beam.Conn(i,:),3) + beam_model.Beam.Offset(i,[3,6,9])','r-','LineWidth',1.5);
end

BeamNodes = sort(unique(BeamNodes));
hold on;
plot3(NODE(BeamNodes,1),NODE(BeamNodes,2),NODE(BeamNodes,3),'ko');

% RBE0s
if beam_model.Info.nrbe0>0
    for i = 1:beam_model.Info.nrbe0
        MNodeid = beam_model.RBE0.Master(i);
        for j = 1:length(beam_model.RBE0.Node(i).data)
            SNodeid = beam_model.RBE0.Node(i).data(j);
            hold on;
            plot3(NODE([MNodeid; SNodeid],1),NODE([MNodeid; SNodeid],2),...
               NODE([MNodeid; SNodeid],3),'r--','LineWidth',1.5);
            plot3(NODE([SNodeid],1),NODE([SNodeid],2),...
                NODE([SNodeid],3),'r*');
        end
        
    end
end
% RBE2s
if beam_model.Info.nrbe2>0
    for i = 1:beam_model.Info.nrbe2
        MNodeid = find(beam_model.Node.ID == beam_model.RBE2.IDM(i));
        for j = 1:length(beam_model.RBE2.IDS(i).data)
            SNodeid = find(beam_model.Node.ID == beam_model.RBE2.IDS(i).data(j));
            hold on;
            plot3(NODE([MNodeid; SNodeid],1),NODE([MNodeid; SNodeid],2),...
                NODE([MNodeid; SNodeid],3),'b--','LineWidth',2);
        end
        
    end
end
end