function outline = defoutline(outline,bounds)
% Define the outline that will be used to form the polygons based on
% planform outline. Stretch the outline in all directions to ensure that no
% points are excluded at special regions such as kinks.

halfSpan           = max(outline(:,2));
outline([1,end],2) = outline([1,end],2) + bounds(1) - 0.1;
idx                = find(outline(:,2) >= halfSpan);
ptsTop             = outline(1:idx(1),:);
ptsTop(:,1)        = ptsTop(:,1) - 0.1;
ptsBot             = outline(idx(2):end,:);
ptsBot(:,1)        = ptsBot(:,1) + 0.1;
outline            = [ptsTop;ptsTop(end,:) + [0,0.5,0];ptsBot(1,:) + [0,0.5,0];ptsBot];
outline            = [outline;outline(1,:)];
end
