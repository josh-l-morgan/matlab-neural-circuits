            

[TFN TPN] = GetMyFile;
TX = load([TPN TFN]);
Point=[mean([TX(:,17),TX(:,23)],2) mean([TX(:,18),TX(:,24)],2)];
%Point = Point * yxum; %scale for resolution
gains= find(TX(:,4)==65280);  %find all gains (green arrows)
losses = find(TX(:,4) == 255); %find all losses (red arrows)
stable = find(TX(:,4) ==16776960) ; %find all non gain loss puncta

Pos = [(1:size(Point,1))' Point]
