%%-Draw Autocorrelegram for DRP

TPN = GetMyDir

load([TPN 'Dat.mat'])




Pos=Dat.Raw.Pos;
scatter(Pos(:,1),Pos(:,2))

maxpos=max(Pos);
R=50;
hold off
for i = 1:size(Pos,1)
   pos=Pos(i,:); 
   if (pos(1)>R) && (pos(2)>R) && (pos(1)< maxpos(1)-R) && (pos(2)<maxpos(2)-R); 
       Dist=dist(Pos,pos);
       Show=Dist<=R;
       scatter(Pos(Show,1)-pos(1),Pos(Show,2)-pos(2),'o','k')
       hold on
    
   end
    
end

xlim([-R R])
ylim([-R R])