%% Draw BI

TPN = GetMyDir

load([TPN 'Use.mat'])
load([TPN 'ONOFF.mat'])

Mids=Use.Mids;
DPos=Use.DPos;
MidLayer=ONOFF.MidLayer;
DotLayer=ONOFF.DotLayer;

Mids=fix(Mids*2)+1;
DPos=fix(DPos*2)+1;

maxMids=max(max(Mids),max(DPos));


ShowBI=zeros(maxMids,'uint8');
for i = 1:size(Mids,1)
    ShowBI(Mids(i,1),Mids(i,2),Mids(i,3))=1+MidLayer(i);
end
for i = 1:size(DPos,1)
    ShowBI(DPos(i,1),DPos(i,2),DPos(i,3))=4+DotLayer(i);
end
imwriteNp(TPN,ShowBI,'ShowBI')