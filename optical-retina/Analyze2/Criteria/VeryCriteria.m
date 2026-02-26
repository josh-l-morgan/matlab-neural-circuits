

TPN=GetMyDir


%% Read Data
if exist([TPN 'data\BigIT.mat'])
    load([TPN 'data\BigIT.mat'])
else
    'Redo the damn dot find to make BigIT'
end

if exist([TPN 'data\DeltaFill.mat'])
    load([TPN 'data\DeltaFill.mat'])
else
    DFdir=dir([TPN 'pics\DeltaFill']); DFdir=DFdir(3:size(DFdir));
    for i = 1: size(DFdir,1)
        DeltaFill(:,:,i)=imread([TPN 'pics\DeltaFill\' DFdir(i).name]);
        PD=(i/size(DFdir,1))*100
    end
    save([TPN 'data\DeltaFill.mat'],'DeltaFill')
end
load([TPN 'data\Results.mat'])
load([TPN 'data\DotStats.mat'])

%% Look at some stuff
image(sum(BigIT,3)*300), pause
image(sum(DeltaFill,3)*300), pause


%% Do something
Zstart = Results.Arbor(1).Top /.3;
Zstop=Results.Arbor(1).Bottom/.3;

DeltaFills=DeltaFill(:,:,Zstart:Zstop);
%clear DeltaFill
Deltas=DotStats(:,3,1);
Delt=Deltas((DotStats(:,3,2)>Zstart) & (DotStats(:,3,2)<Zstop));
hist(Delt,0:.1:10)
for t = .1:.1:10
   CumulativeThresh(round(t*10))=sum(Delt>t);
end
plot(CumulativeThresh)











