function[] = anaGroup(TPN)
colormap colorcube

%%load data
load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])

%%extract info
Pass = SG.passF;
IDs=1:Dots.Num;
IDs=IDs(Pass);
Vox=Dots.Vox(Pass);
ImSize=Dots.ImSize;

%% Group
field=zeros(ImSize);
%%Draw data
for i = 1: length(IDs)
    field(Vox(i).Ind)=IDs(i);
end
[Lab,numLab]=bwlabeln(field,6);
image(max(Lab,[],3)*20)
%%Get IDs

memIds=[];
for i = 1 : numLab
    ids{i,1}=unique(field(Lab==i)); %record source ID from dots
    memIds=[memIds ; ids{i,1}];
end

NumberObjects=numLab
NumberPeaks=length(memIds)
NumberUniquePeaks=length(unique(memIds))

if NumberPeaks~=NumberUniquePeaks
    'Danger !!! peaks are being missed or counted multiple times.',...
        ' Please rerun dot finder (anaDFc)'
end

Grouped.ids=ids;

Voxs=cell(length(ids),1);
DPos=zeros(length(ids),3);
for i = 1:length(ids)
    gid =ids{i};
    for g= 1:length(gid)
        Voxs{i}=[Voxs{i} ; Dots.Vox(gid(g)).Pos];
    end
    DPos(i,:)=mean(Voxs{i},1);
    
end

Grouped.Vox=Voxs;
Grouped.DPos=DPos;

save([TPN 'Grouped.mat'],'Grouped')



%% Draw 
for f = 1:10
maxSum=zeros(Dots.ImSize(1),Dots.ImSize(2));
maxID1=maxSum; maxID2 = maxSum; maxID3 = maxSum;

%P=361
YXsize=Dots.ImSize(1)*Dots.ImSize(2);

for i = 1:length(Grouped.DPos)
    Pos=Grouped.Vox{i,:};
    PosI=sub2ind(ImSize(1:2),Pos(:,1),Pos(:,2));
    maxID1(PosI)=rand*200+50;  
    maxID2(PosI)=rand*200+50;  
    maxID3(PosI)=rand*200+50;   
    maxSum(PosI)=maxSum(PosI)+1;
end

maxPassed2=uint8((maxSum>0)*200);

MaxC=maxID1+(maxSum>1)*1000;
MaxC(:,:,2)=maxID2+(maxSum>1)*1000;
MaxC(:,:,3)=maxID3+(maxSum>1)*1000;
MaxC=uint8(MaxC);
image(MaxC),pause(.01)
GroupI=MaxC;
end

%% load Raw
% load([TPN 'images\maxRaw.mat']) 
% GroupI=maxRaw;
% GroupI(:,:,3)=maxSum*200;
imwrite(GroupI,[TPN 'find\GroupI.tif'],'Compression','none')
image(GroupI),pause(.01)





