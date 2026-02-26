function[] = anaGroup(TPN)
colormap colorcube

%%load data
load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])

%%collect information from dots that passed SG thresholding
Pass = SG.passF;
IDs=find(Pass);
ImSize=Dots.ImSize;
Vox = Dots.Vox;

nearby = 10;  %Distance at which faces are looked at
maxContact = 0.3; % maximum faces/volume to be independent dot

%% Determine distance from dot to dot
Dists = zeros(length(IDs),length(IDs));
for i = 1: length(IDs)
   Dists(i,:) = dist(Dots.Pos(IDs,:),Dots.Pos(IDs(i),:));
  
end
[x y] = find((Dists > 0) & (Dists < nearby));

idX=IDs(x);
idY=IDs(y);

%% Count shared faces of nearby dots
faces = [];
for i = 1: length(idX)
   cDotX = Dots.Vox(idX(i)).Pos;
   cDotY = Dots.Vox(idY(i)).Pos;
   face = zeros(size(cDotX,1),1);
   for v = 1: size(cDotX,1)
    diff = [cDotY(:,1)-cDotX(v,1)  cDotY(:,2)-cDotX(v,2) cDotY(:,3)-cDotX(v,3)];
    face(v) = sum(sum(abs(diff),2)<2);
   end  
    faces(i) = sum(face);
end

%% form new groups is too much contact
newGroup = [];
for i = 1:length(faces);
   Vol = Dots.Vol(idX(i));
   contactRat(i) = faces(i)/Vol;
   if contactRat(i) > maxContact 
    newGroup(size(newGroup,1)+1,:) = [idX(i) idY(i)];
   end   
end

%% Organize groups
depleteIDs= newGroup;
groups = {};
gNum = 0;
while sum(depleteIDs(:))
   gNum = gNum+1;
   newest = depleteIDs(find(depleteIDs>0,1));
   foundNew = 1;
   groups{gNum} = newest;
   while foundNew
       foundNew = 0;
       grby = [];
       for n = 1: length(newest)
            [grby grbx] = find(depleteIDs == newest(n));
            graby = [grby;grby];
       end
       if ~isempty(graby)
           grabbed = newGroup(graby,:);
           grabbed = unique(grabbed(:));
           grabbed = setdiff(grabbed, newest);
           groups{gNum} = [groups{gNum} grabbed'];   
           depleteIDs(graby,:)=0;
           newest = grabbed;
           foundNew = 1;
       end
   end
end

nonGroupedIDs = setdiff(IDs,newGroup(:));
for i = 1:length(nonGroupedIDs)
   groups{length(groups)+1} = nonGroupedIDs(i); 
end



%% Create merged dots 
Grouped.ids=groups;
Voxs=cell(length(groups),1);
DPos=zeros(length(groups),3);
for i = 1:length(groups)
    gid =groups{i};
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

for i = 1:size(Grouped.DPos,1)
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





