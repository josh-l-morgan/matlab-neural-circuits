%function[] = anaGroup(TPN)
colormap colorcube

load([TPN 'Dots.mat'])
load([TPN 'find\SG.mat'])

Pass = SG.passF;

IDs=1:Dots.Num;
IDs=IDs(Pass);
Vox=Dots.Vox(Pass);
ImSize=Dots.ImSize;

%% draw all
field=zeros(ImSize);
for i = 1: length(IDs)
    field(Vox(i).Ind)=IDs(i);
end
[Lab,numLab]=bwlabeln(field,6);
image(max(Lab,[],3))
for i = 1 : numLab
   ids=unique(field(Lab==i));
   if length(ids)>1 % more then one connected
       ItSums=Dots.ItSum(ids);
       ok=ItSums==(max(ItSums)) % Arbitrary !!!!!
       Worthy=ids(find(ok,1)); 
       UnWorthy=ids(~ok);   
       Pass(UnWorthy)=0;
   end
end

SG.PassG=Pass;
save([TPN 'find\SG.mat'],'SG')

maxSum=zeros(Dots.ImSize(1),Dots.ImSize(2));
maxID1=maxSum; maxID2 = maxSum; maxID3 = maxSum;

P=find(Pass);
for i = 1:length(P)

   maxID1(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;  
    maxID2(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;  
    maxID3(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=rand*100+50;   
    maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)=maxSum(mod(Dots.Vox(P(i)).Ind-1,YXsize)+1)+1;
end

maxPassed2=uint8((maxSum>0)*200);

MaxC=maxID1+(maxSum>1)*1000;
MaxC(:,:,2)=maxID2+(maxSum>1)*1000;
MaxC(:,:,3)=maxID3+(maxSum>1)*1000;
MaxC=uint8(MaxC);
image(MaxC),pause(.01)


%% load Raw
load([TPN 'images\maxRaw.mat']) 
Grouped=maxRaw;
Grouped(:,:,3)=maxSum*200;
imwrite(Grouped,[TPN 'find\Grouped.tif'],'Compression','none')
image(Grouped),pause(.01)





