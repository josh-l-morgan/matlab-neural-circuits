clear all
TPN = GetMyDir

TPN=[TPN]
Idir=dir(TPN);
In={};
for i = 1: length(Idir)
    name=Idir(i).name;
    LN=length(name);
    if LN>=3
        if length(name)>=4
        if sum(name(LN-3:LN)=='.tif')==4;
            In{length(In)+1}=name;
        end
        end
    end
end

I = imread([TPN In{1}]);
%I = max(I,[],3);
siz = [size(I,1) size(I,2) length(In)];

if siz(3)>1
    I = zeros(siz,'uint8');
for i = 1:siz(3)
    It=imread([TPN In{i}]);
    It=It>0;
    It=sum(It,3);
    I(:,:,i)=It;
end
end


%% Find

for t = 1:3
   Lab=bwlabel(I(:,:,t)>0);
   for i = 1:max(Lab(:))
      [y x] = find(Lab==i);
      Pos(i,:,t) = [mean(y) mean(x)];
   end       
   
   SE = strel('disk',6);
   Er=I(:,:,t)==0; 
   for i = 1:size(I,1)
      ErSave=Er;
      Er=imerode(Er,SE);
      if sum(Er(:))==0,break,end
      %image(Er*1000),pause(.01)
   end
   
   [y x] = find(ErSave>0);
   Cent(t,:)=[mean(y) mean(x)];
   
end

Lab=bwlabel(I(:,:,4)==3);
for i = 1: max(Lab(:))
[y x] = find(Lab==i);
AbPos(i,:) = [mean(y) mean(x)];
end
AbCent=mean(AbPos)


%% Track

%%Get dists
Num=size(Pos,1);
clear Dists
for i = 1:Num
    Dists(i,:)=dist(Pos(:,:,1),Pos(i,:,3));   
end
image(Dists),pause(.01)





%%Find nearest 4
for i = 1:Num
    Dis=Dists(i,:);
    sDis=sort(Dis);
    Pots(i,:)=find(Dis<=sDis(4));
    pDis(i,:)=Dis(Pots(i,:));
end

if length(unique(Pots(:,1)))==size(Pots,1);
    IDs=Pots(:,1);
else
    
    
    
    
    
    
    
end



%% Move
D2Cent1=dist(Pos(:,:,1),AbCent);
D2Cent3=dist(Pos(:,:,3),AbCent);
Mov=D2Cent1-D2Cent3;

hist(Mov)

   Lab1=bwlabel(I(:,:,1)>0);
   Lab3=bwlabel(I(:,:,3)>0);
LabD1=Lab1;
LabD3=Lab3;
for i = 1:max(Lab1(:))
     LabD1(Lab1==i)=D2Cent1(i);
     LabD3(Lab3==i)=D2Cent3(i);
end

image(LabD1+LabD3-200)















