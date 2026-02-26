clear all
TPN = GetMyDir
TPN = [TPN '\'];
% TPN = '\\128.208.64.36\wonglab\Josh\OtherPeoplesData\RachelH\Contract1\'
dTPN=dir(TPN); dTPN=dTPN(3:size(dTPN,1));
dTPN.name

%% Read in time stack

Planes = size(dTPN,1);


for i = 1: Planes
    Is=imread([TPN dTPN(i).name]);
    I(:,:,:,i)=Is;
end

I = I > 100;
I(:,:,2,:)=I(:,:,2,:)|I(:,:,3,:);

%% Find Center
[y x] = find(I(:,:,2,1));
Cent=[mean(y) mean(x)];

%% Draw Well

cmap=jet(100);
cmap(1,:)=0;
cmap(2,:)=.5;
cmap(3,:)=1;
colormap(cmap);
[sy sx sz]=size(Is);
It=zeros(sy,sx);
It(I(:,:,1,1))=2;
for i = 1: Planes  
    It(I(:,:,2,i))=(i-1)*(90/(Planes-1))+10;
    %It(I(:,:,3,i))=(i-1)*(90/(Planes-1))+10;
end

Cent2=round(Cent);
It(Cent2(1)-5:Cent2(1)+5,Cent2(2)-5:Cent2(2)+5)=3;
subplot(2,2,1)
image(It),pause(.1)


%% 
ID = zeros(sy,sx,Planes,'uint16'); %make ID.
for i = 1:Planes
    L1=bwlabel(I(:,:,2,i),8);
    CBnum=max(L1(:));
  
    for l = 1:CBnum
        [y x] = find(L1==l);
        ym=mean(y); xm= mean(x);
        CB(l,1:2,i)=[ym xm];
        Dist=sqrt((ym-Cent(1))^2 + (xm-Cent(2))^2);
        CB(l,3,i)=Dist;
        ID(:,:,i)=ID(:,:,i)+uint16((L1==l)*Dist);
          
    end
   
end


%% ID Cell bodies
CBs=size(CB,1);
Near(1:CBs,1)=1:CBs;
for p = 2: Planes
    for i = 1: CBs  
      if   CB(i,1,p)==0 %if blank spot
          Near(i,p)=0;
      else
          
        Dist=dist(CB(:,1:2,p-1), CB(i,1:2,p));
        Nearest=find(Dist==min(Dist),1);
        Near(i,p)=Near(Nearest,p-1);
      end
   end   
   
end

%% Find dead
%%'Survived' is one for all cells that make it to the final frame
for i = 1: CBs
    Survived(i,1)=~isempty(find(Near(:,Planes)==i,1));    
end

%% Disambiguate
%%Worth writing if nearest neighbor doesnt work out. 

%% Identify controls
Control=logical(zeros(CBs,1));
for i = 1:CBs
   Control(i) = I(round(CB(i,1)),round(CB(i,2)),3,1)==0; 
end
  

%% Track CBs ID
subplot(2,2,2)
CBid=zeros(sy,sx);
CBid2=CBid;
CBid3=CBid2;
SE=strel('disk',20);
for i = 1 : CBs
    CBid3=CBid3*0;
   for p = 1:Planes
        Targ=find(Near(:,p)==i,1);
        if ~isempty(Targ) %& ~(Targ==i)
            CBid2=CBid2*0;
            CBid2(fix(CB(Targ,1,p))+1,fix(CB(Targ,2,p))+1)=1;
            
            CBid2=imdilate(CBid2>0,SE);
            CBid3=CBid3+CBid2;
        end
   end
    CBid3(CBid3>0)=(90/(CBs-1)*(rand*CBs-1))+5;
    CBid=CBid+CBid3;
    image(CBid),pause(.1)
end

image(CBid)


%% Mesure change
for i = 1 : CBs
    for p = 1:Planes-1
        T1=find(Near(:,p)==i);
        T2=find(Near(:,p+1)==i);
    if ~isempty(T1) & ~isempty(T2)
        C=CB(T1,3,p)-CB(T2,3,p+1);
        Change(i,p)=C;
    end % if not empty
    end %Run all Planes p
end % Run all cells = i

ExpChange=Change(~Control); %% Collect change of experimentals
ControlChange=Change(Control); %% Change of Control cells
Changes=sum(Change,2)
subplot(2,2,3)
scatter(Control+1,Changes)
xlim([0 3])


%% Some Stats
dTPN.name
RankSumP=ranksum(Changes(Control & Survived),Changes(~Control & Survived))
%signtest(




