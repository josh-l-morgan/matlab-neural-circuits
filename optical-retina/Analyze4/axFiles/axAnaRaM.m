function[]=axAnaRaM(TPN,DPN)

%% find ratioed puncta brightness

%Puncta variable should contain yxz position, puncta number, brightness of
%dend channel and brightness of green channel.
% Dendrite matrix should include  position and both channels

%record the raw and median filtered

%matched variables
%vDot in the form [Ir Irm Ig Igm]
%DotID in the form [ID]
%DotPos in the form [y x z]
%Puncta in the form [y x z ID]

%Output Variables
%Center = [y x z] in pixels
%Stats = [ID volume(in pixels) Delta-F-over-f]
%DotStats = [ Stats ; Center ; Center in um]

colormap gray(255)
tic

%DPN = GetMyDir
%Get directory name
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)


%% Vectorize BigFilled
load([TPN 'data/Dots.mat']);


%%READ IMAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Get number of planes and spacer zeros
d=dir(DPN); %get number of files in directory
d=d(3:size(d,1)); %find number of planes
maskdir=dir([TPN 'mask']);maskdir=maskdir(3:size(maskdir,1));


%%Figure out channels

checkPlane = round(size(d,1)/2);
Ic(:,:,:)=imread([DPN d(checkPlane).name]); %read
if size(Ic,3)==1, 'crash= only one channel read' %if only one channel
elseif size(Ic,3)==2, rchannel=2; gchannel=1; %if only two channels
elseif sum(sum(Ic(:,:,3)))==0, rchannel=2; gchannel=1; %if third channel is blank
else rchannel=1; gchannel=2; %if third channel is not blank
end


clear I Ic pred
vDend=zeros(1,2); %start counter for vectorized dendrites
vDot=zeros(1,2); %start counter for vectorized dendrites
DotID=0;
DotPos=zeros(1,3); %start matrix for dot positions
for i=1:size(d,1)
    Ic(:,:,:)=imread([DPN d(i).name]); %read
    Ir=Ic(:,:,rchannel);
    Ig=Ic(:,:,gchannel);
    Imask=imread([TPN 'mask\' maskdir(i).name]);
    Irm=medfilt2(Ir,[3,3]); %median filter
    Igm=medfilt2(Ig,[3,3]); %median filter

    %Get info for each dot voxel and draw Dot map (IDs) for plane
    IDs=zeros(Dots.ImSize(1:2));
    for ds=1:Dots.Num
        for v=1:size(Dots.Vox(ds).Pos,1)
            if Dots.Vox(ds).Pos(v,3)==i
                IDs(Dots.Vox(ds).Pos(v,1),Dots.Vox(ds).Pos(v,2))=ds;
                Dots.Vox(ds).Irm(v)=Irm(Dots.Vox(ds).Pos(v,1),Dots.Vox(ds).Pos(v,2));
                Dots.Vox(ds).Igm(v)=Igm(Dots.Vox(ds).Pos(v,1),Dots.Vox(ds).Pos(v,2));
            end
        end
    end


    %Find Dendrite Values
    %     Mask=Irm>=0; % Threshold plane to make Mask
    Mask=Irm>0;
    Mask(IDs>0)=0;
    if any(Mask(:))
        if exist('pred','var')
            next=size(pred,1)+1;
            p = polyfit(double(Irm(Mask)), double(Igm(Mask)),1);
            pred(next,:) = polyval(p,0:255);
        else
            p = polyfit(double(Irm(Mask)), double(Igm(Mask)),1);
            pred = polyval(p,0:255);
        end
    end
    PercentRead=i/size(d,1)*100
end
mPred = mean(pred);
% %% make look up table for predicted Green value for ever red value
% %%find R-G intensity map
% allPreds=sum(NumPreds,2);
% RtoGall=meanPredm.* NumPreds; %weight ratios
% RtoG=sum(RtoGall,2);
% RtoG(RtoG>0)=RtoG(RtoG>0)./allPreds(RtoG>0);
% Rvals=(0:255)';
%
% for i = 1 : 256
%     if allPreds(i)<100
%         for s= 10:200
%             OKpreds=allPreds>100;
%             OKpreds=OKpreds & (Rvals>(i-20) & Rvals< (i+ 20));
%             Amt=allPreds(Rvals>(i-20) & Rvals< (i+ 20));
%             if sum(Amt)>500,break,end
%         end
%         [p, S]=polyfit(Rvals(OKpreds),RtoG(OKpreds),1);
%         RtoG(i)=i*p(1)+p(2);
%     end
% end
%
% %Gfit=Rvals*p(1)+p(2);
% %plot(RtoG)
%
%
% plot(meanPredm,'r'),pause(1)
% Dots.Im.GreenFromRed=meanPredm;
%
% ImMode=find(allPreds==max(allPreds));
% RBG=ImMode-1;
% GBG=RtoG(ImMode);
% Dots.Im.GreenBackGround=GBG;
%
%
% Dots.Im.GreenFromRed=meanPredm;

%% Find values of each puncta

% for i=1:Dots.Num
%     for v=1:size(Dots.Vox(i).Pos,1)
%         Pred=double(meanPredm(Dots.Vox(i).Irm(v)+1));
%         Dots.Vox(i).DF(v)=(double(Dots.Vox(i).Igm(v))-Pred);
%         Dots.Vox(i).DFOf(v)=(double(Dots.Vox(i).Igm(v))-Pred)/(max(1,Pred-GBG));
%     end
%     Dots.DF(i)=mean(Dots.Vox(i).DF);
%     Dots.DFOf(i)=mean(Dots.Vox(i).DFOf);
%     Dots.DFOfTopHalf(i)=mean(Dots.Vox(i).DFOf(Dots.Vox(i).DFOf >= (mean(Dots.Vox(i).DFOf)-.1)));
% end

for i=1:Dots.Num
    for v=1:size(Dots.Vox(i).Pos,1)
        Dots.Vox(i).DF(v) = double(Dots.Vox(i).Igm(v))-mPred(Dots.Vox(i).Irm(v));
        Dots.Vox(i).DFOf(v) = Dots.Vox(i).DF(v) / mPred(Dots.Vox(i).Irm(v));
    end
    Dots.DF(i)=mean(Dots.Vox(i).DF);
    Dots.DFOf(i)=mean(Dots.Vox(i).DFOf);
    Dots.DFOfTopHalf(i)=mean(Dots.Vox(i).DFOf(Dots.Vox(i).DFOf >= (mean(Dots.Vox(i).DFOf)-.1)));
end

save([TPN 'data/Dots.mat'],'Dots')

%% Draw new filled
%clear IDs
%{
DeltaFill=zeros(Dots.ImSize,'uint8');
for i=1:Dots.Num
    for v=1 :size(Dots.Vox(i).DFOf,2)
        DeltaFill(Dots.Vox(i).Ind)=Dots.DFOf(i);
    end
end
image(max(DeltaFill,[],3)*10),pause(.1)

%imwriteNp(TPN,DeltaFill,'DeltaFill')
%clear DeltaFill
%}

%% Finish

TotalHours=toc/60/60
[TPN(size(TPN,2)-6:size(TPN,2)-1)]
RatioedAt=uint16(clock)
save([TPN 'data/RatioedAt.mat'],'RatioedAt')

%clear all
'Done Ratioing'


