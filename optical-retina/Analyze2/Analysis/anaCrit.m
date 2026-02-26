%function[] = anaCrit(TPN, DPN)

%{
KPN=GetMyDir

Kdir=dir(KPN);
Kdir=Kdir(3:size(Kdir,1));

for k = 66:66%size(Kdir,1)
   TPN = [KPN '\' Kdir(k).name '\'];
    DPN = [TPN 'I\']        
%}
%TPN=GetMyDir
 
%% Load File
%{
DPN = GetMyDir
f=find(DPN=='\');
f2=f(size(f,2)-1);
f3=f(size(f,2)-2);
TPN=DPN(1:f2); %Define target folder (one level up from files)
%}
if exist([TPN 'Dots.mat'])
load([TPN 'Dots.mat'])


%% Define Criteria
Vol=6
Cut=2
ItSum=4
DFOf=1.5 % 2 looks good
DF=10
DistToMask=1

%% Select Puncta
clear pass
for i = 1: Dots.Num
   v=Dots.Vol(i)>=Vol;
   c=Dots.Cut(i)<Cut;
   itsum=Dots.ItSum(i)>=ItSum;
   dfof=Dots.DFOf(i)>=DFOf;
   df=Dots.DF(i)>=DF;
   dm=Dots.DistToMask(i)<=DistToMask;
   pass(i)=v  & c & itsum & dm & df & dfof ;
end

P=find(pass')'; %% list of passing puncta
Dots.OK=pass';
save([TPN 'Dots.mat'],'Dots')

%% Draw puncta matrix



%%ID puncta
%{
IDs=zeros(Dots.ImSize(1),Dots.ImSize(2),'uint16');
for i = 1 : size(P,1)
   for v=1:size(Dots.Vox(P(i)).Pos,1)
        IDs(Dots.Vox(P(i)).Pos(v,1),Dots.Vox(P(i)).Pos(v,2))=P(i); 
   end
end
imwrite(IDs,[TPN 'images\IDs.tif'],'Compression','none')
%}


Passed = zeros(Dots.ImSize,'uint8');
%{
for i = 1: size(Dots.Vox,2)
    Passed(Dots.Vox(i).Ind)=50;    
end
%}

for i = 1: size(P,2)
   Passed(Dots.Vox(P(i)).Ind)=200;     
end

colormap gray(255)
maxPassed=max(Passed,[],3);
image(maxPassed),pause(.01)

%imwriteNp(TPN,Passed,'Passed')
imwrite(maxPassed,[TPN 'images\maxPassed.tif'],'Compression','none')
%}



end %if Dots exist
%end %run all cells



 


