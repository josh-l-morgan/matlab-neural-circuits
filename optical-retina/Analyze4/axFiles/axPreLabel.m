clear;clc;

%%Get directory with image folders for reading image
TPN=GetMyDir;
Slash=find(TPN=='\');
GPN=TPN(1:Slash(length(Slash)-1));
Ifolder=TPN(Slash(length(Slash)-1)+1:Slash(length(Slash))-1)


%Get directory and filename for saving color combine image


'reading',pause(.01)
I=oifread([GPN Ifolder '\']);


'writing',pause(.01)
[I R]=axJuggleCh(I);


% write axon channel for segmentation
Iname=[Ifolder(1:find(Ifolder=='.',1)-1) '_R'];
Iwrite([GPN Iname],R)

