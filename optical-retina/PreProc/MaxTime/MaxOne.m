

%%Get directory with image folders
TPN=GetMyDir;
Slash=find(TPN=='\');
GPN=TPN(1:Slash(length(Slash)-1));
Ifolder=TPN(Slash(length(Slash)-1)+1:Slash(length(Slash))-1);

clear Imaxs


'reading',pause(.01)
I=oifread([GPN Ifolder '\']);


'filtering',pause(.01)
   I=MyMedian(I,2); 
Iname=Ifolder(1:find(Ifolder=='.',1)-1);

'writing',pause(.01)
I=juggleCh(I);
Iwrite([GPN Iname],I)

%% collect maxes
Imax=max(I,[],4);


%%Write Max
imwrite(Imax,[GPN Iname '_max.tif'],'Compression','none')




