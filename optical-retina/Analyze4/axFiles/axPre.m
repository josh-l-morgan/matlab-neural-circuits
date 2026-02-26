

%%Get directory with image folders for reading image
TPN=GetMyDir;
Slash=find(TPN=='\');
GPN=TPN(1:Slash(length(Slash)-1));
Ifolder=TPN(Slash(length(Slash)-1)+1:Slash(length(Slash))-1);


%Get directory and filename for saving color combine image


'reading',pause(.01)
I=oifread([GPN Ifolder '\']);


'writing',pause(.01)
[I R]=axJuggleCh(I);

clear TPN
TPN=GetMyDir;
I = axMask(TPN, I);


% write color combine image used by Matlab Ax scripts
CCfolder = ['\\Wongraid2\wonglab\Daniel\Pre_Data\Morphometry\'...
    Ifolder(1:find(Ifolder=='.',1)-1) '\I'];

if ~exist(CCfolder, 'dir')
    mkdir(CCfolder)
    axIwrite(CCfolder,I)
else
    [fileName pathName] = uiputfile;
    Iwrite([pathName fileName],I)
end




