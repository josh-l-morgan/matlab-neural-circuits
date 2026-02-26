

%% Draw figure

TPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\images\'
mkdir(TPN)


load([MPN 'dsObj.mat'])
load([MPN 'obI.mat'])
fsize = findObSize(dsObj);


%% View more cells
Dim = 1
cellList = ([108  129	109	117	162	131	116	137	130	135	106]);

showCellNames = cellList;

colMap = hsv(256);
col = colMap(ceil((1:length(showCellNames))*256/length(showCellNames)),:);
col = col(randperm(size(col,1)),:);
[I_partCell maxHeight] = showCellTop(obI,dsObj,showCellNames,col,Dim,fsize);
image(uint8(1.05.^double(I_partCell)))
image(I_partCell)



%% View more cells
cellList = ([1009]);

cellList = obI.cell.name;
cellList = cellList(cellList>=1000);

cellList = [ 1006  1009  1012  1014  1021 1025 1027 1028  1029  1030 1031 1032 1033  1034   1036  1037 ...
 1041 1050 1051 1053 1054 1055 1056]
bigPre108List = [ 1006        1009        1012        1014        1021        1025        1027        1028 ...
        1029        1030        1031        1032        1033        1034        1036        1037 ...
 1041        1050        1051        1053        1054        1055        1056]



Istack = zeros(size(I_partCell,1),size(I_partCell,2),3,length(cellList));

for i = 1:length(cellList)
    showCellNames = cellList(i)
    [I_ax axHeight] = showCellTop(obI,dsObj,showCellNames,[1 1 1],Dim,fsize);
    useI = axHeight>maxHeight;
            for c = 1:3
                Itemp = I_partCell(:,:,c);
                Itemp(useI) = I_ax(useI)*2;
                combineI(:,:,c) = Itemp;
            end
            combineI = combineI + I_ax;
%     combineI = (I_ax*2 + I_partCell * .8);
%     meanI = (I_ax+I_partCell)/2;
%     combineI(I_partCell>0) = meanI(I_partCell>0);
     image(uint8(combineI))
    imageName = sprintf('net108_ax%d.png',showCellNames(1));
    imwrite(combineI,[TPN imageName])
    Istack(:,:,:,i) = combineI;
    pause(.01)
end


%  useI = newHeight>maxHeight;
%             maxHeight(useI) = newHeight(useI);
%             
%             for c = 1:3
%                 Itemp = Ic(:,:,c);
%                 Itemp(useI) = I(useI)*col(i,c);
%                 Ic(:,:,c) = Itemp;
%             end
%             
% 

%%
Thicken = 10;
thickStack = zeros(size(Istack,1),size(Istack,2),3,size(Istack,4)*Thicken,'uint8');
for i = 1:size(Istack,4)
    Itemp = uint8(Istack(:,:,:,i));
    image(Itemp),pause(.01)
    for t = 1:Thicken
        thickStack(:,:,:,(i-1)*Thicken+t) = Itemp;
    end
end

myVideo = VideoWriter([TPN 'net108_ax3.avi']);
myVideo.FrameRate = 8;  % Default 30
myVideo.Quality = 75;    % Default 75
open(myVideo);
writeVideo(myVideo, uint8(thickStack));
close(myVideo);




