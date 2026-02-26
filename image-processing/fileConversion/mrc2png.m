
fileName = 'E:\other\OPN\3view-stack-final-bin1.mrc'
TPN = 'E:\other\OPN\pngs\'
mkdir(TPN)

addpath('C:\Users\jlmorgan\Documents\MATLAB\otherPeoplesCode\EMIODist2\')

sliceNum = 1000
for i = 1:sliceNum

   I = ReadMRC(fileName,i,1); 
   newName = sprintf('%s%05.0f.png',TPN,i);
   imwrite(I,newName);
    
end













