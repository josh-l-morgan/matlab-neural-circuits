
SPN = GetMyDir
TPN = GetMyDir
dSPN = dir(SPN),dSPN = dSPN(3:end);
%%

for i = 1:length(dSPN)
    
    nam = dSPN(i).name;
    I = 255-double(imread([SPN nam]));
    Ival = sort(I(:));
    minI = Ival(1);
    maxI = Ival(end);
    midI = median(I(:));%Ival(round(length(Ival)/2))
    
    Ifix = 255 - (I-midI+50)*205/(maxI-midI);
    image(Ifix),pause(.1)
    
    filename = [TPN nam];
    imwrite(Ifix,filename);
  
end
