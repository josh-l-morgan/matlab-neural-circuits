function[]=anaNN()
%finds nearest neigbors in xyz space

global DPN DFN TPN

%%load ratioed data
load([TPN 'data\DotStats.mat'])

%%Select dots where DFOf is greater then .5
Dots=DotStats(DotStats(:,3,1)>.5,:,3);

if size(Dots,1) %if there are Dots

%%Run all dots
for i = 1:size(Dots,1)
    Dist=sqrt((Dots(:,1)-(Dots(i,1))).^2 + (Dots(:,2)-(Dots(i,2))).^2 +(Dots(:,3)-(Dots(i,3))).^2);
    NN(i)=min(Dist(Dist~=0)); %%find Nearest neighbor
end %end i , running all dots

%%save Data

    
    save([TPN 'data\NN.mat'],'NN')
    if exist([TPN 'data\Results.mat'])
        load([TPN 'data\Results.mat'])
    end
    Results.CellStats.NN=NN;
    save([TPN 'data\Results.mat'],'Results')


end %if NN exist