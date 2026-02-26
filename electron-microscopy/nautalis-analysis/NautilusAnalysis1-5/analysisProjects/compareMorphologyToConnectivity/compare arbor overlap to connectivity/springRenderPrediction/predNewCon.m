%%Sorts an axon's connectivity according to which dendritic arbor it
%%overlapped the most
function[newCon] = predNewCon(rawCon,rawPredSyn,mask)

if ~exist('mask','var')
    mask = (con*0+1)>0;
else
    mask = mask>0;
end
con = rawCon.*mask;
predSyn = rawPredSyn.*mask;

for i = 1:size(con,1)
    
   conLine = con(i,:);
   predLine = predSyn(i,:);
   
   [sortPred predIDX] = sort(predLine,'descend');
   [sortCon conIDX] = sort(conLine,'descend');
   
   newConLine(predIDX) = conLine(conIDX);
   
   newCon(i,:) = newConLine; 
    
end
