%euclideanator
function[arrayOut] = euclideanate(inputObj)
arrayOut=[];
for vastSub=1:length(inputObj)
    subArray=struct2array(inputObj(vastSub));
    if ~isempty(subArray)
        subArrayID=[repmat(vastSub,length(subArray(:,1)),1) subArray];
        if isempty(arrayOut)
            arrayOut=subArrayID;
        else
            arrayOut=[arrayOut;subArrayID];
        end
    end
end