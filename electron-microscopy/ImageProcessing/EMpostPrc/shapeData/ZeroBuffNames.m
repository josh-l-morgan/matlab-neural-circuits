bufNum = 3;

TPN = GetMyDir;
TPNd = dir(TPN); TPNd = TPNd(3:end);

for i = 1:length(TPNd);
    nam = TPNd(i).name;
    secPos = regexp(nam,'_Sec');
    if (secPos<(length(nam-4)));
        beforeNum = nam(1:secPos(1)+3);
        afterSec = nam(secPos(1)+4:end);
        
        unPos = regexp(afterSec,'_');
        numNam = afterSec(1:unPos(1)-1);
        afterNum = afterSec(unPos(1):end);
        newNum = repmat('0',[1 bufNum]);
        newNum(bufNum - length(numNam)+1:bufNum) = numNam;
        newNam = [beforeNum newNum afterNum];
        system(['rename ' [TPN nam] ' ' newNam]);
        eval(['!rename ' [TPN nam] ' ' newNam]);

    end
end


