%
% clear all
[TFN TPN] = GetMyFile;
LPN = GetMyDir;

[xnum, xtext,xraw] = xlsread([TPN TFN]);
%%
dates = [];
dType= [];
for i = 1:size(xraw,1)
    dates(i) = datenum(xraw{i,1})+xraw{i,2};
    switch xraw{i,5}
        case 'eht'
            dType(i) = 1;'eht'
        case 'current'
            dType(i) = 2;'current'
        case 'align'
            dType(i) = 3;'align'
        case 'check'
            dType(i) = 4;'check'
        case 'gun'
            dType(i) = 5;'gun'
        case 'guncheck'
            dType(i) = 6;'guncheck'
    end
    
end
eht = find(dType==1);
current = find(dType==2);
align = find(dType==3);
check = find(dType== 4);
gun = find(dType== 5);
guncheck = find(dType== 6);

%% logs
dLPN = dir(LPN); dLPN = dLPN(3:end);
allTimes = [];
allQual = [];
for i = 1:length(dLPN)
    nam = dLPN(i).name
    if sum(regexp(nam,'LogBook_'))
        load([LPN nam]);
        
        regexp(nam,'_w')
        
        waferNumber(i) = str2num(nam(11:13));
        qualNames = logBook.sheets.quality.data(:,1);
        qualVals = [logBook.sheets.quality.data{:,3}];
        
        
        iCon.Names = cat(1,logBook.sheets.imageConditions.data(:,1)); %image names
        iCon.Vals = cat(1,logBook.sheets.imageConditions.data{:,5})*1000000; %working distance
        iCon.times = cat(1,logBook.sheets.imageConditions.data(:,32));
        
        nameList = iCon.Names;
        clear qual2iCon
        for n = 1:length(qualNames)
            tNam = qualNames{n};
            if ~mod(n,100), disp(tNam),end
            matched = regexp(nameList,tNam);
            matches = find(~cellfun(@isempty,matched));
            
            for m = 1:length(nameList)
                otherName = nameList{m};
                if strcmp(otherName,tNam)
                    qual2iCon(n) = m;
                    nameList{m} = '';
                    break
                end
            end
        end
        qualTimes = datenum(iCon.times(qual2iCon));
        qualTimes = datenum(iCon.times);
        waferTimes(i) = min(qualTimes);
        allTimes = cat(1,allTimes, qualTimes);
        allQual = cat(1,allQual, qualVals');
    end
end
%%

[sortTimes idT] = sort(allTimes);
sortQual = allQual(idT);
normQual = (sortQual-250)/100;

binWidth = 1/8;

minTime = min(sortTimes);
maxTime = max(sortTimes);
binTimes = [minTime:binWidth:maxTime];
c = 0;
clear meanQual meanTime
for i = 1:length(binTimes)
   useTimes = find((sortTimes > binTimes(i)) & (sortTimes<(binTimes(i)+binWidth))); 
   if ~isempty(useTimes)
      c = c+1;
      meanQual(c) = mean(normQual(useTimes));
      minQual(c) = min(normQual(useTimes));
      meanTime(c) = mean(sortTimes(useTimes));
   end
end



%% Display

plot(meanTime,meanQual,'k')
hold on
scatter(meanTime,meanQual,'k','.')
%scatter(meanTime,minQual,'m','.')

scatter(dates(eht),dates(eht)./dates(eht)-1,'r','x','LineWidth',10)

alignDates = dates(align);
xAlign = [xraw{align,3}];
meanX = mean(xAlign);
xAlign = xAlign - meanX;
yAlign = [xraw{align,4}];
meanY = mean(yAlign);
yAlign = yAlign - meanY;
[alignDates idx] = sort(alignDates);
xAlign = xAlign(idx);
yAlign = yAlign(idx);



checkDates = dates(check);
xcheck = [xraw{check,3}];
xcheck = xcheck - meanX;
ycheck = [xraw{check,4}];
ycheck = ycheck - meanY;
[checkDates idx] = sort(checkDates);
xcheck = xcheck(idx);
ycheck = ycheck(idx);

plot(alignDates,xAlign,'g')
scatter(alignDates,xAlign,'g','.')
plot(alignDates,yAlign,'b')
scatter(alignDates,yAlign,'b','.')
ylim([-3 3])

scatter(checkDates,xcheck,'g','o')
scatter(checkDates,ycheck,'b','o')

scatter(waferTimes,waferTimes./waferTimes-1,'m','o','LineWidth',2)


hold off


