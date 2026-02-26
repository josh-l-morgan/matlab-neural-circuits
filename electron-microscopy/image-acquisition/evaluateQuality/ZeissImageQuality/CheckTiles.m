%% test quality checks
%{

compare against
1) tissue statistics
2) Dwell time nonlinearity - 200 ns dip
3) digital histogram

jm - current map
%}


clear all
SPN = GetMyDir;
dSPN = dir(SPN);  dSPN = dSPN(3:end);

nams = {};
for i = 1:length(dSPN)
    nam = dSPN(i).name
    if sum(regexp(nam,'.tif'))
        if sum(regexp(nam,'Tile'))
            nams{length(nams)+1} = nam;
            rs = regexp(nam,'_r');
            cs = regexp(nam,'-c');
            row(length(nams)) = str2num(nam(rs+2));
            col(length(nams)) = str2num(nam(cs+2));
        end
    end
end
clear qualResults
for i = 1:length(nams);
    nam = nams{i};
    
   iNam = [SPN nam];
   
   startTime = datenum(clock);
   [quality qImage ] = checkFileQualStandalone(iNam);
   stopTime = datenum(clock);
   myQualTime = (stopTime-startTime)* 24 * 60;
   
   startTime = datenum(clock);
   
   Ifo = imfinfo(iNam);
   Width=Ifo.Width;
   Height = Ifo.Height;
   pRange = [round(Width/2) - 499  round(Width/2) + 500];
   
   
   
   I = imread(iNam,'PixelRegion',{ pRange,pRange});
   %I = imread(iNam);
   I = mean(I,3);
   
   
   zqual = zeissQual((uint8(I/5))*5);
   
   stopTime = datenum(clock);
   zeisTime = (stopTime-startTime)* 24 * 60;
   
   zqual
   qualResults(i).quality = quality.quality;
   qualResults(i).myTime = myQualTime;
   qualResults(i).snr = zqual.snr;
   qualResults(i).resolution = zqual.resolution;
   qualResults(i).zeisTime = zeisTime;
   qualResults(i).name = nam;
   
end

%% Display Results
resultsStr = [];
qMat = zeros(max(row),max(col));
numMat = qMat;
for i = 1:length(qualResults)
    resultsStr = cat(2,resultsStr,...
        sprintf('%s           %3.1f        %3.1f       %1.3f \n',...
    qualResults(i).name, qualResults(i).quality,...
        qualResults(i).snr, qualResults(i).resolution));
    
    
    qMat(row(i),col(i)) = qMat(row(i),col(i)) + qualResults(i).quality;
    numMat(row(i),col(i)) = numMat(row(i),col(i)) + 1;
end

qMat = round(qMat./numMat)
%%
myQual = cat(1,qualResults.quality);
myQual = myQual/max(myQual);
resolution = cat(1,qualResults.resolution);
resolution = resolution/max(resolution);
snr = cat(1,qualResults.snr);
snr = snr/max(snr);

subplot(2,1,1)
plot(myQual)
hold on
plot(resolution,'g')
plot(snr,'r')
hold off
subplot(2,1,2)
scatter(myQual,resolution)
hold on
scatter(myQual,snr)
hold off

%resultsStr




