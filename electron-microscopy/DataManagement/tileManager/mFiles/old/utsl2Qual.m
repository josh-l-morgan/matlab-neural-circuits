%%Compare qualities of multiply taken tiles to identify best.




load('..\matFiles\utslTiles.mat')


%% Make list based on logBooks
%if ~exist([SPN 'BadRetakeFromLogBookList.mat'],'file')

if exist('..\matFiles\logQuals.mat','file')
    load('..\matFiles\logQuals.mat');
else
    
    noLogFound = {};
    logQuals.date = datestr(clock);
    for w = 1:length(utslTiles.waf)
        waf = utslTiles.waf(w).name
        logName = ['LogBook_' waf '.mat'];
        
        if ~exist(['..\logBooks\' logName],'file')
            noLogFound{length(noLogFound)+1} = logName;
            logQuals.waf(w).foundBook = 0;
        else  %read log
            logQuals.waf(w).foundBook = 1;
            
            load(['..\logBooks\' logName])
            tileFiles = cat(1,logBook.sheets.quality.data(:,1));
            qualities = cat(1,logBook.sheets.quality.data{:,3});
            
            for i = 1:length(tileFiles)
                tF = tileFiles{i};
                slashes = regexp(tF,'\');
                dots = regexp(tF,'.tif');
                tileNames{i} = tF(slashes(end)+1:dots(end)-1);
            end
            
            
            for s = 1: length(utslTiles.waf(w).sec)
                for t = 1:length(utslTiles.waf(w).sec(s).tileNames)
                    
                    tileName = utslTiles.waf(w).sec(s).tileNames{t};
                    ind =find(~cellfun(@isempty,regexp(tileNames,tileName)));
                    qualVal = qualities(ind);
                    
                    logQuals.waf(w).sec(s).tile(t).tileName = tileName;
                    logQuals.waf(w).sec(s).tile(t).qualInd = ind;
                    logQuals.waf(w).sec(s).tile(t).qualVal = qualVal;
                    
                end % run all tiles
            end %run all sections
            
            
            
            
        end % if log file exists
        
        
    end
    
    logQuals.noLogFound = noLogFound';
    
    save('..\matFiles\logQuals.mat','logQuals')
    
end % if logQuals does not exist



%% Parse logQuals
worseThresh = .2;
minQuality = 150;

g = 0;
clear gotWorse
for w = 1:length(logQuals.waf)
    if logQuals.waf(w).foundBook
        for s = 1: length(logQuals.waf(w).sec)
            
            clear lastQuals bestQuals
            flagWorse = 0;
            for t = 1:length(logQuals.waf(w).sec(s).tile)
                tile = logQuals.waf(w).sec(s).tile(t);
                bestQual = max(tile.qualVal);
                lastQual = tile.qualVal(find(tile.qualInd == max(tile.qualInd),1));
                
                if ~isempty(lastQual)
                    if (lastQual + (worseThresh*lastQual)) < bestQual % check if image got significantly worse
                        flagWorse = 1;
                        g = g+1;
                        gotWorse.tileName{g,1} = tile.tileName
                        gotWorse.lastQual(g) = lastQual;
                        gotWorse.bestQual(g) = bestQual;
                    end
                    
                    lastQuals(t) = lastQual;
                    bestQuals(t) = bestQual;
                else
                    lastQuals(t) = 0;
                    bestQuals(t) = 0;
                    
                        
                end %is empty
            end %run tiles
            
            logQuals.waf(w).sec(s).lastQuals = lastQuals;
            logQuals.waf(w).sec(s).bestQuals = bestQuals;
            logQuals.waf(w).sec(s).passQuals = sum(bestQuals<minQuality)==0;
            logQuals.waf(w).sec(s).gotWorse = flagWorse;

            
        end %run sections
    end %if found book
end %run wafers

save('..\matFiles\logQuals.mat','logQuals')


%%
plot(gotWorse.lastQual)
hold on
plot(gotWorse.bestQual,'r')
hold off

%
%
%     logDir = ['..\logBooks']
%     dLogDir = dir(logDir); dLogDir = dLogDir(3:end);
%
%     logList = {};
%     for i = 1:length(dLogDir)
%         nam = dLogDir(i).name;
%         if sum(regexp(nam,'LogBook_'))
%            logList{length(logList)+1} = nam;
%         end
%     end
%
%         gotWorse = {};
%         badImages = {};
%     for i = 1:length(logList)
%         sprintf('running %s, wafer %d of %d',logList{i},i,length(logList))
%        load([logDir '\' logList{i}]);
%
%         tileNames = cat(1,logBook.sheets.quality.data(:,1));
%         qualities = cat(1,logBook.sheets.quality.data{:,3});
%
%             %%parse logbook to find failed final tiles
%
%             %%find all retakes
%             [uTiles ia ic] = unique(tileNames);
%             for t = 1:length(uTiles)
%                 quals = qualities(ic == t) ;
%                 if max(quals) > (quals(end)+quals(end)/20) %last image is not the best
%                     gotWorse{length(gotWorse)+1} = uTiles{t};
%                     worseQuals{length(gotWorse)} = quals;
%                     quals
%                 end
%                 if quals(end)<100
%                     badImages{length(badImages)+1} = uTiles{t};
%                 end
%             end
%
%
%     end  %run all logbook files
%
%     BadImagesFromLog.badImages = badImages;
%     BadImagesFromLog.gotWorse = gotWorse;
%     BadImagesFromLog.worseQuals = worseQuals;
%     save([SPN 'BadImagesFromLog.mat'],'BadImagesFromLog')
% %
% %
%     %save([SPN 'BadRetakeFromLogBookList.mat'],'gotWorse')
% %
% % else
% %     load([SPN 'BadRetakeFromLogBookList.mat'])
% % end
%
% %%
% for i = 1:length(gotWorse)
%
% end
%
%
% %%
% %{
% missingTile = {};
% for t = 1:length(gotWorse)
%     tNam = gotWorse{t};
%     slashes = regexp(tNam,'\');
%     secNam = tNam(slashes(2)+1:slashes(3)-1);
%     tileNam = tNam(slashes(3)+1:end-4);
%     waf = str2num(secNam(2:4));
%     if sum(find(allWafs == waf))
%
%
%         if exist([SPN secNam],'dir')
%             dSec = dir([SPN secNam])
%             tilesToCheck = {};
%             for f = 1:length(dSec)
%                 nam = dSec(f).name;
%                 if sum(regexp(nam,tileNam)) & sum(regexp(nam,'.tif'))
%                     tilesToCheck{length(tilesToCheck)+1} = nam;
%                 end
%             end
%
%             if ~isempty(tilesToCheck)
%                 qualChecked = [];
%                 for tc = 1:length(tilesToCheck)
%                     qualVal = checkFileQual([SPN secNam '\' tilesToCheck{tc}])
%                     qualChecked(tc) = qualVal.quality;
%                 end
%                 bestQual{t} = tilesToCheck{find(qualChecked == max(qualChecked),1,'last')}
%
%
%             else
%                 missingTile{length(missingTile)+1,1} = tileNam;
%             end
%
%
%         else
%
%
%
%         end
%     end
% end
%
% save([SPN 'BadRetakeFromQualityCheck.mat'],'bestQual')
% %}