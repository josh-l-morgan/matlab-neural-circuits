
%{
%START PARALLEL
distcomp.feature( 'LocalUseMpiexec', false )
matlabpool open local 4
%}

%%
SPN = 'Z:\joshm\LGNs1\rawMontages\';
WaferList = {'w030','w031','w032','w033','w034','w035','w036'}
TPN = 'E:\SEM_Users\joshm\LGNs1\Processed\downSampWithReplace230\';



if ~exist(TPN,'dir');
    mkdir(TPN);
end

dSPN = dir(SPN); dSPN = dSPN(3:end);
dNams = cat(2,{dSPN.name});
isDirectory = cat(1,dSPN.isdir);
dNams = dNams(isDirectory);

getTileName = 'Tile_r3-c2';

montageDir = regexp(dNams,'Montage');
montageDir = ~cellfun('isempty',montageDir);

sortNams = {};
for w = 1:length(WaferList)
    
    waferName = WaferList{w};
    rightWafer = regexp(dNams,waferName);
    rightWafer = ~cellfun('isempty',rightWafer);
    rightSections =dNams(rightWafer& montageDir);
    sortSections = sort(rightSections);
    wafSecs{w} = sortSections;
    
end



shift1 = zeros(8,25600*25600/8);
shift2 = zeros(3200,8,length(sumAll)/3200/8);
shift2a = zeros(3200,8,length(sumAll)/3200/8);


c = 0;
for w = 1:length(wafSecs);
    secNam = wafSecs{w};
    c = c+1;
    
    for s = 1:length(secNam);
        disp(sprintf('Downsampling wafer %d of %d section %d of %d',...
            w,length(wafSecs),s,length(secNam)))
        %function[] = downSampWithReplacement(SPN,TPN,secNam,s,getTileName,c);
        
        nam = secNam{s};
        secDir  = dir([SPN nam]); secDir = secDir(3:end);
        tileNams = cat(1,{secDir.name});
        tileNams = sort(tileNams);
        findTile = regexp(tileNams,getTileName);
        findTile = ~cellfun('isempty',findTile);
        findTif = regexp(tileNams,'.tif');
        findTif = ~cellfun('isempty',findTif);
        targTif = find(findTile&findTif,1,'first');
        if ~isempty(targTif)
            foundName = tileNams{targTif};
            newName = [TPN 'ds' zeroBuf(c,5) '_' foundName];
            
            if exist([SPN nam '\' foundName],'file') & ~exist(newName)
                
                tic
                I = imread([SPN nam '\' foundName]);
                toc
                tic
                replaceWith = mean(I(:));
                replaceWith = mean([replaceWith tooHigh]);
                shift1(:) = I(:);
                tooHigh = 230;
                highI = shift1>=tooHigh;
                shift1(highI) = 0;%replaceWith;
                numHigh = sum(highI,1);
                sumAll = sum(double(shift1),1);
                testShift = zeros(3200,25600);
                testShift(:) = sumAll(:);
                % image(testShift/8);
                
                %shift2a = shift2;
                shift2(:) = sumAll(:);
                shift2a(:) = numHigh(:);
                finalHigh = squeeze(sum(shift2a,2));
                finalSum = squeeze(sum(shift2,2));
                
                dsMeansSc = finalSum./(64);
                dsMeansSc = finalSum./(64-finalHigh/2);
                dsMeansSc(finalHigh>=40) = replaceWith;
                image(256-dsMeansSc),pause(.01)
                toc
                imwrite(uint8(dsMeansSc),newName,'Compression','none')
            end
        end
        
        
        
            end
    end
    
    
    matlabpool close