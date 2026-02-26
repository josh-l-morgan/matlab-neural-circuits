
%{
%START PARALLEL
distcomp.feature( 'LocalUseMpiexec', false )
matlabpool open local 4
%}

%%
SPN = 'Z:\joshm\LGNs1\rawMontages\';
WaferList = {'w030','w031','w032','w033','w034','w035','w036'}
TPN = 'E:\SEM_Users\joshm\LGNs1\Processed\downSampWithReplace\';



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

c = 0;
for w = 1:length(wafSecs);
    secNam = wafSecs{w};
    c = c+1;
    
    parfor s = 1:length(secNam);
        disp(sprintf('Downsampling wafer %d of %d section %d of %d',...
            w,length(wafSecs),s,length(secNam)))
        downSampWithReplacement(SPN,TPN,secNam,s,getTileName,c);
        
       
    end
end


matlabpool close