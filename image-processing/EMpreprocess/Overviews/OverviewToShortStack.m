

SPN = 'E:\Pratyush\masterUTSL\'

TPN = 'E:\Pratyush\EMprocStacks\Overview1\'
mkdir([TPN 'FullStack\'])
mkdir([TPN 'WaferAve\'])

dSPN = dir([SPN 'w*'])
wNam = {dSPN.name}
clear wVal
for i = 1:length(wNam)
    wVal(i) = str2num(wNam{i}(2:end));
end

[sortWaf ia] = sort(wVal,'ascend');
wNam = wNam(ia);
c = 0;
for i = 1:length(wNam)
    
    foldNam = [SPN wNam{i} '\SectionOverviewsAlignedWithTemplateDirectory\'];
    odir = dir([foldNam '*.tif']);
    inams = {odir.name};
    clear tval
    for t = 1:length(inams)
        nam = inams{t};
        u = regexp(nam,'_');
        d = regexp(nam,'.tif');
        tval(t) = str2num(nam(u+1:d-1));
        
    end
    
    [sortSec ia] = sort(tval,'ascend');
    inams = inams(ia);
    
    info = imfinfo([foldNam inams{1}]);
    Isum = zeros(info.Height,info.Width,'double');
    for t = 1:length(inams)
        disp(sprintf('running section %d of wafer %d',t,i))
        
        I = imread([foldNam inams{t}]);
        c = c + 1;
        fileName = sprintf('%sFullStack\\%05.0f_w%s_s%s.tif',TPN,c,wNam{i},inams{t});
        imwrite(I,fileName);
        Isum = Isum + double(I);
    
    end
    
    Isum = Isum/length(inams);
    fileName = sprintf('%sWaferAve\\%05.0f.tif',TPN,i);
    imwrite(uint8(Isum),fileName);
    
end








