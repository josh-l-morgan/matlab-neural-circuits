
%%Define which channels should be used to make images {[1 2 3] [4 5 6]} is
%%default and makes two RGB images side by side with channels in order
stkChan = {[4 3 2]};


%SPN = 'H:\GadGFP_Albino\'
SPN = 'Z:\Active\morganLab\DATA\Albino_Islands\Immuno\LxS\'
SPN = 'Z:\Active\morganLab\DATA\Ablation\LxQ\'
SPN = 'Z:\Active\morganLab\DATA\Anole\AxJ\Lizard6\'
SPN = 'Z:\Active\morganLab\DATA\Anole\DiI_Tests\'
SPN = 'Z:\Active\morganLab\DATA\Anole\AxJ\AnxJ_Anole6\RetinaSeriesWide\'
SPN = 'H:\FV1000\AxW_retina\'
SPN = 'H:\FV1000\AxU\'
SPN = 'H:\FV1000\AxV\Anole3_3colDye+Cocktail\'
SPN = 'H:\FV1000\AxV\Anole5_5-14-2025\'
SPN = 'H:\FV1000\AxV\Anole5_reslice\'
SPN = 'H:\FV1000\AxV\AnoleX\'
SPN = 'H:\FV1000\AxV\'
%% Find oifs

disp('Finding Oifs')
oList = {};
odList = {};

dirNow = SPN;
od = dir([dirNow '*.oif']);
odList = cat(1,odList,od.folder);
oList = cat(1,oList,od.name);
d1 = dir(SPN);
for f1 = 3:  length(d1)
    if d1(f1).isdir
        dirNow = [d1(f1).folder '\' d1(f1).name '\'];
        od = dir([dirNow '*.oif']);
        oList = cat(1,oList,od.name);
        odList = cat(1,odList,od.folder);
        d2 = dir(dirNow);
        for f2 = 3:length(d2)
            dirNow = [d2(f2).folder '\' d2(f2).name '\'];
            od = dir([dirNow '*.oif']);
            odList = cat(1,odList,od.folder);
            oList = cat(1,oList,od.name);
            d3 = dir(dirNow);
            for f3 = 3:length(d3)
                dirNow = [d3(f3).folder '\' d3(f3).name '\'];
                od = dir([dirNow '*.oif']);
                odList = cat(1,odList,od.folder);
                oList = cat(1,oList,od.name);
                d4 = dir(dirNow);
                for f4 = 3:length(d4)
                    dirNow = [d4(f4).folder '\' d4(f4).name '\'];
                    od = dir([dirNow '*.oif']);
                    odList = cat(1,odList,od.folder);
                    oList = cat(1,oList,od.name);
                    d5 = dir(dirNow);
                    for f5 = 3:length(d5)
                        dirNow = [d5(f5).folder '\' d5(f5).name '\'];
                        od = dir([dirNow '*.oif']);
                        odList = cat(1,odList,od.folder);
                        oList = cat(1,oList,od.name);
                        d6 = dir(dirNow);
                        for f6 = 3:length(d6)
                            dirNow = [d6(f6).folder '\' d6(f6).name '\'];
                            od = dir([dirNow '*.oif']);
                            odList = cat(1,odList,od.folder);
                            oList = cat(1,oList,od.name);
                            d7 = dir(dirNow);
                        end
                    end
                end
            end
        end
    end
end


%% run maxes
for d = 1:length(oList)
    try
    disp(sprintf('Running %d of %d, %s%s',d,length(oList),odList{d},oList{d}))
    
    CPN = [odList{d} '\' oList{d} '.files\'];
    try RGBs = returnMaxFromOIF(CPN);
    catch
        RGBs = {};
    end
    if ~isempty(RGBs)
        rgbCat = [];
        for c = 1:length(RGBs)
            rgbCat = cat(2,rgbCat,RGBs{c});
        end


        fileName = sprintf('%s\\%s_three.png',odList{d},oList{d});
        imwrite(rgbCat,fileName,'png');


        after = odList{d}(length(SPN)+1:end);
        slashes = regexp(after,'\');
        after(slashes) = '-';

        fileName = sprintf('%s\\%s_max.png',SPN,after);
        imwrite(rgbCat,fileName,'png')
    end
    end
end


