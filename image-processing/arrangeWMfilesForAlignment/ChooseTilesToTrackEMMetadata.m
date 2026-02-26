



SPN = 'Y:\Active\Merlin\Morgan Lab\MasterRaw\KxS_L_Albino_LGN\MedRes\'

TPN = 'Z:\Active\morganLab\DATA\KxS_L_AlbinoLGN\TrackEM2_projects\Kxs_MedR2C12\'
TPN = 'Z:\\Active\\morganLab\\DATA\\KxS_L_AlbinoLGN\\TrackEM2_projects\\Kxs_MedR2C12\\'

metaFileName = 'Kxs_AllTiles_10Sections.txt';
fileName = [TPN metaFileName];


useC = [1:3]; %choose column
useR = [1:3]; %choose row
rawTileSize = 5120;
overlap = 6;
useSections = [300 : 310];
tileSize = round(rawTileSize * (100-overlap)/100);


dSPN = dir(SPN)
tc = 0;
clear tPath r c d w s l
for i = 1:length(dSPN)
    nam1 = dSPN(i).name;
    disp(sprintf('reading folder %s',nam1))
    if length(nam1)>2 & dSPN(i).isdir
        if strcmp(nam1(1),'w')
            dNam1 = dir([SPN nam1 '\']);
            for m = 1:length(dNam1)
                nam2 = dNam1(m).name;
                if length(nam2)>2 & dNam1(m).isdir
                    if nam2(1:2) == nam1(1:2)
                        dNam2 = dir([SPN nam1 '\' nam2 '\*.tif']);
                        for t = 1:length(dNam2)
                            nam3 = dNam2(t).name;
                            if strcmp(nam3(1:5),'Tile_')
                                r1 = regexp(nam3,'_r');
                                r2 = regexp(nam3,'-c');
                                r3 = regexp(nam3,'_w');
                                r4 = regexp(nam3,'_sec');
                                secNum = str2num(nam3(r4(1)+4:end-4));
                                if ~isempty(secNum)
                                    tc = tc + 1;
                                    tPath{tc} = [SPN nam1 '\' nam2 '\' nam3];
                                    r(tc) = str2num(nam3(r1(1)+2:r2(1)-1));
                                    c(tc) = str2num(nam3(r2(1)+2:r3(1)-1));
                                    w(tc) = str2num(nam3(r3(1)+2:r3(1)+3));
                                    %w(tc) = str2num(nam3(r3(1)+2:r4(1)-1));
                                    s(tc) = str2num(nam3(r4(1)+4:end-4));
                                    l(tc) = length(nam3);
                                    d(tc) = dNam2(t).datenum;
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end


%% make unique tile ids
maxR = max(r);
maxC = max(c);
maxW = max(w);
maxS = max(s);
%tileID = sub2ind([maxW maxS maxR maxC],w,s,r,c);
tileID = sub2ind([maxC maxR maxS maxW],c,r,s,w);

useIds = sort(unique(tileID),'ascend');


clear uR uC uW uS uPath
for i = 1:length(useIds)

    isId = find(tileID==useIds(i));
    if length(isId)>1
        iDnum = d(isId);
        newest = find(iDnum == max(iDnum),1);
        isId = isId(newest);
    end
    uR(i) = r(isId);
    uC(i) = c(isId);
    uW(i) = w(isId);
    uS(i) = s(isId);
    uPath{i} = tPath{isId};
end

%% extract

 minX = (min(useC)-1) * tileSize;
 minY = (min(useR) - 1) * tileSize;

%%Code positions in X, Y, Z, starting at 0
tc = 0;
if exist([TPN metaFileName],'file'),delete([TPN metaFileName]);end
fid = fopen(fileName,'w');


for w1 = min(w):maxW
    sections = unique(s(w==w1));

    for s1 = sections
        %disp(sprintf('copying w%d of %d, s%d',w1,maxW,s1))
        tc = tc+1;

        if sum(useSections==tc) | isempty(useSections) %check if section should be used
        for cC = 1:length(useC)
            for cR  = 1:length(useR)
                nextTile = find( (uW == w1) & (uS == s1) & (uR == useR(cR)) & (uC ==useC(cC)));
                if ~isempty(nextTile)

                    %newName = sprintf('%07.0f.tif',tc);
                    %copyfile(uPath{nextTile},[TPN newName]);

                    x = (useC(cC)-1) * tileSize - minX;
                    y = (useR(cR) - 1) * tileSize - minY;

                    lineText = sprintf('%s,%d,%d,%d\n',uPath{nextTile},x,y,tc-1);

                    fwrite(fid,lineText);

                end
            end

        end
        end
    end
end
fclose(fid)











