



SPN = 'Y:\Active\Merlin\Morgan Lab\MasterRaw\KxS_L_Albino_LGN\MedRes\'

TPN = 'Z:\Active\scalableminds\KxS_L_Albino_LGN\subStackMiddle\'



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
                                secNum = str2num(nam3(r4(1)+4:end-4));
                                if ~isempty(secNum)

                                    tc = tc + 1;
                                    r1 = regexp(nam3,'_r');
                                    r2 = regexp(nam3,'-c');
                                    r3 = regexp(nam3,'_w');
                                    r4 = regexp(nam3,'_sec');

                                    tPath{tc} = [SPN nam1 '\' nam2 '\' nam3];
                                    r(tc) = str2num(nam3(r1(1)+2:r2(1)-1));
                                    c(tc) = str2num(nam3(r2(1)+2:r3(1)-1));
                                    w(tc) = str2num(nam3(r3(1)+2:r4(1)-1));
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
tileID = sub2ind([maxW maxS maxR maxC],w,s,r,c);

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


useR = 2;
useC = 2;
tc = 0;
for w1 = 1:maxW

    for s1 = 1:maxS
        disp(sprintf('copying w%d of %d, s%d',w1,maxW,s1))
        
        nextTile = find( (uW == w1) & (uS == s1) & ( uR == useR) & (uC == useC));
        if ~isempty(nextTile)
            tc = tc+1;
            newName = sprintf('%07.0f.tif',tc);
            copyfile(uPath{nextTile},[TPN newName]);

        end
    end
end











