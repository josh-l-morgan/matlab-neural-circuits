SPN = 'E:\SeanMcCracken\Victoria\AMIRA\Tiff\'
load([SPN 'dat.mat'])

badData = [9]; %remove bsd data

%% Determine experimtent type from segmentation images present
expType = zeros(length(dat),1);
for i = 1:length(dat)

    if sum(dat(i).hasImage) == 4
        expType(i) = 1;
    elseif dat(i).hasImage(2)
        expType(i) = 2;
    elseif dat(i).hasImage(3)
        expType(i) = 3;
    else
        expType(i) = 0;
    end
end

%expType(9) = 100;%ditch Bad image

lgnId = zeros(length(dat),1);
for i = 1:length(dat)
    lgnId(i) = dat(i).f.lgnId;
end


%% Group by months post crush
%%Groups are 0, 3, 6
threeMonth = {'SM19-1' 'SM19-2'};
sixMonth = {'SM13-3' 'SM14-2' 'SM14-3' 'SM13-3' 'SM14-1' 'SM14-2' 'SM14-3'};
postCrush = zeros(length(dat),1)+1;
for i = 1:length(dat)
    e1 = dat(i).f.expNam;

    for m = 1:length(threeMonth);
        e2 = threeMonth{m};
        if strcmp(e1,e2)
            postCrush(i) = 2;
        end
    end
    for m = 1:length(sixMonth);
        e2 = sixMonth{m};
        if strcmp(e1,e2)
            postCrush(i) = 3;
        end
    end
end

[postCrush expType]

postCrush(badData) = 0;
%% Show experiment info and images
for i = 1:length(dat)
    clf
i
    imshow(dat(i).Isum)
    lgnStr = sprintf('%s expType = %d, lgn = %d',dat(i).f.expNam,expType(i),dat(i).f.lgnId);
    disp(lgnStr)
    drawnow
   
end


xOffset = (postCrush) *2 -1;
noCrush = find(postCrush==1);
isCrush = find(postCrush>1);
xOffset(noCrush) = xOffset(noCrush) + lgnId(noCrush)/2;
xOffset(postCrush>1) = xOffset(postCrush>1) + (expType(postCrush>1)-1)/2;



%% Generate structure for dat grouped by type of experiment
clear g
%% 
for t = 1:3;
    isT = find(postCrush==t);
    g.dat{t} = dat(isT);
    g.areaFCI{t} = cat(1,dat(isT).areaFCI);
    g.countFCI{t} = cat(1,dat(isT).countFCI);
    g.n(t) = length(isT);
    g.expType{t} = expType(isT);
    g.xOffset{t} = xOffset(isT);
    g.lgnId{t} = lgnId(isT);
end

%%Find paired experiments
for t = 1:3;
    g.expIds{t} = [];
    for i = 1:length(g.dat{t});
        g.expIds{t} = [g.expIds{t} g.dat{t}(i).f.expId];
    end
    g.uExpIds{t} = unique(g.expIds{t} );

end


for t = 1:3
    g.isPaired{t} = [];
    g.notPaired{t} = [];

    for i = 1:length(g.uExpIds{t});
        if (sum(g.expIds{t}==g.uExpIds{t}(i))==2)
            g.isPaired{t} = [g.isPaired{t} g.uExpIds{t}(i)];
        else
            g.notPaired{t} = [g.isPaired{t} g.uExpIds{t}(i)];
        end
    end
end
% 
% g.isPaired{2} = [];
% g.notPaired{2} = [];
% for i = 1:length(g.uExpIds{2});
%     if  (sum(g.expIds{3}==g.uExpIds{2}(i))==1);
%         g.isPaired{2} =[g.isPaired{2} g.uExpIds{2}(i)];
%     else
%         g.notPaired{2} =[g.notPaired{2} g.uExpIds{2}(i)];
%     end
% end
% g.notPaired{2} = setdiff(g.expIds{3},g.isPaired{2});
    
%% Show results
clf

regionName = {'Full' 'Contra' 'Ipsi'};
measureName = {'Area' 'Count' 'Density'};

clear showRes
showRes{1} = g.areaFCI;
showRes{2} = g.countFCI;
for t = 1:3
    showRes{3}{t} = g.countFCI{t}./g.areaFCI{t};
end



for sr = 1:length(showRes)

    res = showRes{sr};
    maxY = 0;
    for a = 1:3
        maxY = max([maxY max(res{a}(:))]);
    end

    for a = 1:3
        subplot(3,3,a+length(showRes)*(sr-1))

        for t = 1:3

            gType = g.expType{t};
            mCol = zeros(g.n(t),3); % make marker color
            
            if a==1
                mCol(gType == 2,1) = .5;
                mCol(gType == 3,1) = .9;
            elseif a == 2
                mCol(gType == 3,1) = .9;
            else
                mCol(gType == 2,1) = .9;
            end


            scat = scatter( g.xOffset{t} , res{t}(:,a),'filled');
            scat.CData = mCol;
            hold on
            for p = 1:length(g.isPaired{t})
                e = g.isPaired{t}(p);
                targ = find(g.expIds{t}==e);
                plot([ g.xOffset{t}(targ)   ],[res{t}(targ(1),a) res{t}(targ(2),a)],'k');
            end
        end

        ylim([0 maxY + maxY/5])

        xlim([0 7])
        titleStr = sprintf('%s %s',regionName{a},measureName{sr})
        title(titleStr)
    end

end



%% Show normalized results
clf

regionName = {'Full' 'Contra' 'Ipsi'};
measureName = {'Area' 'Count' 'Density'};

clear showRes
showRes{1} = g.areaFCI;
showRes{2} = g.countFCI;
for t = 1:3
    showRes{3}{t} = g.countFCI{t}./g.areaFCI{t};
end



for sr = 1:length(showRes) % each measurment type

    res = showRes{sr};
    maxY = 0;
    for a = 1:3

        maxY = max([maxY max(res{a}(:))]);
    end

    for a = 1:3 % each region
        subplot(3,3,a+length(showRes)*(sr-1))

        for t = 1:3 % each time point

            if t ==1
               normRef = median(res{t}(:,a));
            end

            gType = g.expType{t};
            mCol = zeros(g.n(t),3); % make marker color
            
            if a==1
                mCol(gType == 2,1) = .5;
                mCol(gType == 3,1) = .9;
            elseif a == 2
                mCol(gType == 3,1) = .9;
            else
                mCol(gType == 2,1) = .9;
            end


            scat = scatter( g.xOffset{t} , res{t}(:,a)/normRef,'filled');
            scat.CData = mCol;
            hold on
            for p = 1:length(g.isPaired{t})
                e = g.isPaired{t}(p);
                targ = find(g.expIds{t}==e);
                plot([ g.xOffset{t}(targ)   ],[res{t}(targ(1),a) res{t}(targ(2),a)]/normRef,'k');
            end
        end

        %ylim([0 3])

        xlim([0 7])
        titleStr = sprintf('%s %s\n%0.5f',regionName{a},measureName{sr},normRef)
        title(titleStr)
    end

end


%% get the intervals

standardDifCI(res{2}


for t = 1:3 % for each time
    for a = 1:3 % each region
        res2{t,a}


%         for p = 1:length(g.isPaired{t})
%             e = g.isPaired{t}(p);
%             targ = find(g.expIds{t}==e);
%             [res{t}(targ(1),a) res{t}(targ(2),a)];
%         end
    end
end


