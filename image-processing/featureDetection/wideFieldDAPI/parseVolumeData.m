


%%Paste data from 'ReorderVol' sheet of 'Copy of LGN Data - Dapi Data' google spreadsheet
%{
expTypesTags = {}; % Rows 1:7
volumeDataRaw = {} % Rows 11:21
%}

%% Parse tags
expTypes = expTypesTags(:,1)
clear expNames
for i = 1:size(expTypes)

    expNames{i} = {};
    for c = 2:size(expTypesTags,2)
        nam = expTypesTags{i,c};
        if ~isempty(nam)
            expNames{i} = cat(2,expNames{i},{nam});
        end
    end   
end

%% Switch NAs
datVc = volumeDataRaw(2:end,2:end);
pick1 = [1 5 6 4 2 3 ];
pick2 = [];

for x = 1:size(datVc,2)
    c = datVc{2,x};
    if strcmp(class(c),'char');
        datVc(1:6,x) = datVc(pick1,x);
    end

    c = datVc{1,x};
    if strcmp(class(c),'char');
        datVc(1,x) = datVc(4,x);
        datVc{4,x} = c;
    end

end

datV = zeros(size(datVc));
for y = 1:size(datVc,1)
    for x = 1:size(datVc,2)
        c = datVc{y,x};
        if ~strcmp(class(c),'char');
            if ~isempty(c)
            datV(y,x) = c;
            end
        end
       
    end
end

%% Group experiments in data

dataNames = volumeDataRaw(1,2:end);
dataTypes = zeros(length(dataNames),1);
for i = 1:length(dataNames)
    nam = dataNames{i};

    for t = 1:length(expTypes)
        hit = sum(strcmp(expNames{t},nam));
        if hit >0
            if hit>1
                'multi hit'
            end
            dataTypes(i) = t;
            break
        end

    end
end

noType = dataNames(dataTypes==0)



%% Vol only
useT = [1 4 5];
v1 = 3;
v2 = 6;


clear g
g{1} = [datV(v1,dataTypes==useT(1)) datV(v2,dataTypes==useT(1))];
g{2} = datV(v1,dataTypes==useT(2));
g{3} = datV(v2,dataTypes==useT(2));
g{4} = datV(v1,dataTypes==useT(3));
g{5} = datV(v2,dataTypes==useT(3));

gCol = [0 0 0; .5 0 0; 0 0 .5; .5 0 0; 0 0 .5];
gX = [1 1.9 2.1 2.9 3.1];

gRef = g{1};
gRef = gRef(gRef>0);
gNorm = mean(gRef);

clf
subplot(1,4,1:2)
hold on

for i = 1:length(g)
    
    gs = g{i}/gNorm * 100;
    gs = gs(gs>0);
    m = mean(gs);
    se = std(gs)/sqrt(length(gs));
    
    scatter(gs*0+gX(i),gs,'markerEdgeColor','none','markerFaceColor',gCol(i,:),'markerfacealpha',.5)
    plot([gX(i)-.1 gX(i)+.1],[m+se m+se],'k','lineWidth',1.5)
    plot([gX(i)-.1 gX(i)+.1],[m-se m-se],'k','lineWidth',1.5)

end
xlim([0 4])
ylim([0 max(datV(:)/gNorm)*120])

param.useMean = 1;



gComp = [g{2} g{4}];
gComp = gComp(gComp>0);
%res = differenceRangeA1(gRef,gComp,param)
res = standardDifCI(gRef,gComp)


gComp = [g{5} g{3}];
gComp = gComp(gComp>0);
%res = differenceRangeA1(gRef,gComp,param)
res = standardDifCI(gRef,gComp)




%% Projections
useT = [1 4 5];
v1 = 2;

clear g
g{1} = [datV(v1,dataTypes==useT(1)) datV(v1+3,dataTypes==useT(1))];
g{2} = datV(v1,dataTypes==useT(2));
g{3} = datV(v1,dataTypes==useT(3));

gCol = [0 0 0 ; 0 0 0; 0 0 0];
gX = [1 2 3];

gRef = g{1};
gRef = gRef(gRef>0);
gNorm = mean(gRef);

subplot(1,4,3)
cla
hold on
for i = 1:length(g)
    
    gs = g{i}/gNorm * 100;
    gs = gs(gs>0);
    m = mean(gs);
    se = std(gs)/sqrt(length(gs));
    
    scatter(gs*0+gX(i),gs,'markerEdgeColor','none','markerFaceColor',gCol(i,:),'markerfacealpha',.5)
    plot([gX(i)-.2 gX(i)+.2],[m+se m+se],'k','lineWidth',1.5)
    plot([gX(i)-.2 gX(i)+.2],[m-se m-se],'k','lineWidth',1.5)

end
xlim([0 4])
ylim([0 max(gRef/gNorm)*120])



gComp = [g{2} g{3}];
gComp = gComp(gComp>0);
%res = differenceRangeA1(gRef,gComp,param)
res = standardDifCI(gRef,gComp)



%% Projections Ipsilateral
useT = [1 4 5];
v1 = 1;

clear g
g{1} = [datV(v1,dataTypes==useT(1)) datV(v1+3,dataTypes==useT(1))];
g{2} = datV(v1,dataTypes==useT(2));
g{3} = datV(v1,dataTypes==useT(3));

gCol = [0 0 0; 0 0 0; 0 0 0];
gX = [1 2 3];

gRef = g{1};
gRef = gRef(gRef>0);
gNorm = mean(gRef);

subplot(1,4,4)
cla
hold on
for i = 1:length(g)
    
    gs = g{i}/gNorm * 100;
    gs = gs(gs>0);
    m = mean(gs);
    se = std(gs)/sqrt(length(gs));
    
    scatter(gs*0+gX(i),gs,'markerEdgeColor','none','markerFaceColor',gCol(i,:),'markerfacealpha',.5)
    plot([gX(i)-.2 gX(i)+.2],[m+se m+se],'k','lineWidth',1.5)
    plot([gX(i)-.2 gX(i)+.2],[m-se m-se],'k','lineWidth',1.5)

end
xlim([0 4])
ylim([0 max(gRef/gNorm)*120])


gComp = [g{2} g{3}];
gComp = gComp(gComp>0);
%res = differenceRangeA1(gRef,gComp,param)
res = standardDifCI(gRef,gComp)




if 0
    epsName = 'D:\WorkDocs\Presentation\2023_RReSToRe\shadedScat.eps';
    print(gcf, epsName, '-depsc2','-painters','-r300')
end







