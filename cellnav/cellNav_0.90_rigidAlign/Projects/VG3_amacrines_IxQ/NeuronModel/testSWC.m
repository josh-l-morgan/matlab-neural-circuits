
%% Load data
if 0
    sm2nmFile = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\sm2nrn_cid11.mat'
    load(sm2nmFile)
    swcS = sm2nrn.nep.swcS;
    pred = swcS.pred +1;
    pos = swcS.pos;

elseif 1
    smFile = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\SMs\sm_cid11.mat';
    load(smFile);
    swcS = sm.nep.swcS;
    pred = swcS.pred +1;
    pos = swcS.pos;

else
    swcFile = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\cid11.swc';
    load(swcFile);
end

%% Make figure
try,close(f), end
f = figure;
axis equal
f.Clipping = 0;

%% render
clf
scatter3(pos(:,1),pos(:,2),pos(:,3),'.','k')
hold on




root = find(pred==0);
tips = setdiff(1:length(pred),unique(pred));

for t = 1:length(tips)
    c = tips(t);
    for p = 1:length(pred);
        scatter3(pos(c,1),pos(c,2),pos(c,3),'r')
        drawnow
        c = pred(c);
        if ~c
            break
        end
    end
end

hold off


%% cleck ordering
for p = 1:length(pred)
    orderOK(p) = pred(p)<p;
end
