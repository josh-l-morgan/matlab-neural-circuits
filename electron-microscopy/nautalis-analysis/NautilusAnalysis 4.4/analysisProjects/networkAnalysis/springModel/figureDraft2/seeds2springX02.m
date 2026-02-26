

itag  = 'giant';

springDir = 'D:\LGNs1\Analysis\springDat\figDraft1b\';
springRes = 'D:\LGNs1\Analysis\springDat\results\';
if ~exist(springDir,'dir'), mkdir(springDir), end
if ~exist(springRes,'dir'), mkdir(springRes), end

load([springRes 'res_all_Phage_edit3.mat'])
load([springRes 'fourSeeds_edit1.mat'])
load([springRes 'fourSeeds_edit1.mat'])



%% Get data

springIn = springUse_cladeGroups(); % get data
springDat = springParameters_X01_figDraft(springIn); % set spring parameters
springDat = springParameters_1D1d_X1_figDraft(springIn);

%  springDat.edges.ew = springDat.edges.ew>0;

%% Run springs

for rerun = 1: 1
    
    % allResults{rerun} = runSprings(springDat)
    allResults{rerun} = runSprings1D1D(springDat)
end

%% Save

if 0
%% Save image

    tag = 'cladeNet';
    
    
       set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])

     %runSprings(springDat,allResults{1})
    set(gcf, 'InvertHardCopy', 'off');
    imageName = sprintf('%sspringRun_%s.png',springDir,tag);
    %print(gcf,imageName,'-dpng','-r1024','-opengl','-noui')
    
    epsName = sprintf('%sspringRun_%s.eps',springDir,tag);
    print(gcf, epsName, '-depsc2','-painters','-r300')
    
    
    
    
    
    

%%
end




%{ 

%% Save Results

result = allResults{rerun};
save([springRes 'res_all_05.mat'],'result')

%}


if 0
%% Scale bar

    clf
    set(gcf,'Position',[1100 100 800 800],'visible','on')
    set(gca,'Position',[.5 .05 .05 .9])
    set(gcf,'color','k')
    colormap(colTable)
    image([0:size(colTable,1)]')
    ylim([1 100])
    set(gca,'Ydir','reverse')
    ticPos(1) = 1;
    set(gca,'YTick',ticPos,'YTickLabel',showRange,'XTick',[],'Ycolor','w')

    epsName = sprintf('%s%s_scaleBar.eps',springDir,itag);
    print(gcf, epsName, '-depsc2','-painters','-r300')


end
