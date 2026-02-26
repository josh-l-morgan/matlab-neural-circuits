function[] = saveSpring(pathName)


    slash = find(pathName == '\');
    itag  = pathName(slash(end)+1:end);
    if isempty(itag)
        itag = 'spring';
    end
    
    springDir = pathName(1:slash(end));
    if ~exist(springDir,'dir'),mkdir(springDir),end
    
    inDir = dir([springDir '*.png']);
    newNum = length(inDir)+1;
    
    set(gcf, 'InvertHardCopy', 'off');
    
     set(gcf,'PaperUnits','points','PaperPosition',[1 1 700 700])
    
    
    imageName = sprintf('%s%s_%07.0f.png',springDir,itag,newNum);
    print(gcf,imageName,'-dpng','-r150','-opengl','-noui')
%     
%     epsName = sprintf('%s%s_%07.0f.eps',springDir,itag,newNum);
%     print(gcf, epsName, '-depsc2','-painters','-r300')