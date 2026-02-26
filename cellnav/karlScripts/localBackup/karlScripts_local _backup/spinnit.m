function out=spinnit(stepAng,figureHandle,fileNamePrefix,saveDir)
    out=1;
    fileDir=saveDir;
    figure(figureHandle);
    hold on
    view(-180,0);
    drawnow
    set(gca,'visible','off')
    for i=1:360/stepAng
        ang=(-180+(i-1)*stepAng);
        fileName=[fileDir fileNamePrefix num2str(i,'%03.f') '.jpg'];
        view(ang,0);
        drawnow
        set(gca,'visible','off')
        fim = getframe(figureHandle);
        testIm = frame2im(fim);
        imwrite(testIm,fileName);
    end
end