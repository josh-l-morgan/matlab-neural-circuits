function[] = runWatWatch(handles)

global watProps
wP = boundProps(handles);

timePosX = 10;
timePosY = 10;

set(handles.axes_Plot,'ylim',[0 255]);
winHeight = wP.y2 - wP.y1 + 1;

intMeanI = ones(winHeight,3,wP.intTime) * 256;

while 1
    wP = boundProps(handles);
    
    if wP.watching
        
        %% Show video and box
        
        if watProps.useCam 
            trigger(watProps.cam);
            I = getdata(watProps.cam);
        else
            I = uint8(rand(512)*256);
            iDir = dir(wP.sourceDir);
            dateNums = [iDir.datenum];
            bytes = [iDir.bytes];
            targ = find((dateNums == max(dateNums))& (bytes>0));
            if ~isempty(targ)
                
                targ = targ(end);
                targ = ceil(rand*(length(iDir)-10))+3; %!!!!!!!!!!!!!!!!
                iName = iDir(targ).name;
                I = imread([wP.sourceDir '\' iName]);
            else
                I = rand(512,512,3)*.2;
            end
            
        end
        
        if watProps.record
            filename = [watProps.recordDir '\i_' num2str(fix(datenum(datetime)*24*60*60*1000)) '.jpg'];
            imwrite(I,filename,'quality',25);
        end
        
        [ys xs cs] = size(I);
        wP.iHeight = ys;
        wP.iWidth = xs;
        wP.iCol = cs;
        
        %axes(handles.axes_Vid)
        image(handles.axes_Vid,I)
        hold(handles.axes_Vid,'on')
        text(handles.axes_Vid,timePosX, timePosY,datestr(datetime),'color','w')
        
        plot(handles.axes_Vid,[wP.x1 wP.x1], [wP.y1 wP.y2],'g')
        plot(handles.axes_Vid,[wP.x2 wP.x2], [wP.y1 wP.y2],'g')
        plot(handles.axes_Vid,[wP.x1 wP.x2], [wP.y1 wP.y1],'g')
        plot(handles.axes_Vid,[wP.x1 wP.x2], [wP.y2 wP.y2],'g')
              
        plot(handles.axes_Vid,[wP.x1 wP.x2], [wP.y1  wP.y1]+ winHeight * wP.win1,'b')
        plot(handles.axes_Vid,[wP.x1 wP.x2], [wP.y1  wP.y1]+ winHeight * wP.win2,'b')
        
        hold(handles.axes_Vid,'off')
        pause(.01)
        
        %% Cut out box and measure
        cutI = I(wP.y1:wP.y2,wP.x1:wP.x2,:);
        
        
        [cutSizeY cutSizeX cutSizeC] = size(cutI);
        
        if (cutSizeY ~= winHeight) | (wP.intTime ~= size(intMeanI,3))
            'doink'
            winHeight = wP.y2 - wP.y1 + 1;
            intMeanI = ones(winHeight,3,wP.intTime) * 256;
        end
        
        allMeanI = squeeze(mean(cutI,2));
        meanI = mean(allMeanI,2);
        
        intMeanI(:,:,1:end-1) = intMeanI(:,:,2:end);
        intMeanI(:,:,end) = allMeanI;
        meanIntMeanI = mean(intMeanI,3);
        
        
        %%Color mode
        if watProps.colorMode == 1
            colMean = mean(meanIntMeanI,2);
        elseif watProps.colorMode == 2 
            colMean = mean(meanIntMeanI(:,1:2),2);
        else     
           colMean = mean(meanIntMeanI,2);
           ePos = getEdge(colMean); 
        end
        
        %axes(handles.axes_Plot);

        hold(handles.axes_Plot,'off')
        plot(handles.axes_Plot,meanIntMeanI(:,1),'color','r','lineWidth',3)
        hold(handles.axes_Plot,'on')
        plot(handles.axes_Plot,meanIntMeanI(:,2),'color','g','lineWidth',3)
        plot(handles.axes_Plot,meanIntMeanI(:,3),'color','b','lineWidth',3)
        plot(handles.axes_Plot,colMean,'color','w','linewidth',3)
        plot(handles.axes_Plot,meanI','color',[.6 .6 .6])

        set(handles.axes_Plot,'color','k')
        view(handles.axes_Plot,[90 -90])
        
        valBot = -10;
        valTop = max(meanIntMeanI(:)) * 1.1;
        plot(handles.axes_Plot,[wP.win1 wP.win1] * cutSizeY,[valBot valTop],'b')
        plot(handles.axes_Plot,[wP.win2  wP.win2]* cutSizeY,[valBot valTop],'b')
        
        
        if watProps.colorMode == 3
            plot(handles.axes_Plot,[ePos.firstCross ePos.firstCross],[valBot valTop],'y','lineWidth',3)
        else
            
            sampMean = colMean(round(winHeight * wP.win1):round(winHeight * wP.win2));
            currentMean = mean(sampMean);
            plot(handles.axes_Plot,[ wP.win1   wP.win2]* winHeight,[currentMean currentMean],'y','lineWidth',3)
            
        end
        plot(handles.axes_Plot,[0  length(meanI)],[wP.thresh1 wP.thresh1],'r')
        
        set(handles.axes_Plot,'ylim',[valBot valTop]);
        set(handles.axes_Plot,'xlim',[0 length(meanI)]);
        set(handles.axes_Plot,'Xdir','reverse')
        pause(.01)
        
        
        %%Trigger pump
        if watProps.colorMode == 3
            shouldTrigger = ePos.tooLow;
        else
            shouldTrigger = currentMean < wP.thresh1;
        end
        
        if shouldTrigger
            secLast = toc(watProps.lastPump);
            if secLast > watProps.pumpInterval
                
                if get(handles.togglebutton_threshPump,'Value')
                    watProps.threshHist = cat(1,watProps.threshHist,datenum(datetime));
                    set(handles.text_Status,'backgroundcolor',[1 1 0 ]);
                    set(handles.text_Status,'string','Auto PUMPING')
                    triggerPump
                
                else
                    'Pump would trigger if toggled on'
                    if watProps.soundOn
                        sound(sin(1:1000)/10)
                    end
                end
                watProps.lastPump = tic;

            end
        end
        
        %% History
        
        histUnit = get(handles.listbox_histUnit,'Value');
        scaleHist = [24*60*60 24*60 24 24];
        maxHist = [60 60 24 7*24];
        nameHist = {'seconds' 'minutes' 'hours' 'hours'};
        
        %axes(handles.axes_History)
        text(handles.axes_History,1,1,'')
        threshHist = (datenum(datetime) - watProps.threshHist) * scaleHist(histUnit);
        
        scatter(threshHist,ones(length(threshHist),1),'k')
        hold(handles.axes_History,'on')
        manHist = (datenum(datetime)-watProps.manHist) * scaleHist(histUnit);
        scatter(handles.axes_History,manHist,ones(length(manHist),1),'r')
        set(handles.axes_History,'xlim',[0 maxHist(histUnit)])
        hold(handles.axes_History,'off')
        
        
        %% Status
        
        if get(handles.togglebutton_threshPump,'Value')
            watProps.status = 3;
        else
            watProps.status = 2;
        end
       if watProps.watching == 0
            set(handles.text_Status,'backgroundcolor',[.6 .6 .6 ]);
            set(handles.text_Status,'string','OFF and waiting')
       elseif watProps.status == 2
            set(handles.text_Status,'backgroundcolor',[1 .5 .5 ]);
            set(handles.text_Status,'string','Watching, but NOT pumping')
        elseif watProps.status == 3
            set(handles.text_Status,'backgroundcolor',[.5 1 .5]);
            set(handles.text_Status,'string','Auto pump is active')
        else
            set(handles.text_Status,'backgroundcolor',[.7 .7 .4]);
            set(handles.text_Status,'string','confused')
        end
        
        
        pause(.01)
    else
        'done watching'
        break
    end
    
end




