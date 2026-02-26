function[] = startCam

global watProps

%%

try
    stop(watProps.cam)
end
watProps.cam = [];
useCam = 1;
% imaqhelp
% imaqfind
% imaqhwinfo
% info = imaqhwinfo('winvideo');
% dev_info = imaqhwinfo('winvideo',1);


%h = findobj('-regexp', 'Status','START')

imaqreset
try
    disp('starting camera')
    vid = videoinput('winvideo',1,'MJPG_1024x768');
    disp('camera started')
catch err
    disp('No camera detected')
    useCam = 0;
end

        
%flushdata, getdata, getsnapshot, peekdata
if useCam
    watProps.src = getselectedsource(vid);
    watProps.src.Exposure = -5;
    watProps.src.Gain = 1;
    watProps.src.BacklightCompensation = 'off';
    watProps.src.ExposureMode = 'manual';
    watProps.src.WhiteBalanceMode = 'auto';
    watProps.src.Gamma = 1;
    watProps.src.Brightness = 75;
    watProps.src.Contrast = 60;

    
    flushdata(vid)
    stop(vid)
    %wait(vid)
    'pausing for start'
    pause(1)
    src = getselectedsource(vid);
    triggerconfig(vid, 'Manual')
    %set(vid,'disklogger',
    vid.FramesPerTrigger = 1;
    vid.TriggerRepeat = Inf;
    %set(vid,'loggingmode','memory')
    start(vid)
    for testI = 1:10000
        pause(.1),testI
        if strcmp(get(vid,'running'),'on')
            break
        end
    end
    trigger(vid)
    for testI = 1:300
        pause(.1),testI
        if strcmp(get(vid,'logging'),'off')
            break
        end
        if testI == 200
            'Turn camera switch OFF and ON manually, wait for camera beep, and then press ON.'
           return
        end
    end
    I = getdata(vid);
    [ys xs cs] = size(I);
    
    
    
else
    
    ys = 512; xs = 512; cs = 3;
    I = uint8(rand(ys,xs,cs));
    vid = [];
end

image(I)
%%
watProps.useCam = useCam;
watProps.cam = vid;

watProps.iWidth = xs;
watProps.iHeight = ys;
watProps.iCol = cs;














