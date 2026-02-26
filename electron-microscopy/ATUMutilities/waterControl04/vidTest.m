imaqhwinfo
info = imaqhwinfo('winvideo');
dev_info = imaqhwinfo('winvideo',1);
%vid = videoinput('winvideo',2);
vid = videoinput('winvideo',1,'MJPG_1024x768');

%preview(vid)
start(vid)
tic
for i=1:10
   disp(i)
   curImg=getsnapshot(vid);
   image(curImag)
   pause(.01);
end
toc
preview(vid)

%%
imaqreset
frames = 5
vid = videoinput('winvideo',2);
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
start(vid)
tic;
for j=1:frames
    data = getdata(vid,1);
    %data = getsnapshot(vid);
    image(data)
    pause(.01)
end
tstop = toc;
framerate = frames / tstop
stop(vid)

%%
imaqreset
frames = 50
%vid = videoinput('winvideo',2);
vid = videoinput('winvideo',1,'MJPG_1024x768');

src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
%vid.ROIPosition = [1 1 10 10];
triggerconfig(vid,'manual');
start(vid)
tic;
for j=1:frames
    trigger(vid);
    data = getdata(vid);
    image(data)
    pause(.01)
end
tstop = toc;
framerate = frames / tstop
stop(vid)


flushdata(vid)




