function[outtrack,currentAx] = getMyIn();

%% from mousetrack2
set(0, 'units', 'pixels'); 
set(gcf, 'units','pixels'); 
set(gca, 'units','pixels');

set(gca, 'ydir', 'reverse');   % 'normal'  Y axis direction 
  
set(gcf, 'KeyPressFcn',{@keyfig});
set(gcf, 'windowbuttondownfcn', {@starttrack}); 
set(gcf, 'windowbuttonupfcn', {@stoptrack}, 'userdata', 10); 
set(gcf, 'Interruptible', 'on'); set(gcf, 'BusyAction', 'queue')    %  queue not loose  actions ...
set(gca, 'Interruptible', 'on'); set(gca, 'BusyAction', 'cancel')
set(gca, 'DrawMode', 'fast') 
set(gca, 'XLimMode', 'manual','YLimMode','manual', 'ZLimMode', 'manual')
set(gcf, 'Renderer', 'painters'); set(gcf, 'DoubleBuffer','on')   % speeds up render, prevents blinking
%set(gcf,'Menubar','none')

pause(.03)
set(gcf,'userdata',0)
waitfor(gcf,'userdata',1);
set(gcf,'userdata',0)
outtrack = get(gca,'userdata');
currentAx = gca;
set(gca,'userdata',[])

function starttrack(imagefig, varargins) 
CurPnt = get(gca, 'CurrentPoint');
coords = CurPnt(2,1:2); 
set(gca, 'userdata', coords );   % disp('tracking started') 
%set(gcf, 'Pointer', 'crosshair') 
%set(gcf, 'windowbuttondownfcn', {@stoptrack}, 'userdata', []);
set(gcf, 'windowbuttonmotionfcn', {@followtrack}); 


function followtrack(imagefig, varargins) 
CurPnt = get(gca, 'CurrentPoint'); 
coords = CurPnt(2,1:2); 
set(gca, 'userdata', [get(gca,'userdata'); coords]); 
% hold on
% plot(coords(1), coords(2), 'r.', 'MarkerSize', 10); 
% hold off

%--------------------------------------------------------------------------
function stoptrack(imagefig, varargins)
%set(gcf, 'Pointer', 'arrow')
CurPnt = get(gca, 'CurrentPoint'); 
coords = CurPnt(2,1:2); 
set(gca, 'userdata', [get(gca,'userdata'); coords]);
set(gcf, 'windowbuttonmotionfcn', []);
set(gcf,'userdata',1)


%------------------------------------------------------------------
function keyfig(src,evnt)
    set(gca,'userdata',evnt)
    set(gcf,'userdata',1)


