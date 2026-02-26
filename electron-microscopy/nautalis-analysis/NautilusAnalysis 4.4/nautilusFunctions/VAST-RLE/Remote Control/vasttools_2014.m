% VastTools.m
% Additional tools for VAST in Matlab
% Version 1.01 by Daniel Berger, August 2014 - October 02 2015

function [] = vasttools()
  global vdata;
  if isfield(vdata,'state')
    warndlg('VastTools is already running. Close other instance. If necessary use "clear -global vdata" to fix.','Error starting VastTools');
    return;
  end;
  
  p = mfilename('fullpath');
  p = p(1:end-size(mfilename(),2));
  addpath(p);
  addpath([p 'VASTControl/']);
  vast=VASTControlClass();  %CAUTION: This contains javaaddpath which clears all globals. WTF, Matlab !

  global vdata;
  vdata.vast=vast;
  vdata.state.isconnected=0;
  vdata.state.connect.ip='127.0.0.1';
  vdata.state.connect.port=22081;
  vdata.state.lastcancel=1;
  vdata.state.guiblocked=0; 
  vdata.data.nroftargetlists=0;
  vdata.data.nrofsimplenavigators=0;
  vdata.etfh=[];
  
  scrsz = get(0,'ScreenSize');
  vdata.fh = figure('units','pixels',...
    'outerposition',[10 scrsz(4)-639-10 640 200],...
    'menubar','none',...
    'numbertitle','off',...
    'name','VastTools Version 1.01 (c) 2015 by Daniel Berger');

  set(vdata.fh,'CloseRequestFcn',{@callback_quit});
  set(vdata.fh,'ResizeFcn',{@callback_resize},'resize','on');
  
    %%%%%% MAIN MENU %%%%%%%
  
  vdata.ui.menu.connectmenu = uimenu(vdata.fh,'Label','Connect');
  vdata.ui.menu.connectvast = uimenu(vdata.ui.menu.connectmenu,'Label','Connect to VAST','Callback',{@callback_connect});
  
  vdata.ui.connectoptions = uimenu(vdata.ui.menu.connectmenu,'Label','Connection Options','Callback',{@callback_connectionoptions});
  vdata.ui.quit = uimenu(vdata.ui.menu.connectmenu,'Label','Quit','Callback',{@callback_quit},'Separator','on','Accelerator','Q');
  vdata.ui.menu.exportmenu = uimenu(vdata.fh,'Label','Export');
  vdata.ui.menu.exportobj = uimenu(vdata.ui.menu.exportmenu,'Label','Export 3D Surfaces as OBJ Files ...','Callback',{@callback_exportobj},'Enable','off');
  vdata.ui.menu.exportprojection = uimenu(vdata.ui.menu.exportmenu,'Label','Export Projection Image ...','Callback',{@callback_exportprojection},'Enable','off');
  vdata.ui.menu.measuremenu = uimenu(vdata.fh,'Label','Measure');
  vdata.ui.menu.measurevolumes = uimenu(vdata.ui.menu.measuremenu,'Label','Measure Segment Volumes ...','Callback',{@callback_measurevol},'Enable','off');
  vdata.ui.menu.measurelengths = uimenu(vdata.ui.menu.measuremenu,'Label','Measure Segment Lengths ...','Callback',{@callback_measurelength},'Enable','off');
  vdata.ui.menu.euclidiantool = uimenu(vdata.ui.menu.measuremenu,'Label','Euclidian Distance Measurement Tool','Callback',{@callback_euclidiantool},'Enable','off');
  vdata.ui.menu.targetlistmenu = uimenu(vdata.fh,'Label','Target List');
  vdata.ui.menu.newtargetlist = uimenu(vdata.ui.menu.targetlistmenu,'Label','New Target List ...','Callback',{@callback_newtargetlist});
  vdata.ui.menu.loadtargetlist = uimenu(vdata.ui.menu.targetlistmenu,'Label','Load Target List ...','Callback',{@callback_loadtargetlist});
  vdata.ui.menu.navigatemenu = uimenu(vdata.fh,'Label','Navigate');
  vdata.ui.menu.newsimplenavigator = uimenu(vdata.ui.menu.navigatemenu,'Label','New Simple Navigator Image From Last Projection Image ...','Callback',{@callback_newsimplenavigator});
  vdata.ui.menu.loadsimplenavigator = uimenu(vdata.ui.menu.navigatemenu,'Label','Load Simple Navigator Image From File ...','Callback',{@callback_loadsimplenavigator});
  
  pos=get(vdata.fh,'Position');
  vdata.ui.cancelbutton = uicontrol('style','push','units','pixels','position',[12 pos(4)-62 55 25],...
    'string','Cancel','Enable','off','callback',{@callback_canceled});
  vdata.ui.message = uicontrol('style','text','unit','pix','position',[85 pos(4)-110 pos(3)-10 100],'fontsize',11,'string','Idle','backgroundcolor',[0.75 0.75 0.65]);
  set(vdata.ui.message,'String',{'DISCLAIMER', 'VastTools is provided as-is. You are using the functions herein, especially the functions for measurement and analysis, at your own risk. This software may contain bugs and produce wrong results. Please make sure your data set is scaled correctly in VAST (check Info / Volume Properties). In case you encounter any bugs, please let me know.'});
  
  
function [] = callback_quit(varargin)
  global vdata;
  
  %Check if there are unsaved target lists
  if (vdata.data.nroftargetlists>0)
    changedtlexists=0;
    
    for instance=1:1:vdata.data.nroftargetlists
      if ishandle(vdata.data.tl(instance).fh)
        if (vdata.data.tl(instance).ischanged==1)
          changedtlexists=changedtlexists+1;
        end;
      end;
    end;

    if (changedtlexists==1)
      res = questdlg('You have a changed target list open. Are you sure you want to quit without saving?','Quit VastTools','Yes','No','Yes');
      if strcmp(res,'No')
        return;
      end
    end;
    if (changedtlexists>1)
      res = questdlg('You have changed target lists open. Are you sure you want to quit without saving?','Quit VastTools','Yes','No','Yes');
      if strcmp(res,'No')
        return;
      end
    end;
  end;
    
  try
    %%%% CLEANUP
    % Disconnect if XTLibServer is connected
    if (vdata.state.isconnected==1)
      vdata.vast.disconnect();
    end;

    % Close simple navigator windows
    if (vdata.data.nrofsimplenavigators>0)
      for instance=1:1:vdata.data.nrofsimplenavigators
        if ishandle(vdata.data.sn(instance).fh)
          delete(vdata.data.sn(instance).fh);
        end
        vdata.data.sn(instance).open=0;
      end;
    end;
    
    % Close target list windows
    if (vdata.data.nroftargetlists>0)
      for instance=1:1:vdata.data.nroftargetlists
        if ishandle(vdata.data.tl(instance).fh)
          delete(vdata.data.tl(instance).fh);
        end
        vdata.data.tl(instance).open=0;
      end;
    end;
    
    % Close euclidian tool window if open
    if ishandle(vdata.etfh)
      delete(vdata.etfh);
    end;
    
    % Close main window
    if ishandle(vdata.fh)
      delete(vdata.fh);
    end
  catch err
    %If something went wrong, delete the current figure.
    delete(gcf);
  end;
  clear -global vdata;
  
  
function [] = callback_resize(varargin)
  global vdata;
  set(vdata.fh,'Units','pixels');
  pos = get(vdata.fh,'OuterPosition');
  hpos=pos(3)+(-1024+560);
  vpos=pos(4)-100;
  pos=get(vdata.fh,'Position');
  set(vdata.ui.cancelbutton,'position',[5 pos(4)-35 55 25]);
  set(vdata.ui.message,'position',[65 5 pos(3)-70 pos(4)-10]);


function [] = updategui()
  global vdata;
  if (vdata.state.guiblocked)
    set(vdata.ui.menu.exportobj,'Enable','off');
    set(vdata.ui.menu.exportprojection,'Enable','off');
    set(vdata.ui.menu.measurevolumes,'Enable','off');
    set(vdata.ui.menu.euclidiantool,'Enable','off');
  else
    set(vdata.ui.menu.exportobj,'Enable','on');
    set(vdata.ui.menu.exportprojection,'Enable','on');
    set(vdata.ui.menu.measurevolumes,'Enable','on');
    set(vdata.ui.menu.euclidiantool,'Enable','on');
  end;
  
function [] = blockgui()
  global vdata;
  vdata.state.guiblocked=1; 
  updategui();
  
function [] = releasegui()
  global vdata;
  vdata.state.guiblocked=0; 
  updategui();
  

function [] = callback_done(varargin)
  global vdata;
  vdata.state.lastcancel=0;
  vdata.ui.temp.closefig=1;
  uiresume(gcbf);
  
  
function [] = callback_canceled(varargin)
  global vdata;
  vdata.state.lastcancel=1;
  vdata.ui.temp.closefig=1;
  uiresume(gcbf);


function [] = callback_connect(varargin)
  global vdata;
  
  if (vdata.state.isconnected==0)
    %Try to connect
    res=vdata.vast.connect(vdata.state.connect.ip,vdata.state.connect.port,1000);
    if (res==0)
      warndlg(['ERROR: Connecting to VAST at ' vdata.state.connect.ip ' port ' sprintf('%d',vdata.state.connect.port) ' failed. Please enable the Remote Control Server in VAST (in the main menu under "Window", "Remote Control API Server"; click "Enable") and make sure the IP and Port settings are correct!'],'Error connecting to VAST');
      return;
    end;
    vdata.state.isconnected=1;
    set(vdata.ui.menu.connectvast,'Label','Disconnect from VAST');
    set(vdata.ui.menu.connectmenu,'Label','Disconnect');
    set(vdata.ui.menu.exportobj,'Enable','on');
    set(vdata.ui.menu.exportprojection,'Enable','on');
    set(vdata.ui.menu.measurevolumes,'Enable','on');
    %set(vdata.ui.menu.measurelengths,'Enable','on');
    set(vdata.ui.menu.euclidiantool,'Enable','on');
  else
    %Disconnect
    res=vdata.vast.disconnect();
    if (res==0)
      warndlg('ERROR: Disconnecting from VAST failed.','Error disconnecting from VAST');
      return;
    end
    vdata.state.isconnected=0;
    set(vdata.ui.menu.connectvast,'Label','Connect to VAST');
    set(vdata.ui.menu.connectmenu,'Label','Connect');
    set(vdata.ui.menu.exportobj,'Enable','off');
    set(vdata.ui.menu.exportprojection,'Enable','off');
    set(vdata.ui.menu.measurevolumes,'Enable','off');
    %set(vdata.ui.menu.measurelengths,'Enable','off');
    set(vdata.ui.menu.euclidiantool,'Enable','off');
  end;
  
  
function [] = callback_connectionoptions(varargin)
  global vdata;
  
  %blockgui();
  scrsz = get(0,'ScreenSize');
  figheight=160;
  f = figure('units','pixels','position',[50 scrsz(4)-100-figheight 360 figheight],'menubar','none','numbertitle','off','name','VastTools - Connection Options','resize','off');
  vpos=figheight-40;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 160 15],'String','IP address of VAST computer:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos-20 300 15],'String','(use 127.0.0.1 if VAST runs on the same computer)','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e1 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[200 vpos 100 20],'String',vdata.state.connect.ip,'horizontalalignment','left');
  vpos=vpos-60;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 160 15],'String','Port Address (default is 22081):','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e2 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[200 vpos 100 20],'String',sprintf('%d',vdata.state.connect.port),'horizontalalignment','left');
  vpos=vpos-30;
  
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[80 20 60 20], 'String','OK', 'CallBack',{@callback_done});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[220 20 60 20], 'String','Cancel', 'CallBack',{@callback_canceled});

  vdata.state.lastcancel=1;
  vdata.ui.temp.closefig=0;
  uiwait(f);
  
  if (vdata.state.lastcancel==0)
    vdata.state.connect.ip=get(e1,'String');
    vdata.state.connect.port = str2num(get(e2,'String'));
  end;
  
  if (vdata.ui.temp.closefig==1) %to distinguish close on button press and close on window x
    close(f);
  end;
  %releasegui();
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Surface OBJ Exporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = callback_exportobj(varargin)
  global vdata;

  if (vdata.state.isconnected==0)
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
    return;
  end;
  
  vinfo=vdata.vast.getinfo();
  if (min(size(vinfo)==0))
    warndlg('ERROR: Requesting data from VAST failed.','Error with remote connection to VAST');
    return;
  end;
  
  if (min([vinfo.datasizex vinfo.datasizey vinfo.datasizez])==0)
    warndlg('ERROR: No volume open in VAST.','VastTools OBJ exporting');
    return;
  end;
  
  nrofsegments=vdata.vast.getnumberofsegments();
  if (nrofsegments==0)
    warndlg('ERROR: No segmentation available in VAST.','VastTools OBJ exporting');
    return;
  end;
  
  blockgui();
  
  %Display parameter dialog
  if (~isfield(vdata.data,'region'))
    vdata.data.region.xmin=0;
    vdata.data.region.xmax=vinfo.datasizex-1;
    vdata.data.region.ymin=0;
    vdata.data.region.ymax=vinfo.datasizey-1;
    vdata.data.region.zmin=0; %first slice
    vdata.data.region.zmax=vinfo.datasizez-1; %last slice
  else
    if (vdata.data.region.xmin<0) vdata.data.region.xmin=0; end;
    if (vdata.data.region.xmax>(vinfo.datasizex-1)) vdata.data.region.xmax=vinfo.datasizex-1; end;
    if (vdata.data.region.ymin<0) vdata.data.region.ymin=0; end;
    if (vdata.data.region.ymax>(vinfo.datasizey-1)) vdata.data.region.ymax=vinfo.datasizey-1; end;
    if (vdata.data.region.zmin<0) vdata.data.region.zmin=0; end; %first slice
    if (vdata.data.region.zmax>(vinfo.datasizez-1)) vdata.data.region.zmax=vinfo.datasizez-1; end;
  end;
  if (~isfield(vdata.data,'exportobj'))
    vdata.data.exportobj.miplevel=0;
    vdata.data.exportobj.slicestep=1;      %4 means for example that every 4th slice exists (0, 4, 8, 12, ...)

    vdata.data.exportobj.blocksizex=1024; %Data block size for processing. For small data sets, make this a bit larger than the data (otherwise objects may be open)
    vdata.data.exportobj.blocksizey=1024;
    vdata.data.exportobj.blocksizez=64;
    vdata.data.exportobj.overlap=1;     %Leave this at 1
    vdata.data.exportobj.xscale=0.001;  %Use these to scale the exported models
    vdata.data.exportobj.yscale=0.001;
    vdata.data.exportobj.zscale=0.001;
    vdata.data.exportobj.xunit=vinfo.voxelsizex;%6*4;  %in nm
    vdata.data.exportobj.yunit=vinfo.voxelsizey; %6*4;  %in nm
    vdata.data.exportobj.zunit=vinfo.voxelsizez; %30; %in nm
    vdata.data.exportobj.outputoffsetx=0; %to translate the exported models in space
    vdata.data.exportobj.outputoffsety=0;
    vdata.data.exportobj.outputoffsetz=0;
    vdata.data.exportobj.invertz=1;

    vdata.data.exportobj.extractwhich=2;
    vdata.data.exportobj.objectcolors=1;
    vdata.data.exportobj.targetfileprefix='Segment_';
    vdata.data.exportobj.targetfolder=pwd;
    vdata.data.exportobj.includefoldernames=1;
    vdata.data.exportobj.closesurfaces=1;
    vdata.data.exportobj.skipmodelgeneration=0;
    vdata.data.exportobj.write3dsmaxloader=1;
    vdata.data.exportobj.savesurfacestats=0;
    vdata.data.exportobj.surfacestatsfile='surfacestats.txt';
  else
    if (vdata.data.exportobj.miplevel>(vinfo.nrofmiplevels-1)) vdata.data.exportobj.miplevel=vinfo.nrofmiplevels-1; end;
  end;
  
  scrsz = get(0,'ScreenSize');
  figheight=630;
  f = figure('units','pixels','position',[50 scrsz(4)-100-figheight 500 figheight],'menubar','none','numbertitle','off','name','VastTools - Export 3D Objects as OBJ Files','resize','off');

  vpos=figheight-40;
 
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15], 'Tag','t1','String','Render at resolution:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(vinfo.nrofmiplevels,1);
  vx=vinfo.voxelsizex;
  vy=vinfo.voxelsizey;
  vz=vinfo.voxelsizez;
  for i=1:1:vinfo.nrofmiplevels
    str{i}=sprintf('Mip %d - (%.2f nm, %.2f nm, %.2f nm) voxels',i-1,vx,vy,vz);
    vx=vx*2; vy=vy*2;
  end;
  pmh = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportobj.miplevel+1,'Position',[170 vpos 290 20]);
  vpos=vpos-30;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Use every nth slice:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e1 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.exportobj.slicestep),'horizontalalignment','left');
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15],'String','Render from area:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos 120 20], 'String','Set to full', 'CallBack',{@callback_region_settofull,0});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-25 120 20], 'String','Set to current voxel', 'CallBack',{@callback_region_settocurrentvoxel,0});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-50 120 20], 'String','Extend to current voxel', 'CallBack',{@callback_region_extendtocurrentvoxel,0});
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','X min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d', vdata.data.region.xmin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','X max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.xmax),'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Y min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Y max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymax),'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Z min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Z max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmax),'horizontalalignment','left');
  vpos=vpos-40;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Voxel size (full res)  X:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e8 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%f', vdata.data.exportobj.xunit),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[240 vpos 150 15],'String','Y:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e9 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[260 vpos 50 20],'String',sprintf('%f', vdata.data.exportobj.yunit),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[330 vpos 150 15],'String','Z:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e10 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[350 vpos 50 20],'String',sprintf('%f', vdata.data.exportobj.zunit),'horizontalalignment','left');
  vpos=vpos-20;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[60 vpos 400 15], 'Tag','t1','String',sprintf('[VAST reports the voxel size to be: (%.2f nm, %.2f nm, %.2f nm)]',vinfo.voxelsizex,vinfo.voxelsizey,vinfo.voxelsizez),'backgroundcolor',get(f,'color'),'horizontalalignment','left');
  set(t,'tooltipstring','To change, enter the values in VAST under "Info / Volume properties" and save to your EM stack file.');
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Scale models by   X:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e11 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.xscale),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[240 vpos 150 15],'String','Y:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e12 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[260 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.yscale),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[330 vpos 150 15],'String','Z:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e13 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[350 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.zscale),'horizontalalignment','left');
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Model output offset   X:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e14 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.outputoffsetx),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[240 vpos 150 15],'String','Y:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e15 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[260 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.outputoffsety),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[330 vpos 150 15],'String','Z:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e16 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[350 vpos 50 20],'String',sprintf('%f',vdata.data.exportobj.outputoffsetz),'horizontalalignment','left');
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Processing block size   X:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e17 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.exportobj.blocksizex),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[240 vpos 150 15],'String','Y:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e18 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[260 vpos 50 20],'String',sprintf('%d',vdata.data.exportobj.blocksizey),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[330 vpos 150 15],'String','Z:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e19 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[350 vpos 50 20],'String',sprintf('%d',vdata.data.exportobj.blocksizez),'horizontalalignment','left');
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15], 'Tag','t1','String','Export what:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(4,1);
  str{1}='All segments individually, uncollapsed';
  str{2}='All segments, collapsed as in VAST';
  str{3}='Selected segment and children, uncollapsed';
  str{4}='Selected segment and children, collapsed as in VAST';
  pmh2 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportobj.extractwhich,'Position',[120 vpos 290 20]);
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','File name prefix:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e20 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[120 vpos 290 20],'String',vdata.data.exportobj.targetfileprefix,'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Object Colors:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(2,1);
  str{1}='Object colors from VAST';
  str{2}='Object volumes as JET colormap';
  pmh3 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportobj.objectcolors,'Position',[120 vpos 290 20]);
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Target folder:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e21 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[120 vpos 290 20],'String',vdata.data.exportobj.targetfolder,'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[420 vpos 60 20], 'String','Browse...', 'CallBack',{@callback_exportobj_browse});
  vpos=vpos-30;
  
  c1 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 250 15],'Value',vdata.data.exportobj.includefoldernames,'string','Include Vast folder names in file names','backgroundcolor',get(f,'color')); 
  c2 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[300 vpos 200 15],'Value',vdata.data.exportobj.invertz,'string','Invert Z axis','backgroundcolor',get(f,'color')); 
  vpos=vpos-25;
  c3 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[300 vpos 200 15],'Value',vdata.data.exportobj.closesurfaces,'string','Close surface sides','backgroundcolor',get(f,'color')); 
  c4 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 250 15],'Value',vdata.data.exportobj.write3dsmaxloader,'string','Write 3dsMax bulk loader script to folder','backgroundcolor',get(f,'color')); 
  vpos=vpos-25;
  c5 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 250 15],'Value',vdata.data.exportobj.skipmodelgeneration,'string','Skip model file generation','backgroundcolor',get(f,'color')); 
  vpos=vpos-30;
  
  c6 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 250 15],'Value',vdata.data.exportobj.savesurfacestats,'string','Save surface statistics to file:','backgroundcolor',get(f,'color')); 
  e21 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[210 vpos 200 20],'String',vdata.data.exportobj.surfacestatsfile,'horizontalalignment','left');
  
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[150 20 60 20], 'String','OK', 'CallBack',{@callback_done});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[290 20 60 20], 'String','Cancel', 'CallBack',{@callback_canceled});

  vdata.state.lastcancel=1;
  vdata.ui.temp.closefig=0;
  uiwait(f);
  
  if (vdata.state.lastcancel==0)
    vdata.data.exportobj.miplevel=get(pmh,'value')-1;
    vdata.data.exportobj.slicestep = str2num(get(e1,'String'));
    vdata.data.region.xmin = str2num(get(vdata.temp.e_xmin,'String'));
    vdata.data.region.xmax = str2num(get(vdata.temp.e_xmax,'String'));
    vdata.data.region.ymin = str2num(get(vdata.temp.e_ymin,'String'));
    vdata.data.region.ymax = str2num(get(vdata.temp.e_ymax,'String'));
    vdata.data.region.zmin = str2num(get(vdata.temp.e_zmin,'String'));
    vdata.data.region.zmax = str2num(get(vdata.temp.e_zmax,'String'));
    
    vdata.data.exportobj.xunit = str2num(get(e8,'String'));
    vdata.data.exportobj.yunit = str2num(get(e9,'String'));
    vdata.data.exportobj.zunit = str2num(get(e10,'String'));
    vdata.data.exportobj.xscale = str2num(get(e11,'String'));
    vdata.data.exportobj.yscale = str2num(get(e12,'String'));
    vdata.data.exportobj.zscale = str2num(get(e13,'String'));
    vdata.data.exportobj.outputoffsetx = str2num(get(e14,'String'));
    vdata.data.exportobj.outputoffsety = str2num(get(e15,'String'));
    vdata.data.exportobj.outputoffsetz = str2num(get(e16,'String'));
    vdata.data.exportobj.blocksizex = str2num(get(e17,'String'));
    vdata.data.exportobj.blocksizey = str2num(get(e18,'String'));
    vdata.data.exportobj.blocksizez = str2num(get(e19,'String'));
    
    vdata.data.exportobj.extractwhich=get(pmh2,'value');
    vdata.data.exportobj.exportmodestring=get(pmh2,'string');
    vdata.data.exportobj.exportmodestring=vdata.data.exportobj.exportmodestring{vdata.data.exportobj.extractwhich};
    vdata.data.exportobj.objectcolors=get(pmh3,'value');
    vdata.data.exportobj.targetfileprefix=get(e20,'String');
    vdata.data.exportobj.targetfolder=get(vdata.temp.e21,'String');
    vdata.data.exportobj.includefoldernames = get(c1,'value');
    vdata.data.exportobj.invertz = get(c2,'value');
    vdata.data.exportobj.closesurfaces = get(c3,'value');
    vdata.data.exportobj.write3dsmaxloader = get(c4,'value');
    vdata.data.exportobj.skipmodelgeneration = get(c5,'value');
    
    vdata.data.exportobj.savesurfacestats = get(c6,'value');
    vdata.data.exportobj.surfacestatsfile = get(e21,'String');
  end;
  
  if (vdata.ui.temp.closefig==1) %to distinguish close on button press and close on window x
    close(f);
  end;

  if (vdata.state.lastcancel==0)
    
    if (vdata.data.exportobj.objectcolors==2)
      if (~isfield(vdata.data,'measurevol'))
        warndlg('ERROR: To use segment volume colors, please compute volumes first ("Measure / Measure Segment Volumes" in the main menu)!','Export 3D Surfaces as OBJ Files');
        releasegui();
        return;
      end;
      if (~isfield(vdata.data.measurevol,'lastvolume'))
        warndlg('ERROR: To use segment volume colors, please compute volumes first ("Measure / Measure Segment Volumes" in the main menu)!','Export 3D Surfaces as OBJ Files');
        releasegui();
        return;
      end;
    end;
    
    if ((vdata.data.exportobj.targetfolder(end)~='/')&&(vdata.data.exportobj.targetfolder(end)~='\'))
      vdata.data.exportobj.targetfolder=[vdata.data.exportobj.targetfolder '/'];
    end;
    
    if ((vdata.data.exportobj.xunit==0)||(vdata.data.exportobj.yunit==0)||(vdata.data.exportobj.zunit==0))
      res = questdlg(sprintf('Warning: The voxel size is set to (%f,%f,%f) which will result in collapsed models. Are you sure you want to continue?',vdata.data.exportobj.xunit,vdata.data.exportobj.yunit,vdata.data.exportobj.zunit),'Export 3D Surfaces as OBJ Files','Yes','No','Yes');
      if strcmp(res,'No')
        releasegui();
        return; 
      end
    end;
    
    extractsurfaces();
  end;  
  releasegui();  
  
  
function [] = callback_region_settofull(varargin)
  global vdata;
  vinfo=vdata.vast.getinfo();
  vdata.data.region.xmin=0;
  vdata.data.region.xmax=vinfo.datasizex-1;
  vdata.data.region.ymin=0;
  vdata.data.region.ymax=vinfo.datasizey-1;
  vdata.data.region.zmin=0; %first slice
  vdata.data.region.zmax=vinfo.datasizez-1; %last slice
  set(vdata.temp.e_xmin,'String',sprintf('%d', vdata.data.region.xmin));
  set(vdata.temp.e_xmax,'String',sprintf('%d', vdata.data.region.xmax));
  set(vdata.temp.e_ymin,'String',sprintf('%d', vdata.data.region.ymin));
  set(vdata.temp.e_ymax,'String',sprintf('%d', vdata.data.region.ymax));
  set(vdata.temp.e_zmin,'String',sprintf('%d', vdata.data.region.zmin));
  set(vdata.temp.e_zmax,'String',sprintf('%d', vdata.data.region.zmax));
  if (varargin{3}==1) callback_update_targetimagesize(); end;
  
  
function [] = callback_region_settocurrentvoxel(varargin)
  global vdata;
  vinfo=vdata.vast.getinfo();
  vdata.data.region.xmin=vinfo.currentviewx;
  vdata.data.region.xmax=vinfo.currentviewx;
  vdata.data.region.ymin=vinfo.currentviewy;
  vdata.data.region.ymax=vinfo.currentviewy;
  vdata.data.region.zmin=vinfo.currentviewz; %first slice
  vdata.data.region.zmax=vinfo.currentviewz; %last slice
  set(vdata.temp.e_xmin,'String',sprintf('%d', vdata.data.region.xmin));
  set(vdata.temp.e_xmax,'String',sprintf('%d', vdata.data.region.xmax));
  set(vdata.temp.e_ymin,'String',sprintf('%d', vdata.data.region.ymin));
  set(vdata.temp.e_ymax,'String',sprintf('%d', vdata.data.region.ymax));
  set(vdata.temp.e_zmin,'String',sprintf('%d', vdata.data.region.zmin));
  set(vdata.temp.e_zmax,'String',sprintf('%d', vdata.data.region.zmax));
  if (varargin{3}==1) callback_update_targetimagesize(); end;
  
  
function [] = callback_region_extendtocurrentvoxel(varargin)
  global vdata;
  vinfo=vdata.vast.getinfo();
  vdata.data.region.xmin=min([vdata.data.region.xmin vinfo.currentviewx]);
  vdata.data.region.xmax=max([vdata.data.region.xmax vinfo.currentviewx]);
  vdata.data.region.ymin=min([vdata.data.region.ymin vinfo.currentviewy]);
  vdata.data.region.ymax=max([vdata.data.region.ymax vinfo.currentviewy]);
  vdata.data.region.zmin=min([vdata.data.region.zmin vinfo.currentviewz]);
  vdata.data.region.zmax=max([vdata.data.region.zmax vinfo.currentviewz]);
  set(vdata.temp.e_xmin,'String',sprintf('%d', vdata.data.region.xmin));
  set(vdata.temp.e_xmax,'String',sprintf('%d', vdata.data.region.xmax));
  set(vdata.temp.e_ymin,'String',sprintf('%d', vdata.data.region.ymin));
  set(vdata.temp.e_ymax,'String',sprintf('%d', vdata.data.region.ymax));
  set(vdata.temp.e_zmin,'String',sprintf('%d', vdata.data.region.zmin));
  set(vdata.temp.e_zmax,'String',sprintf('%d', vdata.data.region.zmax));
  if (varargin{3}==1) callback_update_targetimagesize(); end;

  
function [] = callback_exportobj_browse(varargin)
  global vdata;
  foldername = uigetdir(vdata.data.exportobj.targetfolder,'VastTools - Select target folder for OBJ files:');
  if (foldername~=0)
    set(vdata.temp.e21,'String',foldername);
    vdata.data.exportobj.targetfolder=foldername;
  end;
  
  
function [] = extractsurfaces()
  global vdata;
  
  set(vdata.ui.cancelbutton,'Enable','on');
  set(vdata.ui.message,'String',{'Exporting Surfaces ...','Loading Metadata ...'});
  pause(0.1);
  
  param=vdata.data.exportobj;
  rparam=vdata.data.region;
  [data,res] = vdata.vast.getallsegmentdatamatrix();
  [name,res] = vdata.vast.getallsegmentnames();
  seglayername=getseglayername();
  name(1)=[]; %remove 'Background'
  maxobjectnumber=max(data(:,1));
  
  xmin=bitshift(rparam.xmin,-param.miplevel);
  xmax=bitshift(rparam.xmax,-param.miplevel)-1;
  ymin=bitshift(rparam.ymin,-param.miplevel);
  ymax=bitshift(rparam.ymax,-param.miplevel)-1;
  zmin=rparam.zmin;
  zmax=rparam.zmax;
  
  mipfact=bitshift(1,param.miplevel);
  
  if (((xmin==xmax)||(ymin==ymax)||(zmin==zmax))&&(vdata.data.exportobj.closesurfaces==0))
    warndlg('ERROR: The Matlab surface script needs a volume which is at least two pixels wide in each direction to work. Please adjust "Render from area" values, or enable "Close surface sides".','VastTools OBJ Exporting');
    set(vdata.ui.message,'String','Canceled.');
    set(vdata.ui.cancelbutton,'Enable','off');
    vdata.state.lastcancel=0;
    return;
  end;
  
  % Compute full name (including folder names) from name and hierarchy
  if (param.includefoldernames==1)
    fullname=name;
    for i=1:1:size(data,1)
      j=i;
      while data(j,14)~=0 %Check if parent is not 0
        j=data(j,14);
        fullname{i}=[name{j} '.' fullname{i}];
      end;
    end;
    name=fullname;
  end;
  
  % Compute list of objects to export
  switch param.extractwhich
    case 1  %All segments individually, uncollapsed
      objects=uint32([data(:,1) data(:,2)]); 
      vdata.vast.setsegtranslation([],[]);

    case 2  %All segments, collapsed as in Vast
      %4: Collapse segments as in the view during segment text file exporting
      objects=unique(data(:,18));
      objects=uint32([objects data(objects,2)]);
      vdata.vast.setsegtranslation(data(:,1),data(:,18));
      
    case 3  %Selected segment and children, uncollapsed
      selected=find(bitand(data(:,2),65536)>0);
      if (min(size(selected))==0)
        objects=uint32([data(:,1) data(:,2)]); 
      else
        selected=[selected getchildtreeids(data,selected)];
        objects=uint32([selected' data(selected,2)]);
      end;
      vdata.vast.setsegtranslation(data(selected,1),data(selected,1));
      
    case 4  %Selected segment and children, collapsed as in Vast
      selected=find(bitand(data(:,2),65536)>0);
      if (min(size(selected))==0)
        %None selected: choose all, collapsed
        selected=data(:,1);
        objects=unique(data(:,18));
      else
        selected=[selected getchildtreeids(data,selected)];
        objects=unique(data(selected,18));
      end;

      objects=uint32([objects data(objects,2)]);
      vdata.vast.setsegtranslation(data(selected,1),data(selected,18));
  end;
  
  
  % Compute number of blocks in volume
  nrxtiles=0; tilex1=xmin;
  while (tilex1<=xmax)
    tilex1=tilex1+param.blocksizex-param.overlap;
    nrxtiles=nrxtiles+1;
  end;
  nrytiles=0; tiley1=ymin;
  while (tiley1<=ymax)
    tiley1=tiley1+param.blocksizey-param.overlap;
    nrytiles=nrytiles+1;
  end;
  nrztiles=0; tilez1=zmin;
  if (vdata.data.exportobj.slicestep==1)
    slicenumbers=zmin:zmax;
    while (tilez1<=zmax)
      tilez1=tilez1+param.blocksizez-param.overlap;
      nrztiles=nrztiles+1;
    end;
  else
    slicenumbers=zmin:vdata.data.exportobj.slicestep:zmax;
    nrztiles=ceil(size(slicenumbers,2)/(param.blocksizez-param.overlap));
    j=1;
    for p=1:param.blocksizez-param.overlap:size(slicenumbers,2)
      pe=min([p+param.blocksizez-1 size(slicenumbers,2)]);
      blockslicenumbers{j}=slicenumbers(p:pe);
      j=j+1;
    end;
  end;
  param.nrxtiles=nrxtiles; param.nrytiles=nrytiles; param.nrztiles=nrztiles;
  
  
  %Go through all blocks and extract surfaces
  param.farray=cell(maxobjectnumber,param.nrxtiles,param.nrytiles,param.nrztiles);
  param.varray=cell(maxobjectnumber,param.nrxtiles,param.nrytiles,param.nrztiles);
  
  param.objects=objects;
  param.objectvolume=zeros(size(objects,1),1);

  tilez1=zmin; tz=1;
  while ((tz<=nrztiles)&&(vdata.state.lastcancel==0))
    tilez2=tilez1+param.blocksizez-1;
    if (tilez2>zmax) tilez2=zmax; end;
    tilezs=tilez2-tilez1+1;
    tiley1=ymin; ty=1;
    while ((ty<=nrytiles)&&(vdata.state.lastcancel==0))
      tiley2=tiley1+param.blocksizey-1;
      if (tiley2>ymax) tiley2=ymax; end;
      tileys=tiley2-tiley1+1;
      tilex1=xmin; tx=1;
      while ((tx<=nrxtiles)&&(vdata.state.lastcancel==0))
        tilex2=tilex1+param.blocksizex-1;
        if (tilex2>xmax) tilex2=xmax; end;
        tilexs=tilex2-tilex1+1;
        
        message={'Exporting Surfaces ...',sprintf('Loading Segmentation Cube (%d,%d,%d) of (%d,%d,%d)...',tx,ty,tz,nrxtiles,nrytiles,nrztiles)};
        set(vdata.ui.message,'String',message);
        pause(0.01);
        %Read this cube
        if (vdata.data.exportobj.slicestep==1)
          [segimage,values,numbers,bboxes,res] = vdata.vast.getsegimageRLEdecodedbboxes(param.miplevel,tilex1,tilex2,tiley1,tiley2,tilez1,tilez2,0);
        else
          bs=blockslicenumbers{tz};
          segimage=uint16(zeros(tilex2-tilex1+1,tiley2-tiley1+1,size(bs,2)));
          numarr=int32(zeros(maxobjectnumber,1));
          bboxarr=zeros(maxobjectnumber,6)-1;
          firstblockslice=bs(1);
          for i=1:1:size(bs,2)
            [ssegimage,svalues,snumbers,sbboxes,res] = vdata.vast.getsegimageRLEdecodedbboxes(param.miplevel,tilex1,tilex2,tiley1,tiley2,bs(i),bs(i),0);
            segimage(:,:,i)=ssegimage;
            snumbers(svalues==0)=[];
            sbboxes(svalues==0,:)=[];
            sbboxes(:,[3 6])=sbboxes(:,[3 6])+i-1;
            svalues(svalues==0)=[];
            if (min(size(svalues))>0)
              numarr(svalues)=numarr(svalues)+snumbers;
              bboxarr(svalues,:)=vdata.vast.expandboundingboxes(bboxarr(svalues,:),sbboxes);
            end;
          end;
          values=find(numarr>0);
          numbers=numarr(values);
          bboxes=bboxarr(values,:);
        end;
        
        message={'Exporting Surfaces ...',sprintf('Processing Segmentation Cube (%d,%d,%d) of (%d,%d,%d)...',tx,ty,tz,nrxtiles,nrytiles,nrztiles)};
        set(vdata.ui.message,'String',message);
        pause(0.01);        
        
        numbers(values==0)=[];
        bboxes(values==0,:)=[];
        values(values==0)=[];
        
        if (min(size(values))>0)
          % VAST now translates the voxel data before transmission because Matlab is too slow.

          %Close surfaces
          xvofs=0; yvofs=0; zvofs=0; ttxs=tilexs; ttys=tileys; ttzs=tilezs;
          if (vdata.data.exportobj.closesurfaces==1)
            if (tx==1)
              segimage=cat(1,zeros(1,size(segimage,2),size(segimage,3)),segimage);
              bboxes(:,1)=bboxes(:,1)+1;
              bboxes(:,4)=bboxes(:,4)+1;
              xvofs=-1; 
              ttxs=ttxs+1;
            end;
            if (ty==1)
              segimage=cat(2,zeros(size(segimage,1),1,size(segimage,3)),segimage);
              bboxes(:,2)=bboxes(:,2)+1;
              bboxes(:,5)=bboxes(:,5)+1;
              yvofs=-1;
              ttys=ttys+1;
            end;
            if (tz==1)
              segimage=cat(3,zeros(size(segimage,1),size(segimage,2),1),segimage);
              bboxes(:,3)=bboxes(:,3)+1;
              bboxes(:,6)=bboxes(:,6)+1;
              zvofs=-1;
              ttzs=ttzs+1;
            end;
            if (tx==nrxtiles)
              segimage=cat(1,segimage,zeros(1,size(segimage,2),size(segimage,3)));
              ttxs=ttxs+1;
            end;
            if (ty==nrytiles)
              segimage=cat(2,segimage,zeros(size(segimage,1),1,size(segimage,3)));
              ttys=ttys+1;
            end;
            if (tz==nrztiles)
              segimage=cat(3,segimage,zeros(size(segimage,1),size(segimage,2),1));
              ttzs=ttzs+1;
            end;
          end;
          
          %Extract all segments
          segnr=1;
          while ((segnr<=size(values,1))&&(vdata.state.lastcancel==0))
            seg=values(segnr);

            if (mod(segnr,10)==1)
              set(vdata.ui.message,'String',[message sprintf('Objects %d-%d of %d ...',segnr,min([segnr+9 size(values,1)]),size(values,1))]);
              pause(0.01);
            end;
            
            bbx=bboxes(segnr,:);
            bbx=bbx+[-1 -1 -1 1 1 1];
            if (bbx(1)<1) bbx(1)=1; end;
            if (bbx(2)<1) bbx(2)=1; end;
            if (bbx(3)<1) bbx(3)=1; end;
            if (bbx(4)>ttxs) bbx(4)=ttxs; end;
            if (bbx(5)>ttys) bbx(5)=ttys; end;
            if (bbx(6)>ttzs) bbx(6)=ttzs; end;

            %Adjust extracted subvolumes to be at least 2 pixels in each direction
            if (bbx(1)==bbx(4))
              if (bbx(1)>1)
                bbx(1)=bbx(1)-1;
              else
                bbx(4)=bbx(4)+1;
              end;
            end;
            if (bbx(2)==bbx(5))
              if (bbx(2)>1)
                bbx(2)=bbx(2)-1;
              else
                bbx(5)=bbx(5)+1;
              end;
            end;
            if (bbx(3)==bbx(6))
              if (bbx(3)>1)
                bbx(3)=bbx(3)-1;
              else
                bbx(6)=bbx(6)+1;
              end;
            end;
            
            subseg=segimage(bbx(1):bbx(4),bbx(2):bbx(5),bbx(3):bbx(6)); %(ymin:ymax,xmin:xmax,zmin:zmax);
            subseg=double(subseg==seg);
            
            [f,v]=isosurface(subseg,0.5);
            if (size(v,1)>0)
              %adjust coordinates for bbox and when we added empty slices at beginning
              v(:,1)=v(:,1)+bbx(2)-1+yvofs;
              v(:,2)=v(:,2)+bbx(1)-1+xvofs;
              v(:,3)=v(:,3)+bbx(3)-1+zvofs;
              
              v(:,1)=v(:,1)+tiley1-1;
              v(:,2)=v(:,2)+tilex1-1;
              if (vdata.data.exportobj.slicestep==1)
                v(:,3)=v(:,3)+tilez1-1;
              else
                v(:,3)=((v(:,3)-0.5)*vdata.data.exportobj.slicestep)+0.5+firstblockslice-1;
              end;
              v(:,1)=v(:,1)*param.yscale*param.yunit*mipfact;
              v(:,2)=v(:,2)*param.xscale*param.xunit*mipfact;
              v(:,3)=v(:,3)*param.zscale*param.zunit;
            end;
            param.farray{seg,tx,ty,tz}=f;
            param.varray{seg,tx,ty,tz}=v;
            
            segnr=segnr+1;
          end;
        end;
        
        tilex1=tilex1+param.blocksizex-param.overlap;
        tx=tx+1;
      end;
      tiley1=tiley1+param.blocksizey-param.overlap;
      ty=ty+1;
    end;
    tilez1=tilez1+param.blocksizez-param.overlap;
    tz=tz+1;
  end;
  
  vdata.vast.setsegtranslation([],[]);
  
  if (vdata.state.lastcancel==0)
    message={'Exporting Surfaces ...', 'Merging meshes...'};
    set(vdata.ui.message,'String',message);
    pause(0.01);
    
    param.objectsurfacearea=zeros(size(objects,1),1);
    
    switch vdata.data.exportobj.objectcolors
      case 1  %actual object colors
        colors=zeros(size(param.objects,1),3);
        for segnr=1:1:size(param.objects,1)
          seg=param.objects(segnr,1);
          %Get color from where the color is currently inherited from
          inheritseg=data(seg,18);
          colors(seg,:)=data(inheritseg, 3:5);
        end;
      case 2  %colors from volume
        j=jet(256);
        vols=1+255*vdata.data.measurevol.lastvolume/max(vdata.data.measurevol.lastvolume);
        cols=j(round(vols),:);
        objs=vdata.data.measurevol.lastobjects(:,1);
        colors=zeros(size(param.objects,1),3); %vcols=zeros(nro,3);
        colors(objs,:)=cols*255;
    end;


    %Write 3dsmax bulk loader script
    if (vdata.data.exportobj.write3dsmaxloader==1)
      save3dsmaxloader([param.targetfolder 'loadallobj_here.ms']);
    end;
    
    %Merge full objects from components
    segnr=1;
    while ((segnr<=size(param.objects,1))&&(vdata.state.lastcancel==0))
      seg=param.objects(segnr,1);
      set(vdata.ui.message,'String',{'Exporting Surfaces ...', ['Merging parts of ' name{seg} '...']});
      pause(0.01);
      
      cofp=[];
      covp=[];
      vofs=0;

      for z=1:1:param.nrztiles
        for y=1:1:param.nrytiles
          for x=1:1:param.nrxtiles
            if (x==1)
              f=param.farray{seg,x,y,z};
              v=param.varray{seg,x,y,z};
            else
              %disp(sprintf('Merging object %d, cube (%d,%d,%d)...',seg,x,y,z));
              [f,v]=mergemeshes(f,v,param.farray{seg,x,y,z},param.varray{seg,x,y,z});
            end;
          end;
          if (y==1)
            fc=f;
            vc=v;
          else
            %disp(sprintf('Merging object %d, row (%d,%d)...',seg,y,z));
            [fc,vc]=mergemeshes(fc,vc,f,v);
          end;
        end;
        if (z==1)
          fp=fc;
          vp=vc;
        else
          %disp(sprintf('Merging object %d, plane %d...',seg,z));
          [fp,vp]=mergemeshes(fp,vp,fc,vc);
          
          %Take out non-overlapping part of matrices to speed up computation
          if ((size(vp,1)>1)&&(size(fp,1)>1))
            vcut=find(vp(:,3)==max(vp(:,3)),1,'first')-1;
            fcutind=find(fp>vcut,1,'first');
            [fcut,j]=ind2sub(size(fp),fcutind); fcut=fcut-1;
          
            covp=[covp; vp(1:vcut,:)]; vp=vp(vcut+1:end,:);
            ovofs=vofs;
            vofs=vofs+vcut;
            cofp=[cofp; fp(1:fcut,:)+ovofs]; fp=fp(fcut+1:end,:)-vcut;
          end;
        end;
      end;
      
      vp=[covp; vp];
      fp=[cofp; fp+vofs];

      %invert Z axis if requested
      if (vdata.data.exportobj.invertz==1)
        if (size(vp,1)>0)
          vp(:,3)=-vp(:,3);
        end;
      end;
      
      on=name{find(data(:,1)==seg)};
      on(on==' ')='_';
      on(on=='?')='_';
      on(on=='*')='_';
      on(on=='\')='_';
      on(on=='/')='_';
      on(on=='|')='_';
      on(on==':')='_';
      on(on=='"')='_';
      on(on=='<')='_';
      on(on=='>')='_';
      filename=[param.targetfolder param.targetfileprefix sprintf('_%04d_%s.obj',seg,on)];

      if ((vdata.data.exportobj.skipmodelgeneration==0)&&(max(size(vp))>0))
        objectname=[param.targetfileprefix sprintf('_%04d_%s',seg,name{seg})];
        mtlfilename=[param.targetfileprefix sprintf('_%04d_%s.mtl',seg,on)];
        mtlfilenamewithpath=[filename(1:end-3) 'mtl'];
        materialname=[param.targetfileprefix sprintf('_%04d_material',seg)];

        set(vdata.ui.message,'String',{'Exporting Surfaces ...', ['Saving ' filename ' as Wavefront OBJ.....']});
        pause(0.01);

        if (vdata.data.exportobj.invertz==1)
          vertface2obj_mtllink(vp,fp,filename,objectname,mtlfilename,materialname);
        else
          vertface2obj_mtllink_invnormal(vp,fp,filename,objectname,mtlfilename,materialname);
        end;

        savematerialfile(mtlfilenamewithpath,materialname,colors(seg,:)/255);
        
      else
        %disp([filename ' is empty. Not saving.']);
      end;
      
      param.vparray{seg}=vp;
      param.fparray{seg}=fp;
      
      %%%%%% Compute surface size
      if (vdata.data.exportobj.savesurfacestats==1)
        %disp('Evaluating surface area...');
        set(vdata.ui.message,'String',{'Exporting Surfaces ...', ['Evaluating surface area of ' name{seg} ' ...']});
        pause(0.01);
        if (min(size(vp))>0)
          tnr=segnr;
          for tri=1:1:size(fp,1)
            v0=vp(fp(tri,1),:);
            v1=vp(fp(tri,2),:);
            v2=vp(fp(tri,3),:);
            a=cross(v1-v0,v2-v0); %abs not necessary because the values are squared later
            param.objectsurfacearea(tnr)=param.objectsurfacearea(tnr)+sqrt(sum(a.*a))/2;
          end;
        end;
      end;
      
      segnr=segnr+1;
    end;
  end;
  
  if ((vdata.state.lastcancel==0)&&(vdata.data.exportobj.savesurfacestats==1))
    %write surface area values to text file
    fid = fopen([param.targetfolder vdata.data.exportobj.surfacestatsfile], 'wt');
    if (fid>0)
      fprintf(fid,'%% VastTools Surface Area Export\r\n');
      fprintf(fid,'%% Provided as-is, no guarantee for correctness!\r\n');
      fprintf(fid,'%% %s\r\n\r\n',get(vdata.fh,'name'));
      
      fprintf(fid,'%% Source File: %s\r\n',seglayername);
      fprintf(fid,'%% Mode: %s\r\n', vdata.data.exportobj.exportmodestring);
      fprintf(fid,'%% Area: (%d-%d, %d-%d, %d-%d)\r\n',rparam.xmin,rparam.xmax,rparam.ymin,rparam.ymax,rparam.zmin,rparam.zmax);
      fprintf(fid,'%% Computed at voxel size: (%f,%f,%f)\r\n',param.xscale*param.xunit*mipfact,param.yscale*param.yunit*mipfact,param.zscale*param.zunit*vdata.data.exportobj.slicestep);
      fprintf(fid,'%% Columns are: Object Name, Object ID, Surface Area in Export\r\n\r\n');
      for segnr=1:1:size(param.objects,1)
        seg=param.objects(segnr,1);
        fprintf(fid,'"%s"  %d  %f\r\n',name{seg},seg,param.objectsurfacearea(segnr));
      end;
      fprintf(fid,'\r\n');
      fclose(fid);
    end;
  end;
  
  if (vdata.state.lastcancel==0)
    set(vdata.ui.message,'String','Done.');
  else
    set(vdata.ui.message,'String','Canceled.');
  end;
  set(vdata.ui.cancelbutton,'Enable','off');
  vdata.state.lastcancel=0;
  

function name=getseglayername()
  global vdata;
  nrl=vdata.vast.getnroflayers();
  found=0; l=0; name=[];
  while ((found==0)&&(l<nrl))
    linf=vdata.vast.getlayerinfo(l);
    if (linf.type==1)
      name=linf.name;
      found=1;
    end;
    l=l+1;
  end;
  

function [f,v]=mergemeshes(f1,v1,f2,v2)
  %mergemeshes.m
  %merges meshes defined by (f1,v1) and (f2,v2)
  %by Daniel Berger for MIT-BCS Seung, August 2011
  
  if (min(size(v1)))==0
    f=f2;
    v=v2;
    return;
  end;
  
  if (min(size(v2)))==0
    f=f1;
    v=v1;
    return;
  end;

  nrofvertices1=size(v1,1);
  nrofvertices2=size(v2,1);
  f2=f2+nrofvertices1;
  
  %find overlapping vertex region
  minv1=min(v1);
  maxv1=max(v1);
  minv2=min(v2);
  maxv2=max(v2);
  ovmin=max(minv1,minv2);
  ovmax=min(maxv1,maxv2);
  
  ov1=[(1:size(v1,1))' v1];
  ov1=ov1(((ov1(:,2)>=ovmin(1))&(ov1(:,2)<=ovmax(1))&(ov1(:,3)>=ovmin(2))&(ov1(:,3)<=ovmax(2))&(ov1(:,4)>=ovmin(3))&(ov1(:,4)<=ovmax(3))),:);
  ov2=[(1:size(v2,1))' v2];
  ov2=ov2(((ov2(:,2)>=ovmin(1))&(ov2(:,2)<=ovmax(1))&(ov2(:,3)>=ovmin(2))&(ov2(:,3)<=ovmax(2))&(ov2(:,4)>=ovmin(3))&(ov2(:,4)<=ovmax(3))),:);
  
  if (min(size(ov2))==0)
    %Non-overlapping objects
    f=[f1; f2];
    v=[v1; v2];
    return;
  end;

  %Link all vertices in v2 which have corresponding vertices in v1 to the v1 vertex
  deletevertex=zeros(nrofvertices2,1);
  facetouched=zeros(size(f2,1),3);
  oldcomparison=0;
  if (oldcomparison==1)
    for oi=1:1:size(ov1,1) %nrofvertices1
      i=ov1(oi,1);
      %This vertex is in the overlap zone. Find corresponding vertex in v2
      r2=find(ismember(ov2(:,2:4),v1(i,:),'rows'));
      i2=ov2(r2,1);
      
      %i is the index (row number) of the corresponding vertex in v1
      %i2 is the index (row number) of the corresponding vertex in v2
      %Exchange all occurences of this vertex in f2 with the vertex number used in f1
      facetouched(f2==i2+nrofvertices1)=1;
      f2(f2==i2+nrofvertices1)=i;
      deletevertex(i2)=1;
    end;
  else
    [c,i1a,i2a]=intersect(ov1(:,2:4),ov2(:,2:4),'rows');
    for oi=1:1:size(i1a,1)
      i=ov1(i1a(oi),1);
      i2=ov2(i2a(oi),1);
      facetouched(f2==i2+nrofvertices1)=1;
      f2(f2==i2+nrofvertices1)=i;
      deletevertex(i2)=1;
    end;
  end;

  %Delete the unused vertices in v2 and re-label the vertices in v2 accordingly
  %compute list of old and new vertex numbers
  z=[[1:nrofvertices1]'; zeros(size(deletevertex,1),1)];
  zp=nrofvertices1+1;
  for sp=1:1:size(deletevertex,1)
    if (deletevertex(sp)==0)
      z(nrofvertices1+sp)=zp;
      zp=zp+1;
    end;
  end;
  
  f2d=z(f2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Delete duplicated faces from f2d
  facetouched=max(facetouched,[],2); %Flags for faces which had corners edited; only those might be duplicated
  pf2d=[(1:size(f2d,1))' sort(f2d,2)]; %Prefix for indexing
  pf2d=pf2d(facetouched==1,:); %Pick out only those faces from f2 which were changed
  sf2d=sortrows(pf2d,[2 3 4]);
  sf1=sortrows(sort(f1,2));
  
  r=ismember(sf2d(:,2:end),sf1,'rows');
  removeface=sf2d(r==1,1);
  v2d=v2(deletevertex==0,:);
  
  f=[f1; f2d];
  v=[v1; v2d];

  
function vertface2obj_mtllink_invnormal(v,f,filename,objectname,mtlfilename,materialname)
  % VERTFACE2OBJ Save a set of vertex coordinates and faces as a Wavefront/Alias Obj file
  % VERTFACE2OBJ(v,f,fname)
  %     v is a Nx3 matrix of vertex coordinates.
  %     f is a Mx3 matrix of vertex indices.
  %     fname is the filename to save the obj file.
  
  fid = fopen(filename,'w');
  
  fprintf(fid,'mtllib %s\n',mtlfilename);
  fprintf(fid,'usemtl %s\n',materialname);
  
  for i=1:size(v,1)
    fprintf(fid,'v %f %f %f\n',v(i,1),v(i,2),v(i,3));
  end
  
  fprintf(fid,'g %s\n',objectname);
  
  for i=1:size(f,1);
    %2 1 3 order to flip normal
    fprintf(fid,'f %d %d %d\n',f(i,2),f(i,1),f(i,3));
  end
  fprintf(fid,'g\n');
  
  fclose(fid);
  
  
function vertface2obj_mtllink(v,f,filename,objectname,mtlfilename,materialname)
  % VERTFACE2OBJ Save a set of vertex coordinates and faces as a Wavefront/Alias Obj file
  % VERTFACE2OBJ(v,f,fname)
  %     v is a Nx3 matrix of vertex coordinates.
  %     f is a Mx3 matrix of vertex indices.
  %     fname is the filename to save the obj file.
  
  fid = fopen(filename,'w');
  
  fprintf(fid,'mtllib %s\n',mtlfilename);
  fprintf(fid,'usemtl %s\n',materialname);
  
  for i=1:size(v,1)
    fprintf(fid,'v %f %f %f\n',v(i,1),v(i,2),v(i,3));
  end
  
  fprintf(fid,'g %s\n',objectname);
  
  for i=1:size(f,1);
    fprintf(fid,'f %d %d %d\n',f(i,1),f(i,2),f(i,3));
  end
  fprintf(fid,'g\n');
  
  fclose(fid);
  
  
function savematerialfile(filename,materialname,color)
  %Define constant material parameters
  Ns=50.0000;
  Ni=1.5000;
  d=0.4000; %Opacity (0.0 is transparent, 1.0 is opaque)
  Tr=0.0000;
  Tf=[1.0000 1.0000 1.0000];
  illum=2;
  Ka=[0.0000 0.0000 0.0000];
  Kd=color;
  Ks=[1.0000 1.0000 1.0000];
  Ke=[0.0000 0.0000 0.0000];
  
  %disp(sprintf('Saving %s ...',filename));
  fid = fopen(filename, 'wt');
  fprintf(fid,'# MTL writer by Daniel Berger, Jan 2015\r\n');
  fprintf(fid,'newmtl %s\r\n',materialname);
  fprintf(fid,'  Ns %.4f\r\n',Ns);
  fprintf(fid,'  Ni %.4f\r\n',Ni);
  fprintf(fid,'  d %.4f\r\n',d);
  fprintf(fid,'  Tr %.4f\r\n',Tr);
  fprintf(fid,'  Tf %.4f %.4f %.4f\r\n',Tf(1),Tf(2),Tf(3));
  fprintf(fid,'  illum %d\r\n',illum);
  fprintf(fid,'  Ka %.4f %.4f %.4f\r\n',Ka(1),Ka(2),Ka(3));
  fprintf(fid,'  Kd %.4f %.4f %.4f\r\n',Kd(1),Kd(2),Kd(3));
  fprintf(fid,'  Ks %.4f %.4f %.4f\r\n',Ks(1),Ks(2),Ks(3));
  fprintf(fid,'  Ke %.4f %.4f %.4f\r\n',Ke(1),Ke(2),Ke(3));
  fprintf(fid,'\r\n');
  fclose(fid);


function save3dsmaxloader(filename)
  fid = fopen(filename, 'wt');
  fprintf(fid,'files = getfiles ".\\*.obj"\r\n');
  fprintf(fid,'for f in files do (importfile (f) #noprompt)\r\n');
  fclose(fid);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Projection Exporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = callback_exportprojection(varargin)
  global vdata;
  
  if (vdata.state.isconnected==0)
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
    return;
  end;
  
  vinfo=vdata.vast.getinfo();
  if (min(size(vinfo)==0))
    warndlg('ERROR: Requesting data from VAST failed.','Error with remote connection to VAST');
    return;
  end;
  
  if (min([vinfo.datasizex vinfo.datasizey vinfo.datasizez])==0)
    warndlg('ERROR: No volume open in VAST.','VastTools projection image exporting');
    return;
  end;
  
%  blockgui();
  
  %Display parameter dialog
  if (~isfield(vdata.data,'region'))
    vdata.data.region.xmin=0;
    vdata.data.region.xmax=vinfo.datasizex-1;
    vdata.data.region.ymin=0;
    vdata.data.region.ymax=vinfo.datasizey-1;
    vdata.data.region.zmin=0; %first slice
    vdata.data.region.zmax=vinfo.datasizez-1; %last slice
  else
    if (vdata.data.region.xmin<0) vdata.data.region.xmin=0; end;
    if (vdata.data.region.xmax>(vinfo.datasizex-1)) vdata.data.region.xmax=vinfo.datasizex-1; end;
    if (vdata.data.region.ymin<0) vdata.data.region.ymin=0; end;
    if (vdata.data.region.ymax>(vinfo.datasizey-1)) vdata.data.region.ymax=vinfo.datasizey-1; end;
    if (vdata.data.region.zmin<0) vdata.data.region.zmin=0; end; %first slice
    if (vdata.data.region.zmax>(vinfo.datasizez-1)) vdata.data.region.zmax=vinfo.datasizez-1; end;
  end;
  if (~isfield(vdata.data,'exportproj'))
    vdata.data.exportproj.miplevel=0;
    
    vdata.data.exportproj.overlap=0;
    vdata.data.exportproj.slicestep=1;
    vdata.data.exportproj.projaxis=5; %1:+X, 2:-X, 3:+Y, 4:-Y, 5:+Z, 6:-Z
    vdata.data.exportproj.stretchz=2;

    vdata.data.exportproj.savetofile=1;
    vdata.data.exportproj.showinwindow=1;
    vdata.data.exportproj.targetfilename='projection.png';
    vdata.data.exportproj.targetfoldername=pwd;
    
    vdata.data.exportproj.segpreprocess=2; %extractwhich
    vdata.data.exportproj.expandsegmentation=0; %number of pixels to expand segmentation map (negative values shrink)
    vdata.data.exportproj.blurdistance=0; %in pixels; blurs inward
    vdata.data.exportproj.imagesource=1; %1: Segmentation; 2: Selected EM layer; 3: Screenshots
    vdata.data.exportproj.opacitysource=1; %1: Segmented, 2: Unsegmented, 3: All
    vdata.data.exportproj.blendmode=1; %1:Alpha-blend, 2: Additive, 3: Maximum projection
    vdata.data.exportproj.objectopacity=1; %opacity strength 0..1
    vdata.data.exportproj.useshadows=0;
    vdata.data.exportproj.shadowcone=2;
    vdata.data.exportproj.depthattenuation=1;

    vdata.data.exportproj.finalnormalize=1;
  else
    if (vdata.data.exportproj.miplevel>(vinfo.nrofmiplevels-1)) vdata.data.exportproj.miplevel=vinfo.nrofmiplevels-1; end;
  end;
  
  nrofsegments=vdata.vast.getnumberofsegments();
  if (nrofsegments==0)
    vdata.data.exportproj.rendermode=2;
  end;
  
  scrsz = get(0,'ScreenSize');
  figheight=600;
  f = figure('units','pixels','position',[50 scrsz(4)-100-figheight 500 figheight],'menubar','none','numbertitle','off','name','VastTools - Export Projection Image','resize','off');

  vpos=figheight-40;
 
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15], 'Tag','t1','String','Render at resolution:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(vinfo.nrofmiplevels,1);
  vx=vinfo.voxelsizex;
  vy=vinfo.voxelsizey;
  vz=vinfo.voxelsizez;
  for i=1:1:vinfo.nrofmiplevels
    str{i}=sprintf('Mip %d - (%.2f nm, %.2f nm, %.2f nm) voxels',i-1,vx,vy,vz);
    vx=vx*2; vy=vy*2;
  end;
  vdata.temp.pmh = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.miplevel+1,'Position',[170 vpos 300 20],'Callback',{@callback_update_targetimagesize});
  vpos=vpos-30;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15],'String','Render from area:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos 130 20], 'String','Set to full', 'CallBack',{@callback_region_settofull,1});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-25 130 20], 'String','Set to current voxel', 'CallBack',{@callback_region_settocurrentvoxel,1});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-50 130 20], 'String','Extend to current voxel', 'CallBack',{@callback_region_extendtocurrentvoxel,1});
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','X min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d', vdata.data.region.xmin),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','X max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.xmax),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Y min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymin),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Y max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymax),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Z min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmin),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Z max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmax),'horizontalalignment','left','Callback',{@callback_update_targetimagesize});
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Use every nth slice:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_slicestep = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.exportproj.slicestep),'horizontalalignment','left');
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15], 'Tag','t1','String','Projection axis:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str={'+X (lowest X in front)','-X (highest X in front)','+Y (lowest Y in front)','-Y (highest Y in front)','+Z (lowest Z in front)','-Z (highest Z in front)'};
  vdata.temp.pmaxis = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.projaxis,'Position',[170 vpos 150 20],'Callback',{@callback_update_targetimagesize});
  str={'No stretching','Stretch Z (nearest)','Stretch Z (interpolated)'};
  vdata.temp.pmstretch = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.stretchz,'Position',[330 vpos 140 20]);
  vpos=vpos-30;

  vdata.temp.t_targetsize= uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 440 16],'backgroundcolor',[0.75 0.75 0.65],'horizontalalignment','left');
  callback_update_targetimagesize();
  vpos=vpos-40;
  
  c1 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 150 15],'Value',vdata.data.exportproj.savetofile,'string','Save to file','backgroundcolor',get(f,'color')); 
  c2 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[170 vpos 200 15],'Value',vdata.data.exportproj.showinwindow,'string','Show render progress','backgroundcolor',get(f,'color')); 

  vpos=vpos-30;
    
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','File name:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  e20 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[120 vpos 290 20],'String',vdata.data.exportproj.targetfilename,'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Target folder:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e21 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[120 vpos 290 20],'String',vdata.data.exportproj.targetfoldername,'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[420 vpos 60 20], 'String','Browse...', 'CallBack',{@callback_exportproj_browse});
  vpos=vpos-40;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15], 'Tag','t1','String','Segmentation preprocessing:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(4,1);
  str{1}='All segments individually, uncollapsed';
  str{2}='All segments, collapsed as in VAST';
  str{3}='Selected segment and children, uncollapsed';
  str{4}='Selected segment and children, collapsed as in VAST';
  pmh2 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.segpreprocess,'Position',[180 vpos 290 20]);
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Expand segments by n pixels:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_expandseg = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[190 vpos 60 20],'String',sprintf('%d', vdata.data.exportproj.expandsegmentation),'horizontalalignment','left');

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[280 vpos 120 15],'String','Blur edges by n pixels:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_blurdist = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[410 vpos 60 20],'String',sprintf('%d', vdata.data.exportproj.blurdistance),'horizontalalignment','left');
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Image source:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(4,1);
  str{1}='Segmentation';
  str{2}='Selected EM layer';
  str{3}='Screenshots';
  str{4}='Seg Volume Colors';

  pmh3 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.imagesource,'Position',[110 vpos 120 20]);

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[250 vpos 120 15],'String','Opacity source:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(3,1);
  str{1}='Segmented areas';
  str{2}='Unsegmented areas';
  str{3}='All';
  pmh4 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.opacitysource,'Position',[340 vpos 130 20]);
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Blending mode:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(3,1);
  str{1}='Alpha blending';
  str{2}='Additive';
  str{3}='Max projection';
  pmh5 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.exportproj.blendmode,'Position',[120 vpos 110 20]);

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[290 vpos 120 15],'String','Object opacity [0..1]:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_objectopacity = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[410 vpos 60 20],'String',sprintf('%f', vdata.data.exportproj.objectopacity),'horizontalalignment','left');
  vpos=vpos-30;
  
  c3 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 100 15],'Value',vdata.data.exportproj.useshadows,'string','Use shadows','backgroundcolor',get(f,'color')); 
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[150 vpos-1 160 15],'String','Shadow cone angle (pix/slice):','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_shadowcone = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[310 vpos-2 100 20],'String',sprintf('%f', vdata.data.exportproj.shadowcone),'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 200 15],'String','Depth attenuation (far brightness) [0..1]:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_depthattenuation = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[240 vpos 60 20],'String',sprintf('%f', vdata.data.exportproj.depthattenuation),'horizontalalignment','left');
  vpos=vpos-30;
  c4 = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[30 vpos 250 15],'Value',vdata.data.exportproj.finalnormalize,'string','Normalize projection image','backgroundcolor',get(f,'color')); 
  vpos=vpos-30;
  
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[150 20 60 20], 'String','OK', 'CallBack',{@callback_done});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[290 20 60 20], 'String','Cancel', 'CallBack',{@callback_canceled});

  vdata.state.lastcancel=1;
  vdata.ui.temp.closefig=0;
  uiwait(f);
  
  if (vdata.state.lastcancel==0)
    vdata.data.exportproj.miplevel=get(vdata.temp.pmh,'value')-1;
    
    vdata.data.region.xmin = str2num(get(vdata.temp.e_xmin,'String'));
    vdata.data.region.xmax = str2num(get(vdata.temp.e_xmax,'String'));
    vdata.data.region.ymin = str2num(get(vdata.temp.e_ymin,'String'));
    vdata.data.region.ymax = str2num(get(vdata.temp.e_ymax,'String'));
    vdata.data.region.zmin = str2num(get(vdata.temp.e_zmin,'String'));
    vdata.data.region.zmax = str2num(get(vdata.temp.e_zmax,'String'));
    
    vdata.data.exportproj.slicestep = str2num(get(vdata.temp.e_slicestep,'String'));
    vdata.data.exportproj.projaxis = get(vdata.temp.pmaxis,'value');
    vdata.data.exportproj.stretchz = get(vdata.temp.pmstretch,'value'); %get(c1,'value');
    vdata.data.exportproj.savetofile = get(c1,'value');
    vdata.data.exportproj.showinwindow = get(c2,'value');
    vdata.data.exportproj.targetfilename=get(e20,'String');
    vdata.data.exportproj.targetfoldername=get(vdata.temp.e21,'String');

    vdata.data.exportproj.segpreprocess=get(pmh2,'value');
    vdata.data.exportproj.expandsegmentation = str2num(get(vdata.temp.e_expandseg,'String'));
    vdata.data.exportproj.blurdistance = str2num(get(vdata.temp.e_blurdist,'String'));
    vdata.data.exportproj.imagesource=get(pmh3,'value');
    vdata.data.exportproj.opacitysource=get(pmh4,'value');
    vdata.data.exportproj.blendmode=get(pmh5,'value');
    vdata.data.exportproj.objectopacity = str2num(get(vdata.temp.e_objectopacity,'String'));
    if (vdata.data.exportproj.objectopacity<0) vdata.data.exportproj.objectopacity=0; end;
    if (vdata.data.exportproj.objectopacity>1) vdata.data.exportproj.objectopacity=1; end;
    
    vdata.data.exportproj.useshadows = get(c3,'value');
    vdata.data.exportproj.shadowcone = str2num(get(vdata.temp.e_shadowcone,'String'));
   
    vdata.data.exportproj.depthattenuation = str2num(get(vdata.temp.e_depthattenuation,'String'));
    if (vdata.data.exportproj.depthattenuation<0) vdata.data.exportproj.depthattenuation=0; end;
    if (vdata.data.exportproj.depthattenuation>1) vdata.data.exportproj.depthattenuation=1; end;
    vdata.data.exportproj.finalnormalize = get(c4,'value');
  end;
  
  if (vdata.ui.temp.closefig==1) %to distinguish close on button press and close on window x
    close(f);
  end;

  if (vdata.state.lastcancel==0)

    if ((nrofsegments==0) && (vdata.data.exportproj.imagesource~=2))
      warndlg('ERROR: No segmentation available in VAST. Cannot use segmentation during exporting!','VastTools projection image exporting');
      releasegui();
      return;
    end;
    
    if (vdata.data.exportproj.imagesource==4)
      if (~isfield(vdata.data,'measurevol'))
        warndlg('ERROR: To use segment volume colors, please compute volumes first ("Measure / Measure Segment Volumes" in the main menu)!','VastTools projection image exporting');
        releasegui();
        return;
      end;
      if (~isfield(vdata.data.measurevol,'lastvolume'))
        warndlg('ERROR: To use segment volume colors, please compute volumes first ("Measure / Measure Segment Volumes" in the main menu)!','VastTools projection image exporting');
        releasegui();
        return;
      end;
    end;
    
    %reevaluate target image size
    xmin=bitshift(vdata.data.region.xmin,-vdata.data.exportproj.miplevel);
    xmax=bitshift(vdata.data.region.xmax,-vdata.data.exportproj.miplevel);
    ymin=bitshift(vdata.data.region.ymin,-vdata.data.exportproj.miplevel);
    ymax=bitshift(vdata.data.region.ymax,-vdata.data.exportproj.miplevel);
    zval=vdata.data.region.zmin:vdata.data.exportproj.slicestep:vdata.data.region.zmax;
    switch vdata.data.exportproj.projaxis
      case {1,2} %X
        timagewidth=(ymax-ymin+1);
        timageheight=max(size(zval));
      case {3,4} %Y
        timagewidth=(xmax-xmin+1);
        timageheight=max(size(zval));
      case {5,6} %Z
        timagewidth=(xmax-xmin+1);
        timageheight=(ymax-ymin+1);
    end;
    
    if ((timagewidth>2000)||(timageheight>2000))
      res = questdlg(sprintf('With these settings you will render an image which is %d by %d pixels large. Are you sure?',timagewidth, timageheight),'VastTools projection image exporting','Yes','No','Yes');
      if strcmp(res,'No')
        releasegui();
        return; 
      end;
    end;

    renderprojection();
  end;  
  releasegui();  

  
function [] = callback_update_targetimagesize(varargin)
  global vdata;
  
  vdata.data.exportproj.miplevel=get(vdata.temp.pmh,'value')-1;
  vdata.data.region.xmin = str2num(get(vdata.temp.e_xmin,'String'));
  vdata.data.region.xmax = str2num(get(vdata.temp.e_xmax,'String'));
  vdata.data.region.ymin = str2num(get(vdata.temp.e_ymin,'String'));
  vdata.data.region.ymax = str2num(get(vdata.temp.e_ymax,'String'));
  vdata.data.region.zmin = str2num(get(vdata.temp.e_zmin,'String'));
  vdata.data.region.zmax = str2num(get(vdata.temp.e_zmax,'String'));
  vdata.data.exportproj.slicestep = str2num(get(vdata.temp.e_slicestep,'String'));
  vdata.data.exportproj.projaxis = get(vdata.temp.pmaxis,'value');
  
  xmin=bitshift(vdata.data.region.xmin,-vdata.data.exportproj.miplevel);
  xmax=bitshift(vdata.data.region.xmax,-vdata.data.exportproj.miplevel);
  ymin=bitshift(vdata.data.region.ymin,-vdata.data.exportproj.miplevel);
  ymax=bitshift(vdata.data.region.ymax,-vdata.data.exportproj.miplevel);
  zval=vdata.data.region.zmin:vdata.data.exportproj.slicestep:vdata.data.region.zmax;
  switch vdata.data.exportproj.projaxis
    case {1,2} %X
      timagewidth=(ymax-ymin+1);
      timageheight=max(size(zval));
    case {3,4} %Y
      timagewidth=(xmax-xmin+1);
      timageheight=max(size(zval));
    case {5,6} %Z
      timagewidth=(xmax-xmin+1);
      timageheight=(ymax-ymin+1);
  end;

  set(vdata.temp.t_targetsize,'String',sprintf('Target image size with these settings: (%d, %d)',timagewidth,timageheight));


function [] = callback_exportproj_browse(varargin)
  global vdata;
  foldername = uigetdir(vdata.data.exportproj.targetfoldername,'VastTools - Select target folder for projection image:');
  if (foldername~=0)
    set(vdata.temp.e21,'String',foldername);
    vdata.data.exportproj.targetfoldername=foldername;
  end;
  
function [] = renderprojection()
  global vdata;
  
  set(vdata.ui.cancelbutton,'Enable','on');
  set(vdata.ui.message,'String',{'Generating Projection Image ...','Loading Metadata ...'});
  pause(0.1);
  
  %releasegui();
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Colors from Volume
  if (vdata.data.exportproj.imagesource==4)
    nro=vdata.vast.getnumberofsegments();
    j=jet(256);
    vols=1+255*vdata.data.measurevol.lastvolume/max(vdata.data.measurevol.lastvolume);
    cols=j(round(vols),:);
    objs=vdata.data.measurevol.lastobjects(:,1);
    vcols=zeros(nro,3);
    vcols(objs,:)=cols*255;
  end;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  param=vdata.data.exportproj;
  rparam=vdata.data.region;
  vinfo=vdata.vast.getinfo();
  [selectedlayernr, selectedemlayernr, selectedsegmentlayernr, res]=vdata.vast.getselectedlayernr();
  if (selectedsegmentlayernr>-1)
    [data,res] = vdata.vast.getallsegmentdatamatrix();
    [name,res] = vdata.vast.getallsegmentnames();
    name(1)=[]; %remove 'Background'
    maxobjectnumber=max(data(:,1));
  end;
  
  xmin=bitshift(rparam.xmin,-param.miplevel);
  xmax=bitshift(rparam.xmax,-param.miplevel)-1;
  ymin=bitshift(rparam.ymin,-param.miplevel);
  ymax=bitshift(rparam.ymax,-param.miplevel)-1;
  zmin=rparam.zmin;
  zmax=rparam.zmax;
  
  mipfact=bitshift(1,param.miplevel);
  
  slicestep=vdata.data.exportproj.slicestep;
  
  segpreprocess=vdata.data.exportproj.segpreprocess; %extractwhich
  collapsesegments=0; if ((segpreprocess==2)||(segpreprocess==4)) collapsesegments=1; end;
  expandsegmentation=vdata.data.exportproj.expandsegmentation; %number of pixels to expand segmentation map (negative values shrink)
  blurdistance=vdata.data.exportproj.blurdistance; %in pixels; blurs inward
  imagesource=vdata.data.exportproj.imagesource; %1: Segmentation; 2: Selected EM layer; 3: Screenshots
  opacitysource=vdata.data.exportproj.opacitysource; %1: Segmented, 2: Unsegmented, 3: All
  blendmode=vdata.data.exportproj.blendmode; %1:Alpha-blend, 2: Additive, 3: Maximum projection
  objectopacity=vdata.data.exportproj.objectopacity; %opacity strength 0..1
  useshadows=vdata.data.exportproj.useshadows;
  shadowcone=vdata.data.exportproj.shadowcone;
  depthattenuation=vdata.data.exportproj.depthattenuation;
  finalnormalize=vdata.data.exportproj.finalnormalize;
  
  if (vdata.data.exportproj.showinwindow==1)
    lfig=figure;
  end;
 
  % Compute list of objects to export
  if ((imagesource==1)||(opacitysource<3))
    switch segpreprocess
      case 1  %All segments individually, uncollapsed
        vdata.vast.setsegtranslation([],[]);
        
      case 2  %All segments, collapsed as in Vast
        %4: Collapse segments as in the view during segment text file exporting
        vdata.vast.setsegtranslation(data(:,1),data(:,18));
        
      case 3  %Selected segment and children, uncollapsed
        selected=find(bitand(data(:,2),65536)>0);
        if (min(size(selected))~=0)
          selected=[selected getchildtreeids(data,selected)];
        end;
        vdata.vast.setsegtranslation(data(selected,1),data(selected,1));
        
      case 4  %Selected segment and children, collapsed as in Vast
        selected=find(bitand(data(:,2),65536)>0);
        if (min(size(selected))==0)
          %None selected: choose all, collapsed
          selected=data(:,1);
        else
          selected=[selected getchildtreeids(data,selected)];
        end;
        vdata.vast.setsegtranslation(data(selected,1),data(selected,18));
    end;
  end;
  
  %%%%%%%%%%%%%%%%%
  %Set up mapping between source and target volumes based on projection axis
  
  slabthickness=16;
  
  swidth=xmax-xmin+1;
  sheight=ymax-ymin+1;
  sdepth=zmax-zmin+1;
  
  switch vdata.data.exportproj.projaxis
    case 1 %'+X (lowest X in front)'
      twidth=ymax-ymin+1;
      theight=zmax-zmin+1;
      tdepth=xmax-xmin+1;
      twmin=ymin;
      thmin=zmin;
      tdmin=xmin;
      tdir=1;

      projectback=[[0 0 mipfact xmin*mipfact]; [mipfact 0 0 ymin*mipfact]; [0 slicestep 0 zmin]];
      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabboxes=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1]) ymin ymax zmin zmax];
        end;
      else
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[slabc(i) slabc(i) ymin ymax zmin zmax];
        end;
      end;
    case 2 %'-X (highest X in front)'
      twidth=ymax-ymin+1;
      theight=zmax-zmin+1;
      tdepth=xmax-xmin+1;
      twmin=ymin;
      thmin=zmin;
      tdmin=xmin;
      tdir=-1;

      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabbox=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        slabstart=fliplr(slabstart);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1]) ymin ymax zmin zmax];
        end;
      else
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        slabc=fliplr(slabc);
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[slabc(i) slabc(i) ymin ymax zmin zmax];
        end;
      end;
    case 3 %'+Y (lowest Y in front)'
      twidth=xmax-xmin+1;
      theight=zmax-zmin+1;
      tdepth=ymax-ymin+1;
      twmin=xmin;
      thmin=zmin;
      tdmin=ymin;
      tdir=1;
      projectback=[[mipfact 0 0 xmin*mipfact]; [0 0 mipfact ymin*mipfact]; [0 slicestep 0 zmin]];
      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabbox=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1]) zmin zmax];
        end;
      else
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax slabc(i) slabc(i) zmin zmax];
        end;
      end;
    case 4 %'-Y (highest Y in front)'
      twidth=xmax-xmin+1;
      theight=zmax-zmin+1;
      tdepth=ymax-ymin+1;
      twmin=xmin;
      thmin=zmin;
      tdmin=ymin;
      tdir=-1;
      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabbox=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        slabstart=fliplr(slabstart);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1]) zmin zmax];
        end;
      else
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        slabc=fliplr(slabc);
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax slabc(i) slabc(i) zmin zmax];
        end;
      end;
    case 5 %'+Z (lowest Z in front)'
      twidth=xmax-xmin+1;
      theight=ymax-ymin+1;
      tdepth=zmax-zmin+1;
      twmin=xmin;
      thmin=ymin;
      tdmin=zmin;
      tdir=1;
      projectback=[[mipfact 0 0 xmin*mipfact]; [0 mipfact 0 ymin*mipfact]; [0 0 slicestep zmin]];
      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabbox=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax ymin ymax slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1])];
        end;
      else %load each slice individually
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax ymin ymax slabc(i) slabc(i)];
        end;
      end;
    case 6 %'-Z (highest Z in front)'
      twidth=xmax-xmin+1;
      theight=ymax-ymin+1;
      tdepth=zmax-zmin+1;
      twmin=xmin;
      thmin=ymin;
      tdmin=zmin;
      tdir=-1;
      projectback=[[mipfact 0 0 xmin*mipfact]; [0 -mipfact 0 ymax*mipfact]; [0 0 -slicestep zmax]];
      if (slicestep==1) %load in slabs
        nrofslabs=ceil(tdepth/slabthickness);
        sslabbox=zeros(nrofslabs,6);
        slabstart=tdmin:slabthickness:tdmin+tdepth-1;
        slabstart=fliplr(slabstart);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax ymin ymax slabstart(i) min([slabstart(i)+slabthickness-1 tdmin+tdepth-1])];
        end;
      else %load each slice individually
        slabc=tdmin:slicestep:tdmin+tdepth-1;
        slabc=fliplr(slabc);
        nrofslabs=size(slabc,2);
        sslabbox=zeros(nrofslabs,6);
        for i=1:1:nrofslabs
          sslabbbox(i,:)=[xmin xmax ymin ymax slabc(i) slabc(i)];
        end;
      end;
  end;
      
  %%%%%%%%%%%%%%%%%
  
  rtim=zeros(theight,twidth);
  gtim=zeros(theight,twidth);
  btim=zeros(theight,twidth);
  topcolor=zeros(theight,twidth);
  mask=zeros(theight,twidth);
  shadow=ones(theight,twidth);
  ttranspmap=ones(theight,twidth);
  zmap=zeros(theight,twidth)-1;
  
  depth=0;
  
  shadowK = fspecial('disk',shadowcone);
  if (blurdistance>0)
    edgeblurK = fspecial('disk',blurdistance);
  end;
  
  d=1;
  while ((d<=nrofslabs)&&(vdata.state.lastcancel==0))
   
    if ((imagesource==1)||(opacitysource<3))
      %Load segmentation
      [segimage,res] = vdata.vast.getsegimageRLEdecoded(param.miplevel,sslabbbox(d,1),sslabbbox(d,2),sslabbbox(d,3),sslabbbox(d,4),sslabbbox(d,5),sslabbbox(d,6),0);
    else
      segimage=[];
    end;
    if (imagesource==2)
      %Load selected EM layer
      [emimage,res] = vdata.vast.getemimage(selectedemlayernr,param.miplevel,sslabbbox(d,1),sslabbbox(d,2),sslabbbox(d,3),sslabbbox(d,4),sslabbbox(d,5),sslabbbox(d,6));
      if (size(size(emimage),2)==3)
        emimage=permute(emimage,[2 1 3]);
      else
        emimage=permute(emimage,[2 1 3 4]);
      end;
    else
      emimage=[];
    end;
    if (imagesource==3)
      %Load screenshots
      [scsimage,res] = vdata.vast.getscreenshotimage(param.miplevel,sslabbbox(d,1),sslabbbox(d,2),sslabbbox(d,3),sslabbbox(d,4),sslabbbox(d,5),sslabbbox(d,6),collapsesegments);
      scsimage=permute(scsimage,[2 1 3 4]);
    else
      scsimage=[];
    end;
    
    %Rotate to target orientation
    switch vdata.data.exportproj.projaxis
    case {1 2} %'+X (lowest X in front)' %'-X (highest X in front)'
      %y->x, z->y, x->z
      if (min(size(segimage))>0)
        segimage=permute(segimage,[2 3 1]);
      end;
      if (min(size(emimage))>0)
        if (size(size(emimage),2)==3)
          emimage=permute(emimage,[2 3 1]);
        else
          emimage=permute(emimage,[2 3 1 4]);
        end;
      end;
      if (min(size(scsimage))>0)
        scsimage=permute(scsimage,[2 3 1 4]);
      end;
    case {3 4} %'+Y (lowest Y in front)' %'-Y (highest Y in front)'
      if (min(size(segimage))>0)
        segimage=permute(segimage,[1 3 2]);
      end;
      if (min(size(emimage))>0)
        if (size(size(emimage),2)==3)
          emimage=permute(emimage,[1 3 2]);
        else
          emimage=permute(emimage,[1 3 2 4]);
        end;
      end;
      if (min(size(scsimage))>0)
        scsimage=permute(scsimage,[1 3 2 4]);
      end;
    case {5 6} %'+Z (lowest Z in front)' %'-Z (highest Z in front)'
      %All fine
    end;
    
    %Go through slab slice by slice
    if (min(size(segimage))>0)
      zlist=1:1:size(segimage,3);
    else
      if (min(size(emimage))>0)
        zlist=1:1:size(emimage,3);
      else
        zlist=1:1:size(scsimage,3);
      end;
    end;
    
    if (tdir<0)
      zlist=fliplr(zlist);
    end;

    message={'Generating Projection Image ...',sprintf('Processing Slice %d...',depth)};
    set(vdata.ui.message,'String',message);
    pause(0.01);
    
    for i=zlist
      depthattenuate=1-(depth/tdepth);
      depthattenuate=depthattenuate*(1-depthattenuation)+depthattenuation;
      
      if ((imagesource==1)||(opacitysource<3))
        simg=squeeze(segimage(:,:,i))';
      end;
        
      if (expandsegmentation~=0)
        %cheap expand segments (uses max)
        if ((imagesource==1)||(opacitysource<3))
          simg2=simg;
          for j=1:1:expandsegmentation
            simg2(1:end-1,:)=max(simg2(1:end-1,:), simg2(2:end,:));
            simg2(2:end,:)=max(simg2(1:end-1,:), simg2(2:end,:));
            simg2(:,1:end-1)=max(simg2(:,1:end-1), simg2(:,2:end));
            simg2(:,2:end)=max(simg2(:,1:end-1), simg2(:,2:end));
          end;
          simg(simg==0)=simg2(simg==0);
        end;
      end;
      
      %GENERATE IMAGE SOURCE
      rsim=zeros(theight,twidth); gsim=zeros(theight,twidth); bsim=zeros(theight,twidth);
      switch imagesource
        case 1 %1: Segmentation
          rsim(simg>0)=data(simg(simg>0),3);
          gsim(simg>0)=data(simg(simg>0),4);
          bsim(simg>0)=data(simg(simg>0),5);
        case 2 %2: Selected EM layer
          if (size(size(emimage),2)==3)
            emimg=double(squeeze(emimage(:,:,i))');
            rsim=emimg;
            gsim=emimg;
            bsim=emimg;
          else
            rsim=double(squeeze(emimage(:,:,i,1))');
            gsim=double(squeeze(emimage(:,:,i,2))');
            bsim=double(squeeze(emimage(:,:,i,3))');
          end;
        case 3 %3: Screenshots
          if (size(size(scsimage),2)==4)
            rsim=double(squeeze(scsimage(:,:,i,1))');
            gsim=double(squeeze(scsimage(:,:,i,2))');
            bsim=double(squeeze(scsimage(:,:,i,3))');
          else
            rsim=double(squeeze(scsimage(:,:,1))');
            gsim=double(squeeze(scsimage(:,:,2))');
            bsim=double(squeeze(scsimage(:,:,3))');
          end;
        case 4 %4: Segment volume colors
          rsim(simg>0)=vcols(simg(simg>0),1);
          gsim(simg>0)=vcols(simg(simg>0),2);
          bsim(simg>0)=vcols(simg(simg>0),3);
      end;
      
      %GENERATE ALPHA SOURCE
      stranspmap=ones(theight,twidth);
      switch opacitysource
        case 1  %Segmented
          stranspmap(simg>0)=1-objectopacity;
        case 2  %Unsegmented
          stranspmap(simg==0)=1-objectopacity;
        case 3  %All
          stranspmap=stranspmap*(1-objectopacity);
      end;
      
      if (blurdistance>0)
        stranspmap2 = imfilter(stranspmap,edgeblurK,'same');
        stranspmap2(stranspmap2<0)=0;
        stranspmap2(stranspmap2>1)=1;
        stranspmap(stranspmap2>stranspmap)=stranspmap2(stranspmap2>stranspmap);
      end;
      
      %APPLY DEPTH ATTENUATION
      if (depthattenuate~=1)
        rsim=rsim*depthattenuate;
        gsim=gsim*depthattenuate;
        bsim=bsim*depthattenuate;
      end;
      
      %APPLY SHADOW MAP
      if (vdata.data.exportproj.useshadows)
        rsim=rsim.*shadow;
        gsim=gsim.*shadow;
        bsim=bsim.*shadow;
      end;
      
      %COMBINE WITH TARGET
      switch blendmode
        case 1 %1: Alpha-blend
          rsim=rsim.*(1-stranspmap);
          gsim=gsim.*(1-stranspmap);
          bsim=bsim.*(1-stranspmap);

          rtim=rtim + rsim.*(ttranspmap);  %transpmap is 1 for transparent and 0 for opaque
          gtim=gtim + gsim.*(ttranspmap);
          btim=btim + bsim.*(ttranspmap);
          ttranspmap=ttranspmap.*stranspmap;
        case 2 %2: Additive
          rtim=rtim + rsim.*(1-stranspmap);
          gtim=gtim + gsim.*(1-stranspmap);
          btim=btim + bsim.*(1-stranspmap);
        case 3 %3: Maximum projection
          rim=rsim.*(1-stranspmap); rtim(rim>rtim)=rim(rim>rtim);
          gim=gsim.*(1-stranspmap); gtim(gim>gtim)=gim(gim>gtim);
          bim=bsim.*(1-stranspmap); btim(bim>btim)=bim(bim>btim);
      end;
      
      shadow = imfilter(shadow,shadowK,'same');
      shadow=shadow.*stranspmap;
      
      zmap((zmap==-1)&(stranspmap<1))=depth;
      
      depth=depth+1;
    end;

    if (vdata.data.exportproj.showinwindow==1)
      targetimage=zeros(theight,twidth,3);
      targetimage(:,:,1)=rtim;
      targetimage(:,:,2)=gtim;
      targetimage(:,:,3)=btim;
      set(0, 'currentfigure', lfig);
      imshow(targetimage/255);
      title(sprintf('Rendering in progress; %d of %d ...',d,nrofslabs));
    end;
    d=d+1;
  end;
  
  targetimage=zeros(theight,twidth,3);
  targetimage(:,:,1)=rtim;
  targetimage(:,:,2)=gtim;
  targetimage(:,:,3)=btim;
  
  %%%%%%%%%%%%%%%%%
  
  if (vdata.data.exportproj.finalnormalize)
    maxval=max(targetimage(:));
    targetimage=targetimage/maxval*255;
  end;
  
  if ((vdata.data.exportproj.stretchz>1)&&(vdata.data.exportproj.projaxis~=5)&&(vdata.data.exportproj.projaxis~=6))
    mipfakt=bitshift(1,param.miplevel);
    vx=vinfo.voxelsizex*mipfakt;
    aspect=vinfo.voxelsizez/vx;
    txs=size(targetimage,2);
    tys=round(size(targetimage,1)*aspect);
    if (vdata.data.exportproj.stretchz==2)
      targetimage=imresize(targetimage,[tys txs],'nearest');
      zmap=imresize(zmap,[tys txs],'nearest');
    else
      targetimage=imresize(targetimage,[tys txs]);
      zmap=imresize(zmap,[tys txs],'nearest');
    end;
  else
    aspect=1.0;
  end;
  
  switch vdata.data.exportproj.projaxis
    case 1 %'+X (lowest X in front)'
    case 2 %'-X (highest X in front)'
      targetimage=flipdim(targetimage,2);
      zmap=flipdim(zmap,2);
    case 3 %'+Y (lowest Y in front)'
    case 4 %'-Y (highest Y in front)'
      targetimage=flipdim(targetimage,2);
      zmap=flipdim(zmap,2);
    case 5 %'+Z (lowest Z in front)'
    case 6 %'-Z (highest Z in front)'
      targetimage=flipdim(targetimage,2);
      zmap=flipdim(zmap,2);
  end;
  
  if (vdata.data.exportproj.savetofile==1)
    filename =[vdata.data.exportproj.targetfoldername '/' vdata.data.exportproj.targetfilename];
    imwrite(targetimage/255,filename);
  end;

  if (vdata.data.exportproj.showinwindow==0)
    lfig=figure;
  end;
  figure(lfig);
  imshow(targetimage/255);
  if (vdata.state.lastcancel==0)
    title('Final');
  else
    title('Canceled');
  end;
  
  %Save parameters of last render for simple navigator
  vdata.data.exportproj.lastimage.image=targetimage;
  vdata.data.exportproj.lastimage.zmap=zmap;
  vdata.data.exportproj.lastimage.stretchz=aspect;
  vdata.data.exportproj.lastimage.region=vdata.data.region;
  vdata.data.exportproj.lastimage.projectback=projectback;
  
  vdata.vast.setsegtranslation([],[]);

  if (vdata.state.lastcancel==0)
    set(vdata.ui.message,'String','Done.');
  else
    set(vdata.ui.message,'String','Canceled.');
  end;
  set(vdata.ui.cancelbutton,'Enable','off');
  vdata.state.lastcancel=0;
  pause(0.1);

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Volume Measurement Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [] = callback_measurevol(varargin)
  global vdata;

  if (vdata.state.isconnected==0)
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
    return;
  end;
  
  vinfo=vdata.vast.getinfo();
  if (min(size(vinfo)==0))
    warndlg('ERROR: Requesting data from VAST failed.','Error with remote connection to VAST');
    return;
  end;
  
  if (min([vinfo.datasizex vinfo.datasizey vinfo.datasizez])==0)
    warndlg('ERROR: No volume open in VAST.','VastTools volume measurement');
    return;
  end;
  
  nrofsegments=vdata.vast.getnumberofsegments();
  if (nrofsegments==0)
    warndlg('ERROR: No segmentation available in VAST.','VastTools volume measurement');
    return;
  end;
  
  blockgui();
  
  %Display parameter dialog
  if (~isfield(vdata.data,'region'))
    vdata.data.region.xmin=0;
    vdata.data.region.xmax=vinfo.datasizex-1;
    vdata.data.region.ymin=0;
    vdata.data.region.ymax=vinfo.datasizey-1;
    vdata.data.region.zmin=0; %first slice
    vdata.data.region.zmax=vinfo.datasizez-1; %last slice
  end;
  if (~isfield(vdata.data,'measurevol'))
    vdata.data.measurevol.miplevel=0;
    vdata.data.measurevol.xunit=vinfo.voxelsizex;
    vdata.data.measurevol.yunit=vinfo.voxelsizey;
    vdata.data.measurevol.zunit=vinfo.voxelsizez;
    vdata.data.measurevol.analyzewhich=1;
    vdata.data.measurevol.targetfoldername=[pwd '\'];
    vdata.data.measurevol.targetfilename='volumestats.txt';
  end;
  
  scrsz = get(0,'ScreenSize');
  figheight=360;
  f = figure('units','pixels','position',[50 scrsz(4)-100-figheight 500 figheight],'menubar','none','numbertitle','off','name','VastTools - Measure Segment Volumes','resize','off');

  vpos=figheight-30;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 400 15], 'Tag','t1','String',sprintf('VAST reports the voxel size to be: (%.2f nm, %.2f nm, %.2f nm)',vinfo.voxelsizex,vinfo.voxelsizey,vinfo.voxelsizez),'backgroundcolor',get(f,'color'),'horizontalalignment','left');
  set(t,'tooltipstring','To change, enter the values in VAST under "Info / Volume properties" and save to your EM stack file.');
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15], 'Tag','t1','String','Analyze at resolution:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(vinfo.nrofmiplevels,1);
  vx=vinfo.voxelsizex;
  vy=vinfo.voxelsizey;
  vz=vinfo.voxelsizez;
  for i=1:1:vinfo.nrofmiplevels
    str{i}=sprintf('Mip %d - (%.2f nm, %.2f nm, %.2f nm) voxels',i-1,vx,vy,vz);
    vx=vx*2; vy=vy*2;
  end;
  pmh = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.measurevol.miplevel+1,'Position',[170 vpos 300 20]);
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 120 15],'String','Analyze in area:','backgroundcolor',get(f,'color'),'horizontalalignment','left');

  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos 120 20], 'String','Set to full', 'CallBack',{@callback_region_settofull,0});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-25 120 20], 'String','Set to current voxel', 'CallBack',{@callback_region_settocurrentvoxel,0});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[340 vpos-50 120 20], 'String','Extend to current voxel', 'CallBack',{@callback_region_extendtocurrentvoxel,0});
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','X min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d', vdata.data.region.xmin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','X max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_xmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.xmax),'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Y min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Y max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_ymax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.ymax),'horizontalalignment','left');
  vpos=vpos-30;
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[130 vpos 100 15],'String','Z min:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmin = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[170 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmin),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[230 vpos 100 15],'String','Z max:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e_zmax = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[270 vpos 50 20],'String',sprintf('%d',vdata.data.region.zmax),'horizontalalignment','left');
  vpos=vpos-40;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 150 15],'String','Voxel size (full res)  X:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e7 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[150 vpos 50 20],'String',sprintf('%f', vdata.data.measurevol.xunit),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[210 vpos 150 15],'String','Y:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e8 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[230 vpos 50 20],'String',sprintf('%f', vdata.data.measurevol.yunit),'horizontalalignment','left');
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[300 vpos 150 15],'String','Z:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.e9 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[320 vpos 50 20],'String',sprintf('%f', vdata.data.measurevol.zunit),'horizontalalignment','left');
  vpos=vpos-40;

  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 100 15],'String','Analyze what:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  str=cell(4,1);
  str{1}='All segments individually, uncollapsed';
  str{2}='All segments, collapsed as in VAST';
  str{3}='Selected segment and children, uncollapsed';
  str{4}='Selected segment and children, collapsed as in VAST';
  pmh2 = uicontrol('Style','popupmenu','String',str,'Value',vdata.data.measurevol.analyzewhich,'Position',[120 vpos 290 20]);
  vpos=vpos-30;
  
  t = uicontrol('Style','text', 'Units','Pixels', 'Position',[30 vpos 280 15],'String','Target text file for volume measurement results:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vpos=vpos-25;
  vdata.temp.e10 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[30 vpos 380 20],'String',[vdata.data.measurevol.targetfoldername vdata.data.measurevol.targetfilename],'horizontalalignment','left');
  p1 = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[420 vpos 60 20], 'String','Browse ...','CallBack',{@callback_measurevol_browse});
  vpos=vpos-30;
  
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[150 20 60 20], 'String','OK', 'CallBack',{@callback_done});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[300 20 60 20], 'String','Cancel', 'CallBack',{@callback_canceled});

  vdata.state.lastcancel=1;
  vdata.ui.temp.closefig=0;
  uiwait(f);

  if (vdata.state.lastcancel==0)
    vdata.data.measurevol.miplevel=get(pmh,'value')-1;
    vdata.data.region.xmin = str2num(get(vdata.temp.e_xmin,'String'));
    vdata.data.region.xmax = str2num(get(vdata.temp.e_xmax,'String'));
    vdata.data.region.ymin = str2num(get(vdata.temp.e_ymin,'String'));
    vdata.data.region.ymax = str2num(get(vdata.temp.e_ymax,'String'));
    vdata.data.region.zmin = str2num(get(vdata.temp.e_zmin,'String'));
    vdata.data.region.zmax = str2num(get(vdata.temp.e_zmax,'String'));
    
    vdata.data.measurevol.xunit = str2num(get(vdata.temp.e7,'String'));
    vdata.data.measurevol.yunit = str2num(get(vdata.temp.e8,'String'));
    vdata.data.measurevol.zunit = str2num(get(vdata.temp.e9,'String'));
    vdata.data.measurevol.analyzewhich =get(pmh2,'value');
    vdata.data.measurevol.exportmodestring=get(pmh2,'string');
    vdata.data.measurevol.exportmodestring=vdata.data.measurevol.exportmodestring{vdata.data.measurevol.analyzewhich};
    vdata.data.measurevol.targetfile = get(vdata.temp.e10,'String');
    [vdata.data.measurevol.targetfoldername,vdata.data.measurevol.targetfilename]=splitfilename(vdata.data.measurevol.targetfile);
  end;
  
  if (vdata.ui.temp.closefig==1) %to distinguish close on button press and close on window x
    close(f);
  end;

  if (vdata.state.lastcancel==0)
    measurevolumes();
  end;
  releasegui();
  

function [] = callback_measurevol_browse(varargin)
  global vdata;
  defaultname=[vdata.data.measurevol.targetfoldername vdata.data.measurevol.targetfilename];
  [targetfilename,targetfoldername,filterindex] = uiputfile({'*.txt','Text Files (*.txt)'; '*.*', 'All Files (*.*)'},'Select Target Text File',defaultname);
  if ((~isequal(targetfilename,0)) && (~isequal(targetfoldername,0)))
    set(vdata.temp.e10,'String',[targetfoldername targetfilename]);
    vdata.data.measurevol.targetfilename=targetfilename;
    vdata.data.measurevol.targetfoldername=targetfoldername;
  end;

  
function [foldername,filename]=splitfilename(fullfilename)
  p=size(fullfilename,2);
  while (p>1)&&(fullfilename(p)~='\')&&(fullfilename(p)~='/')
    p=p-1;
  end;
  foldername=fullfilename(1:p);
  filename=fullfilename(p+1:end);

  
function res=measurevolumes()
  global vdata;
  
  set(vdata.ui.cancelbutton,'Enable','on');
  
  miplevel=vdata.data.measurevol.miplevel;   %bitshift(vdata.data.measurevol.xmax,-5)
  areaxmin=bitshift(vdata.data.region.xmin,-miplevel);
  areaxmax=bitshift(vdata.data.region.xmax+1,-miplevel)-1;
  areaymin=bitshift(vdata.data.region.ymin,-miplevel);
  areaymax=bitshift(vdata.data.region.ymax+1,-miplevel)-1;
  areazmin=vdata.data.region.zmin;
  areazmax=vdata.data.region.zmax;
  
  tilesizex=1024;
  tilesizey=1024;
  tilesizez=64;
  
  tilestartx=[areaxmin:tilesizex:areaxmax];
  tilestarty=[areaymin:tilesizey:areaymax];
  tilestartz=[areazmin:tilesizez:areazmax];
  nrxtiles=size(tilestartx,2);
  nrytiles=size(tilestarty,2);
  nrztiles=size(tilestartz,2);
  
  analyzewhich=vdata.data.measurevol.analyzewhich;
  
  data=vdata.vast.getallsegmentdatamatrix();
  name=vdata.vast.getallsegmentnames();
  name(1)=[];  %remove 'Background'
  
  % Compute list of objects to export
  switch analyzewhich
    case 1  %All segments individually, uncollapsed
      objects=uint32([data(:,1) data(:,2)]); 
      vdata.vast.setsegtranslation([],[]);

    case 2  %All segments, collapsed as in Vast
      %4: Collapse segments as in the view during segment text file exporting
      objects=unique(data(:,18));
      objects=uint32([objects data(objects,2)]);
      vdata.vast.setsegtranslation(data(:,1),data(:,18));
      
    case 3  %Selected segment and children, uncollapsed
      selected=find(bitand(data(:,2),65536)>0);
      if (min(size(selected))==0)
        objects=uint32([data(:,1) data(:,2)]); 
      else
        selected=[selected getchildtreeids(data,selected)];
        objects=uint32([selected' data(selected,2)]);
      end;
      vdata.vast.setsegtranslation(data(selected,1),data(selected,1));
      
    case 4  %Selected segment and children, collapsed as in Vast
      selected=find(bitand(data(:,2),65536)>0);
      if (min(size(selected))==0)
        %None selected: choose all, collapsed
        selected=data(:,1);
        objects=unique(data(:,18));
      else
        selected=[selected getchildtreeids(data,selected)];
        objects=unique(data(selected,18));
      end;

      objects=uint32([objects data(objects,2)]);
      vdata.vast.setsegtranslation(data(selected,1),data(selected,18));
  end;
  
  nrvox=zeros(max(objects(:,1)),1);
  
  z=1;
  while ((z<=size(tilestartz,2))&&(vdata.state.lastcancel==0))
    minz=tilestartz(z); maxz=min([tilestartz(z)+tilesizez-1 areazmax]);
    y=1;
    while ((y<=size(tilestarty,2))&&(vdata.state.lastcancel==0))
      miny=tilestarty(y); maxy=min([tilestarty(y)+tilesizey-1 areaymax]);
      x=1;
      while ((x<=size(tilestartx,2))&&(vdata.state.lastcancel==0))
        minx=tilestartx(x); maxx=min([tilestartx(x)+tilesizex-1 areaxmax]);
        
        message={'Measuring Volumes ...',sprintf('Loading Segmentation Cube (%d,%d,%d) of (%d,%d,%d)...',x,y,z,nrxtiles,nrytiles,nrztiles)};
        set(vdata.ui.message,'String',message);
        pause(0.01);
        
        %Volumes
        [v,n,res] = vdata.vast.getRLEcountunique(miplevel,minx,maxx,miny,maxy,minz,maxz,0);
        if (res==1)
          n(v==0)=[];
          v(v==0)=[];
          if (min(size(v))>0)
            nrvox(v)=nrvox(v)+n;
          end;
        end;
        x=x+1;
      end;
      y=y+1;
    end;
    z=z+1;
  end;
  
  nrvox=nrvox(objects(:,1));
  
  vdata.data.measurevol.lastobjects=objects;
  vdata.data.measurevol.lastvolume=nrvox;
  
  vdata.vast.setsegtranslation([],[]);
  
  if (vdata.state.lastcancel==0)
    %write surface area values to text file

    mipfact=bitshift(1,vdata.data.measurevol.miplevel);
    voxsizex=vdata.data.measurevol.xunit*mipfact;
    voxsizey=vdata.data.measurevol.yunit*mipfact;
    voxsizez=vdata.data.measurevol.zunit;
    voxelvol=voxsizex*voxsizey*voxsizez;
    fid = fopen(vdata.data.measurevol.targetfile, 'wt');
    if (fid>0)
      fprintf(fid,'%% VastTools Object Volume Export\r\n');
      fprintf(fid,'%% Provided as-is, no guarantee for correctness!\r\n');
      fprintf(fid,'%% %s\r\n\r\n',get(vdata.fh,'name'));
      
      fprintf(fid,'%% Source File: %s\r\n',getseglayername());
      fprintf(fid,'%% Mode: %s\r\n', vdata.data.measurevol.exportmodestring);
      fprintf(fid,'%% Area: (%d-%d, %d-%d, %d-%d)\r\n',areaxmin,areaxmax,areaymin,areaymax,areazmin,areazmax);
      fprintf(fid,'%% Computed at voxel size: (%f,%f,%f)\r\n',voxsizex,voxsizey,voxsizez);
      fprintf(fid,'%% Columns are: Object Name, Object ID, Voxel Count, Object Volume\r\n\r\n');

      for segnr=1:1:size(objects,1)
        seg=objects(segnr,1);
        fprintf(fid,'"%s"  %d  %d  %f\r\n',name{seg},seg,nrvox(segnr),nrvox(segnr)*voxelvol);
      end;
      fprintf(fid,'\r\n');
      fclose(fid);
    end;
  end;
  
  if (vdata.state.lastcancel==0)
    set(vdata.ui.message,'String','Done.');
  else
    set(vdata.ui.message,'String','Canceled.');
  end;
  set(vdata.ui.cancelbutton,'Enable','off');
  vdata.state.lastcancel=0;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euclidian Measurement Tool

function [] = callback_euclidiantool(varargin)
  global vdata;
  
  %Check if tool already open
  if (ishandle(vdata.etfh))
    figure(vdata.etfh);
    return;
  end;
  
  if (isfield('vdata','temp')==0)
    vdata.temp.et.c1=[0 0 0];
  end;
  if (isfield('vdata.temp','et')==0)
    vdata.temp.et.c1=[0 0 0];
    vdata.temp.et.voxelsize=[1 1 1];
    if (vdata.state.isconnected)
      [x,y,z]=vdata.vast.getviewcoordinates();
      vdata.temp.et.c1=[x y z];
      vinfo=vdata.vast.getinfo();
      vdata.temp.et.voxelsize=[vinfo.voxelsizex vinfo.voxelsizey vinfo.voxelsizez];
    end;
    vdata.temp.et.c2=[0 0 0];
    c1=double(vdata.temp.et.c1); c2=double(vdata.temp.et.c2);
    vdata.temp.et.voxdist=sqrt(sum((c1-c2).*(c1-c2)));
    c1nm=c1.*double(vdata.temp.et.voxelsize); c2nm=c2.*double(vdata.temp.et.voxelsize);
    vdata.temp.et.nmdist=sqrt(sum((c1nm-c2nm).*(c1nm-c2nm)));
  end;
  
  %blockgui();
  scrsz = get(0,'ScreenSize');
  figheight=190;
  f = figure('units','pixels','position',[50 scrsz(4)-100-figheight 440 figheight],'menubar','none','numbertitle','off','name','VastTools - Euclidian Distance Measurement','resize','off');
  vdata.etfh=f;
  vpos=figheight-40;
  
  c1=vdata.temp.et.c1; c2=vdata.temp.et.c2; vs=vdata.temp.et.voxelsize;

  uicontrol('Style','text', 'Units','Pixels', 'Position',[20 vpos 150 15],'String','First Coordinate:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.et.e1 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[180 vpos 140 20],'String',sprintf('(%d, %d, %d)',c1(1),c1(2),c1(3)),'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[330 vpos 40 20], 'String','Get', 'CallBack',{@callback_et_get1});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[380 vpos 40 20], 'String','GO!', 'CallBack',{@callback_et_go1});
  vpos=vpos-30;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[20 vpos 150 15],'String','Second Coordinate:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.et.e2 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[180 vpos 140 20],'String',sprintf('(%d, %d, %d)',c2(1),c2(2),c2(3)),'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[330 vpos 40 20], 'String','Get', 'CallBack',{@callback_et_get2});
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[380 vpos 40 20], 'String','GO!', 'CallBack',{@callback_et_go2}); 
  vpos=vpos-30;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[20 vpos 150 15],'String','Voxel Size:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.et.e3 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[180 vpos 140 20],'String',sprintf('(%.02f, %.02f, %.02f)',vs(1),vs(2),vs(3)),'horizontalalignment','left');
  p = uicontrol('Style','PushButton', 'Units','Pixels', 'Position',[330 vpos 60 20], 'String','Update', 'CallBack',{@callback_et_getvoxelsize});
  vpos=vpos-30;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[20 vpos 170 15],'String','Euclidian Distance in Voxels:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.et.e4 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[180 vpos 130 20],'String',sprintf('%f',vdata.temp.et.voxdist),'horizontalalignment','left');
  vpos=vpos-30;
  
  uicontrol('Style','text', 'Units','Pixels', 'Position',[20 vpos 170 15],'String','Euclidian Distance in nm:','backgroundcolor',get(f,'color'),'horizontalalignment','left');
  vdata.temp.et.e5 = uicontrol('Style','Edit', 'Units','Pixels', 'Position',[180 vpos 130 20],'String',sprintf('%f',vdata.temp.et.nmdist),'horizontalalignment','left');

  set(vdata.ui.menu.euclidiantool,'checked','on');
  uiwait(f);
  
  set(vdata.ui.menu.euclidiantool,'checked','off');

  
function [] = callback_et_go1(varargin)
  global vdata;
  if (vdata.state.isconnected)
    x=vdata.temp.et.c1(1);
    y=vdata.temp.et.c1(2);
    z=vdata.temp.et.c1(3);
    vdata.vast.setviewcoordinates(x,y,z);
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;  

function [] = callback_et_go2(varargin)
  global vdata;
  if (vdata.state.isconnected)
    x=vdata.temp.et.c2(1);
    y=vdata.temp.et.c2(2);
    z=vdata.temp.et.c2(3);
    vdata.vast.setviewcoordinates(x,y,z);
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;  

function [] = callback_et_get1(varargin)
  global vdata;
  
  if (vdata.state.isconnected)
    [x,y,z]=vdata.vast.getviewcoordinates();
    vdata.temp.et.c1=[x y z];
    set(vdata.temp.et.e1,'string',sprintf('(%d, %d, %d)',x,y,z));
    c1=double(vdata.temp.et.c1); c2=double(vdata.temp.et.c2);
    vdata.temp.et.voxdist=sqrt(sum((c1-c2).*(c1-c2)));
    set(vdata.temp.et.e4,'String',sprintf('%f',vdata.temp.et.voxdist));
    c1nm=c1.*double(vdata.temp.et.voxelsize); c2nm=c2.*double(vdata.temp.et.voxelsize);
    vdata.temp.et.nmdist=sqrt(sum((c1nm-c2nm).*(c1nm-c2nm)));
    set(vdata.temp.et.e5,'String',sprintf('%f',vdata.temp.et.nmdist));
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;
  
function [] = callback_et_get2(varargin)
  global vdata;
  
  if (vdata.state.isconnected)
    [x,y,z]=vdata.vast.getviewcoordinates();
    vdata.temp.et.c2=[x y z];
    set(vdata.temp.et.e2,'string',sprintf('(%d, %d, %d)',x,y,z));
    c1=double(vdata.temp.et.c1); c2=double(vdata.temp.et.c2);
    vdata.temp.et.voxdist=sqrt(sum((c1-c2).*(c1-c2)));
    set(vdata.temp.et.e4,'String',sprintf('%f',vdata.temp.et.voxdist));
    c1nm=c1.*double(vdata.temp.et.voxelsize); c2nm=c2.*double(vdata.temp.et.voxelsize);
    vdata.temp.et.nmdist=sqrt(sum((c1nm-c2nm).*(c1nm-c2nm)));
    set(vdata.temp.et.e5,'String',sprintf('%f',vdata.temp.et.nmdist));
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;
  
function [] = callback_et_getvoxelsize(varargin)
  global vdata;
  
  if (vdata.state.isconnected)
    vinfo=vdata.vast.getinfo();
    vdata.temp.et.voxelsize=[vinfo.voxelsizex vinfo.voxelsizey vinfo.voxelsizez];
    vs=vdata.temp.et.voxelsize;
    set(vdata.temp.et.e3,'String',sprintf('(%.02f, %.02f, %.02f)',vs(1),vs(2),vs(3)));
    c1=double(vdata.temp.et.c1); c2=double(vdata.temp.et.c2);
    vdata.temp.et.voxdist=sqrt(sum((c1-c2).*(c1-c2)));
    set(vdata.temp.et.e4,'String',sprintf('%f',vdata.temp.et.voxdist));
    c1nm=c1.*double(vdata.temp.et.voxelsize); c2nm=c2.*double(vdata.temp.et.voxelsize);
    vdata.temp.et.nmdist=sqrt(sum((c1nm-c2nm).*(c1nm-c2nm)));
    set(vdata.temp.et.e5,'String',sprintf('%f',vdata.temp.et.nmdist));
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target List Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = callback_newtargetlist(varargin)
  global vdata;
  
  name = inputdlg({'Enter Target List Name:'},'VastTool - New Target List',1,{'VAST Targets'});
  targetlistname=name{1};
  vdata.data.nroftargetlists=vdata.data.nroftargetlists+1;
  targetlistwindow(targetlistname,vdata.data.nroftargetlists,[]);
  
  
function [] = callback_loadtargetlist(varargin)
  global vdata;  
  
[filename, pathname] = uigetfile({'*.mat';'*.*'},'Select target list file to open...');
  if (filename==0)
    %'Cancel' was pressed. Don't load.
    return;
  end;

  vdata.data.nroftargetlists=vdata.data.nroftargetlists+1;
  targetlistwindow([],vdata.data.nroftargetlists,[pathname filename]);
  
  
function [] = targetlistwindow(targetlistname, instance, inputfilename)
  global vdata;
  
  scrsz = get(0,'ScreenSize');
  vdata.data.tl(instance).fh = figure('units','pixels','outerposition',[300 scrsz(4)-639-300 1024 640],...
    'menubar','none','numbertitle','off','resize','on','name',[targetlistname ' [VastTools Target List]']);
  set(vdata.data.tl(instance).fh,'CloseRequestFcn',{@callback_tlquit, instance});
  set(vdata.data.tl(instance).fh,'ResizeFcn',{@callback_tlresize, instance});
  vdata.data.tl(instance).open=1;
  
  vdata.data.tl(instance).menu.file = uimenu(vdata.data.tl(instance).fh,'Label','File');
  vdata.data.tl(instance).menu.savetargetlist = uimenu(vdata.data.tl(instance).menu.file,'Label','Save Target List ...','Callback',{@callback_savetargetlist, instance});
  vdata.data.tl(instance).menu.close = uimenu(vdata.data.tl(instance).menu.file,'Label','Close Target List','Callback',{@callback_tlquit, instance});
  
  vdata.data.tl(instance).menu.edit = uimenu(vdata.data.tl(instance).fh,'Label','Edit');
  vdata.data.tl(instance).menu.cutselectedrows = uimenu(vdata.data.tl(instance).menu.edit,'Label','Cut Selected Rows','Callback',{@callback_tlcutselectedrows, instance});
  vdata.data.tl(instance).menu.copyselectedrows = uimenu(vdata.data.tl(instance).menu.edit,'Label','Copy Selected Rows','Callback',{@callback_tlcopyselectedrows, instance});
  vdata.data.tl(instance).menu.pasteselectedrows = uimenu(vdata.data.tl(instance).menu.edit,'Label','Paste Rows Below Selected','Callback',{@callback_tlpasteselectedrows, instance});
  vdata.data.tl(instance).menu.insertseparator = uimenu(vdata.data.tl(instance).menu.edit,'Label','Insert Separator Below Selected Row','Separator','on','Callback',{@callback_tlinsertseparator, instance});
  
  columnname =   {'Click', '       Coordinates       ', 'Zoom', 'Sel. Segnr', '   Target Name   ', '     Properties     ', '       Notes       ', '       Comments       '};
  columnformat = {'char', 'char', 'char', 'char', 'char', 'char', 'char', 'char'};
  columneditable =  [false false false false true true true true];
  
  vdata.data.tl(instance).sendzoom=1;
  vdata.data.tl(instance).sendselected=1;

  pos=get(vdata.data.tl(instance).fh,'Position');
  
  vdata.data.tl(instance).ui.addbutton = uicontrol('style','push','units','pixels','position',[150 pos(4)-65 250 30],...
    'fontsize',12,'string','Add Current VAST Location','callback',{@callback_tladdcurrentlocation, instance});
  
  vdata.data.tl(instance).ui.sendzoom = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[10 pos(4)-45 100 15],'Value',vdata.data.tl(instance).sendzoom,'string','Send Zoom','backgroundcolor',get(vdata.data.tl(instance).fh,'color'),'callback',{@callback_tlsendzoom, instance});
  vdata.data.tl(instance).ui.sendselect = uicontrol('Style','checkbox', 'Units','Pixels', 'Position',[10 pos(4)-65 100 15],'Value',vdata.data.tl(instance).sendselected,'string','Send Selected','backgroundcolor',get(vdata.data.tl(instance).fh,'color'),'callback',{@callback_tlsendselected, instance});
  
  vdata.data.tl(instance).ui.table = uitable('Units','pixels','Position', [10 10 pos(3)-20 pos(4)-80],'ColumnName', columnname, ...
    'ColumnFormat', columnformat, 'ColumnEditable', columneditable, 'RowName',[], ...
    'CellSelectionCallback',{@callback_tlcellselection, instance}, 'CellEditCallback',{@callback_tlcelledit, instance}, 'ButtonDownFcn',{@callback_tltablerightclick, instance});

  vdata.data.tl(instance).selected=[];
  vdata.data.tl(instance).filename=inputfilename;
  vdata.data.tl(instance).targetlistname=targetlistname;
  
  if (min(size(inputfilename>0)))
    %Load data
    targetlist=load(inputfilename,'-mat','targetlist');
    if (~isfield(targetlist,'targetlist'))
      warndlg('This is not a valid target list file. A .MAT file with struct "targetlist" is expected.',['Error loading ' inputfilename]);
      return;
    end;
    targetlist=targetlist.targetlist;
    
    if (~isfield(targetlist,'coords'))
      warndlg('"targetlist.coords" missing in target list file.',['Error loading ' inputfilename]);
      callback_tlquit();
      return;
    end;
    
    targetlist.coords=double(targetlist.coords);
    nroftargets=size(targetlist.coords,1);

    if (~isfield(targetlist,'name'))
      for i=1:1:nroftargets
        targetlist.name{i}=sprintf('Target %d',i);
      end;
    end;
    if (~isfield(targetlist,'segmentnr'))
      for i=1:1:nroftargets
        targetlist.segmentnr(i)=-1;
      end;
    end;
    if (~isfield(targetlist,'properties'))
      for i=1:1:nroftargets
        targetlist.properties{i}='';
      end;
    end;
    if (~isfield(targetlist,'notes'))
      for i=1:1:nroftargets
        targetlist.notes{i}='';
      end;
    end;
    if (~isfield(targetlist,'comments'))
      for i=1:1:nroftargets
        targetlist.comments{i}='';
      end;
    end;

    vdata.data.tl(instance).nroftargets=nroftargets;
    vdata.data.tl(instance).coords=targetlist.coords;
    vdata.data.tl(instance).name=targetlist.name;
    vdata.data.tl(instance).segmentnr=targetlist.segmentnr;
    vdata.data.tl(instance).properties=targetlist.properties;
    vdata.data.tl(instance).notes=targetlist.notes;
    vdata.data.tl(instance).comments=targetlist.comments;
    
    set(vdata.data.tl(instance).fh,'name',[inputfilename '  [VastTools Target List]']);
    
    sz=size(vdata.data.tl(instance).segmentnr);
    if (sz(2)>sz(1)) 
      vdata.data.tl(instance).segmentnr = vdata.data.tl(instance).segmentnr'; 
    end;
    
  else
    vdata.data.tl(instance).nroftargets=0;
    vdata.data.tl(instance).coords=[];
    vdata.data.tl(instance).name={};
    vdata.data.tl(instance).segmentnr=[];
    vdata.data.tl(instance).properties={};
    vdata.data.tl(instance).notes={};
    vdata.data.tl(instance).comments={};
  end;
  vdata.data.tl(instance).ischanged=0;
  tl_updatetable(instance);
  
  
function tl_updatetable(instance)
  global vdata;

  separatorstring='<html><body bgcolor="#C0C0A6"><b>........................................</b></body></html>';

  nroftargets=vdata.data.tl(instance).nroftargets;
  if (nroftargets==0)
    dat=[];
    set(vdata.data.tl(instance).ui.table,'Data', dat);
    vdata.data.tl(instance).ui.tabledata=dat;
    return;
  end;
  
  dat=cell(nroftargets,3);
  targetlist.coords=vdata.data.tl(instance).coords;
  targetlist.name=vdata.data.tl(instance).name;
  targetlist.segmentnr=vdata.data.tl(instance).segmentnr;
  targetlist.properties=vdata.data.tl(instance).properties;
  targetlist.notes=vdata.data.tl(instance).notes;
  targetlist.comments=vdata.data.tl(instance).comments;
  
  for i=1:1:nroftargets
    if (~isnan(targetlist.coords(i,1)))
      dat{i,1}='GO!';
      dat{i,2}=sprintf('(%d, %d, %d)',targetlist.coords(i,1),targetlist.coords(i,2),targetlist.coords(i,3));
      dat{i,3}=sprintf('%d',targetlist.coords(i,4));
      dat{i,4}=targetlist.segmentnr(i);
      dat{i,5}=targetlist.name{i};
      dat{i,6}=targetlist.properties{i};
      dat{i,7}=targetlist.notes{i};
      dat{i,8}=targetlist.comments{i};
    else
      dat{i,1}= separatorstring;
      dat{i,2}= separatorstring;
      dat{i,3}= separatorstring;
      dat{i,4}= separatorstring;
      dat{i,5}=targetlist.name{i};
      dat{i,6}=targetlist.properties{i};
      dat{i,7}=targetlist.notes{i};
      dat{i,8}=targetlist.comments{i};
    end;
  end;
  
  set(vdata.data.tl(instance).ui.table,'Data', dat);
  vdata.data.tl(instance).ui.tabledata=dat;
  
  
function [] = callback_tlquit(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.tl(instance).ischanged==1)
    res = questdlg('This target list was changed. Close without saving?','Close Target List','Yes','No','Yes');
    if strcmp(res,'No') 
      return; 
    end
  end;
  
  %%%% CLEANUP
  if ishandle(vdata.data.tl(instance).fh) 
    delete(vdata.data.tl(instance).fh); 
  end
  vdata.data.tl(instance).open=0;
  vdata.data.tl(instance).fh=[];
  
  
function [] = callback_tlresize(varargin)
  global vdata;
  instance=varargin{3};
  
  set(vdata.data.tl(instance).fh,'Units','pixels');
  pos = get(vdata.data.tl(instance).fh,'OuterPosition');
  hpos=pos(3)+(-1024+560);
  vpos=pos(4)-100;
  pos=get(vdata.data.tl(instance).fh,'Position');
  
  set(vdata.data.tl(instance).ui.addbutton,'position',[150 pos(4)-40 250 30]);
  set(vdata.data.tl(instance).ui.sendzoom,'Position',[10 pos(4)-25 100 15]);
  set(vdata.data.tl(instance).ui.sendselect,'Position',[10 pos(4)-40 100 15]);
  set(vdata.data.tl(instance).ui.table,'Position', [10 10 pos(3)-20 pos(4)-60]);
  
  
function [] = callback_savetargetlist(varargin)
  global vdata;
  instance=varargin{3};
  
  if (min(size(vdata.data.tl(instance).filename>0)))
    targetname=[vdata.data.tl(instance).filename];
  else
    targetname=[vdata.data.tl(instance).targetlistname '.mat'];
  end;
  [filename, pathname] = uiputfile({'*.mat';'*.*'},'Select target list file to save...',targetname);
  if (filename==0)
    %'Cancel' was pressed. Don't save.
    return;
  end;
  
  targetlist.coords=vdata.data.tl(instance).coords;
  targetlist.name=vdata.data.tl(instance).name;
  targetlist.segmentnr=vdata.data.tl(instance).segmentnr;
  targetlist.properties=vdata.data.tl(instance).properties;
  targetlist.notes=vdata.data.tl(instance).notes;
  targetlist.comments=vdata.data.tl(instance).comments;

  save([pathname filename],'targetlist');
  vdata.data.tl(instance).filename=[pathname filename];
  set(vdata.data.tl(instance).fh,'name',[vdata.data.tl(instance).filename '  [VastTools Target List]']);
  vdata.data.tl(instance).ischanged=0;
  
function [] = callback_tlsendzoom(varargin)
  global vdata;
  instance=varargin{3};   
  
  if (get(vdata.data.tl(instance).ui.sendzoom,'Value') == get(vdata.data.tl(instance).ui.sendzoom,'Max'))
	  vdata.data.tl(instance).sendzoom=1;
  else
    vdata.data.tl(instance).sendzoom=0;
  end
  
    
function [] = callback_tlsendselected(varargin)
  global vdata;
  instance=varargin{3};
  
  if (get(vdata.data.tl(instance).ui.sendselect,'Value') == get(vdata.data.tl(instance).ui.sendselect,'Max'))
	  vdata.data.tl(instance).sendselected=1;
  else
    vdata.data.tl(instance).sendselected=0;
  end
  

function [] = callback_tlcelledit(varargin)
  global vdata;
  instance=varargin{3};
  
  row=varargin{2}.Indices(1);
  col=varargin{2}.Indices(2);
  if (col==5) %target name
    vdata.data.tl(instance).name{row}=varargin{2}.NewData;
  end;
  if (col==6) %notes
    vdata.data.tl(instance).properties{row}=varargin{2}.NewData;
  end;
  if (col==7) %notes
    vdata.data.tl(instance).notes{row}=varargin{2}.NewData;
  end;
  if (col==8) %comments
    vdata.data.tl(instance).comments{row}=varargin{2}.NewData;
  end;
  vdata.data.tl(instance).ischanged=1;
  
  
function [] = callback_tlcellselection(varargin)
  %Callback for click into target list
  global vdata;
  instance=varargin{3};
  
  selected = varargin{2}.Indices;
  
  vdata.data.tl(instance).selected=selected;
  if (min(size(selected))>0)
    if (selected(1,2)==1)
      tcoords=vdata.data.tl(instance).coords(selected(1,1),:);
      if (~isnan(tcoords(1)))
        if (vdata.state.isconnected)
          vdata.vast.setviewcoordinates(tcoords(1),tcoords(2),tcoords(3));
          if (vdata.data.tl(instance).sendzoom==1)
            vdata.vast.setviewzoom(tcoords(4));
          end;
          if (vdata.data.tl(instance).sendselected==1)
            vdata.vast.setselectedsegmentnr(vdata.data.tl(instance).segmentnr(selected(1,1)));
          end;
        else
          warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
        end;
      end;
    end;
  end;


function [] = callback_tladdcurrentlocation(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.state.isconnected)
    [tcoords(1),tcoords(2),tcoords(3)]=vdata.vast.getviewcoordinates();
    zoom=vdata.vast.getviewzoom();
    selectedsegmentnr=vdata.vast.getselectedsegmentnr();
    
    ins=size(vdata.data.tl(instance).coords,1);
    selected=vdata.data.tl(instance).selected;
    if (min(size(selected))>0)
      rowlist=unique(selected(:,1));
      ins=rowlist(end);
    end;
  
    vdata.data.tl(instance).nroftargets = vdata.data.tl(instance).nroftargets+1;
    vdata.data.tl(instance).coords = [vdata.data.tl(instance).coords(1:ins,:); double([tcoords zoom]); vdata.data.tl(instance).coords(ins+1:end,:)];
    vdata.data.tl(instance).name = [vdata.data.tl(instance).name(1:ins) sprintf('Target %d',vdata.data.tl(instance).nroftargets) vdata.data.tl(instance).name(ins+1:end)];
    vdata.data.tl(instance).segmentnr = [vdata.data.tl(instance).segmentnr(1:ins); selectedsegmentnr; vdata.data.tl(instance).segmentnr(ins+1:end)];
    vdata.data.tl(instance).properties = [vdata.data.tl(instance).properties(1:ins), ' ', vdata.data.tl(instance).properties(ins+1:end)];
    vdata.data.tl(instance).notes = [vdata.data.tl(instance).notes(1:ins), ' ', vdata.data.tl(instance).notes(ins+1:end)];
    vdata.data.tl(instance).comments = [vdata.data.tl(instance).comments(1:ins), ' ', vdata.data.tl(instance).comments(ins+1:end)];
    vdata.data.tl(instance).ischanged=1;
    tl_updatetable(instance);    
    
  else
    warndlg('ERROR: Not connected to VAST. please connect before using this function.','Not connected to VAST');
  end;

  
function [] = callback_tlinsertseparator(varargin)
  global vdata;
  instance=varargin{3};
  
  selected=vdata.data.tl(instance).selected;
  separatorstring='<html><body bgcolor="#C0C0A6"><b>........................................</b></body></html>';
  
  if (min(size(selected))>0)
    rowlist=unique(selected(:,1));
    ins=rowlist(end);
    vdata.data.tl(instance).nroftargets = vdata.data.tl(instance).nroftargets+1;
    vdata.data.tl(instance).coords = [vdata.data.tl(instance).coords(1:ins,:); [nan nan nan nan]; vdata.data.tl(instance).coords(ins+1:end,:)];
    vdata.data.tl(instance).name = [vdata.data.tl(instance).name(1:ins) separatorstring vdata.data.tl(instance).name(ins+1:end)];
    vdata.data.tl(instance).segmentnr = [vdata.data.tl(instance).segmentnr(1:ins); 0; vdata.data.tl(instance).segmentnr(ins+1:end)];
    vdata.data.tl(instance).properties = [vdata.data.tl(instance).properties(1:ins), separatorstring, vdata.data.tl(instance).properties(ins+1:end)];
    vdata.data.tl(instance).notes = [vdata.data.tl(instance).notes(1:ins), separatorstring, vdata.data.tl(instance).notes(ins+1:end)];
    vdata.data.tl(instance).comments = [vdata.data.tl(instance).comments(1:ins), separatorstring, vdata.data.tl(instance).comments(ins+1:end)];
  else
    %nothing selected: append
    vdata.data.tl(instance).nroftargets = vdata.data.tl(instance).nroftargets+1;
    vdata.data.tl(instance).coords = [vdata.data.tl(instance).coords; [nan nan nan nan]];
    vdata.data.tl(instance).name = [vdata.data.tl(instance).name separatorstring];
    vdata.data.tl(instance).segmentnr = [vdata.data.tl(instance).segmentnr; 0];
    vdata.data.tl(instance).properties = [vdata.data.tl(instance).properties, separatorstring];
    vdata.data.tl(instance).notes = [vdata.data.tl(instance).notes, separatorstring];
    vdata.data.tl(instance).comments = [vdata.data.tl(instance).comments, separatorstring];
  end;
  vdata.data.tl(instance).ischanged=1;
  tl_updatetable(instance);
  
  
function [] = callback_tlcutselectedrows(varargin)
  global vdata;
  instance=varargin{3};
  
  selected=vdata.data.tl(instance).selected;

  if (min(size(selected))>0)
    cutlist=unique(selected(:,1));
    
    vdata.data.copy.nroftargets=size(cutlist,1);
    vdata.data.copy.coords=vdata.data.tl(instance).coords(cutlist,:);
    vdata.data.copy.name=vdata.data.tl(instance).name(cutlist);
    vdata.data.copy.segmentnr=vdata.data.tl(instance).segmentnr(cutlist);
    vdata.data.copy.properties=vdata.data.tl(instance).properties(cutlist);
    vdata.data.copy.notes=vdata.data.tl(instance).notes(cutlist);
    vdata.data.copy.comments=vdata.data.tl(instance).comments(cutlist);
    
    vdata.data.tl(instance).nroftargets = vdata.data.tl(instance).nroftargets-size(cutlist,1);
    vdata.data.tl(instance).coords(cutlist,:)=[];
    vdata.data.tl(instance).name(cutlist)=[];
    vdata.data.tl(instance).segmentnr(cutlist)=[];
    vdata.data.tl(instance).properties(cutlist)=[];
    vdata.data.tl(instance).notes(cutlist)=[];
    vdata.data.tl(instance).comments(cutlist)=[];
    
    vdata.data.tl(instance).ischanged=1;
    tl_updatetable(instance);
  end;

  
function [] = callback_tlcopyselectedrows(varargin)
  global vdata;
  instance=varargin{3};
  
  selected=vdata.data.tl(instance).selected;
  if (min(size(selected))>0)
    copylist=unique(selected(:,1));
    
    vdata.data.copy.nroftargets=size(copylist,1);
    vdata.data.copy.coords=vdata.data.tl(instance).coords(copylist,:);
    vdata.data.copy.name=vdata.data.tl(instance).name(copylist);
    vdata.data.copy.segmentnr=vdata.data.tl(instance).segmentnr(copylist);
    vdata.data.copy.properties=vdata.data.tl(instance).properties(copylist);
    vdata.data.copy.notes=vdata.data.tl(instance).notes(copylist);
    vdata.data.copy.comments=vdata.data.tl(instance).comments(copylist);
  end;
  
  
function [] = callback_tlpasteselectedrows(varargin)
  global vdata;
  instance=varargin{3};
  
  if (~isfield(vdata.data,'copy'))
    return; %copy buffer is empty
  end;
  
  selected=vdata.data.tl(instance).selected;
  ins=0;
  if (min(size(selected))>0)
    ins=selected(end,1);
  end;
  
  vdata.data.tl(instance).nroftargets = vdata.data.tl(instance).nroftargets+vdata.data.copy.nroftargets;
  vdata.data.tl(instance).coords = [vdata.data.tl(instance).coords(1:ins,:); vdata.data.copy.coords; vdata.data.tl(instance).coords(ins+1:end,:)];
  vdata.data.tl(instance).name = [vdata.data.tl(instance).name(1:ins) vdata.data.copy.name vdata.data.tl(instance).name(ins+1:end)];
  vdata.data.tl(instance).segmentnr = [vdata.data.tl(instance).segmentnr(1:ins); vdata.data.copy.segmentnr; vdata.data.tl(instance).segmentnr(ins+1:end)];
  vdata.data.tl(instance).properties = [vdata.data.tl(instance).properties(1:ins), vdata.data.copy.properties, vdata.data.tl(instance).properties(ins+1:end)];
  vdata.data.tl(instance).notes = [vdata.data.tl(instance).notes(1:ins), vdata.data.copy.notes, vdata.data.tl(instance).notes(ins+1:end)];
  vdata.data.tl(instance).comments = [vdata.data.tl(instance).comments(1:ins), vdata.data.copy.comments, vdata.data.tl(instance).comments(ins+1:end)];
  
  vdata.data.tl(instance).ischanged=1;
  tl_updatetable(instance);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple Navigator Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [] = callback_newsimplenavigator(varargin)
  global vdata;
  
  if (~isfield(vdata.data,'exportproj'))
    warndlg('Please render a projection image first, using Export / Export Projection Image.','Error generating Simple Navigator Window');
    return;
  end;
  
  if (~isfield(vdata.data.exportproj,'lastimage'))
    warndlg('Please render a projection image first, using Export / Export Projection Image.','Error generating Simple Navigator Window');
    return;
  end;
  
  vdata.data.nrofsimplenavigators=vdata.data.nrofsimplenavigators+1;
  simplenavigatorwindow(vdata.data.nrofsimplenavigators,[]);
  
  
 function [] = callback_loadsimplenavigator(varargin)
  global vdata;
  
  [filename, pathname] = uigetfile({'*.mat';'*.*'},'Select simple navigator image file to open...');
  if (filename==0)
    %'Cancel' was pressed. Don't load.
    return;
  end;
  
  vdata.data.nrofsimplenavigators=vdata.data.nrofsimplenavigators+1;
  simplenavigatorwindow(vdata.data.nrofsimplenavigators,[pathname filename]);  
  
  
function [] = simplenavigatorwindow(instance, inputfilename)
  global vdata;

  %Check if the simple navigator image is loaded from a file and the file contains the correct data

  if (min(size(inputfilename))>0)
    vars = whos('-file',inputfilename);
    if (ismember('sndata', {vars.name})==0)
      warndlg(['ERROR: The file "' inputfilename '" is not a VastTools Simple Navigator Image file.'],'Error loading Simple Navigator Image');
      return;
    end;
  end;
  
  scrsz = get(0,'ScreenSize');
  vdata.data.sn(instance).fh = figure('units','pixels','outerposition',[300 scrsz(4)-639-300 640 640],...
    'menubar','none','numbertitle','off','resize','on','name','Simple Navigator Window');
  set(vdata.data.sn(instance).fh,'CloseRequestFcn',{@callback_snquit, instance});
  set(vdata.data.sn(instance).fh,'ResizeFcn',{@callback_snresize, instance});
  vdata.data.sn(instance).open=1;
  
  vdata.data.sn(instance).menu.file = uimenu(vdata.data.sn(instance).fh,'Label','File');
  vdata.data.sn(instance).menu.savetargetlist = uimenu(vdata.data.sn(instance).menu.file,'Label','Save Simple Navigator Image ...','Callback',{@callback_snsave, instance});
  vdata.data.sn(instance).menu.close = uimenu(vdata.data.sn(instance).menu.file,'Label','Close Simple Navigator Window','Callback',{@callback_snquit, instance});
  
  if (min(size(inputfilename))>0)
    load(inputfilename,'sndata');
    vdata.data.sn(instance).data=sndata;
    clear sndata;
    set(vdata.data.sn(instance).fh,'name',['Simple Navigator Image ' inputfilename]);
  else
    vdata.data.sn(instance).data=vdata.data.exportproj.lastimage;
    set(vdata.data.sn(instance).fh,'name',['New Simple Navigator Image']);
  end;
  
  vdata.data.sn(instance).ui.ax = axes('units','pixels', 'position',[20 640-512-100+28 512 512]);
  vdata.data.sn(instance).ui.XLM = get(vdata.data.sn(instance).ui.ax,'xlim');
  vdata.data.sn(instance).ui.YLM = get(vdata.data.sn(instance).ui.ax,'ylim');
  vdata.data.sn(instance).ui.AXP = get(vdata.data.sn(instance).ui.ax,'pos');
  vdata.data.sn(instance).ui.DFX = diff(vdata.data.sn(instance).ui.XLM);
  vdata.data.sn(instance).ui.DFY = diff(vdata.data.sn(instance).ui.YLM);
  
  vdata.data.sn(instance).ui.R = image(vdata.data.sn(instance).data.image/255,'Parent',vdata.data.sn(instance).ui.ax);
  set(vdata.data.sn(instance).ui.ax,'xtick',[],'ytick',[])  % Get rid of ticks.
  
  set(vdata.data.sn(instance).fh,'windowbuttonmotionfcn',{@callback_simplenavmousemove, instance}); % Set the motion detector.
  set(vdata.data.sn(instance).ui.R,'buttondownfcn',{@callback_simplenavimageclick, instance}); %Set the button press callback
  
  %%%% Toolbar
  icons=imread('vttoolicons.png');
  ht = uitoolbar(vdata.data.sn(instance).fh);
  icon=icons(:,1:19,:);
  vdata.data.sn(instance).ui.toolbar_pointer = uitoggletool(ht,'CData',icon,'TooltipString','Goto','Separator','on','OnCallback',{@callback_sntoolbar_pointer_on, instance},'OffCallback',{@callback_sntoolbar_pointer_off, instance},'state','on');
  icon=icons(:,21:39,:);
  vdata.data.sn(instance).ui.toolbar_zoom = uitoggletool(ht,'CData',icon,'TooltipString','Zoom','OnCallback',{@callback_sntoolbar_zoom_on, instance},'OffCallback',{@callback_sntoolbar_zoom_off, instance});
  icon=icons(:,41:59,:);
  vdata.data.sn(instance).ui.toolbar_pan = uitoggletool(ht,'CData',icon,'TooltipString','Pan','OnCallback',{@callback_sntoolbar_pan_on, instance},'OffCallback',{@callback_sntoolbar_pan_off, instance});
  vdata.data.sn(instance).state.actionmode=0;
  vdata.data.sn(instance).filename=inputfilename;


function callback_sntoolbar_pointer_on(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==1)
    vdata.data.sn(instance).state.actionmode=0;
    set(vdata.data.sn(instance).ui.toolbar_zoom,'State','off');
  end;
  if (vdata.data.sn(instance).state.actionmode==2)
    vdata.data.sn(instance).state.actionmode=0;
    set(vdata.data.sn(instance).ui.toolbar_pan,'State','off');
  end;
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  
function callback_sntoolbar_pointer_off(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==0) %only if switching off pointer mode
    vdata.data.sn(instance).state.actionmode=0;
    set(vdata.data.sn(instance).ui.toolbar_pointer,'State','on'); %switch back on immediately
  end;
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  
function callback_sntoolbar_zoom_on(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==0)
    vdata.data.sn(instance).state.actionmode=1;
    set(vdata.data.sn(instance).ui.toolbar_pointer,'State','off');
  end;
  if (vdata.data.sn(instance).state.actionmode==2)
    vdata.data.sn(instance).state.actionmode=1;
    set(vdata.data.sn(instance).ui.toolbar_pan,'State','off');
  end;
  
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  zoom on;
  h=zoom(vdata.data.sn(instance).fh);
  %set(h,'Motion','horizontal');
  %linkaxes(mndata.ui.ax,'xy'); %does not work?
  set(h,'ActionPostCallback', {@callback_zoom_end, instance});
  
function callback_sntoolbar_zoom_off(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==1) %only if switching off zoom mode
    vdata.data.sn(instance).state.actionmode=0;
    set(vdata.data.sn(instance).ui.toolbar_pointer,'State','on');
  end;
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  zoom off;
  
function callback_zoom_end(varargin)
  % OBJ         handle to the figure that has been clicked on.
  % EVENT_OBJ   handle to event object. The object has the same properties as the EVENT_OBJ of the 'ModePreCallback' callback.
  global vdata;
  instance=varargin{3};
  obj=varargin{1};
  event_obj=varargin{2};
  
  %Adjust zoom region to square for correct aspect ratio
  xlimits=xlim(vdata.data.sn(instance).ui.ax);
  ylimits=ylim(vdata.data.sn(instance).ui.ax);
  if (xlimits(2)-xlimits(1))>(ylimits(2)-ylimits(1))
    xsh=(xlimits(2)-xlimits(1))/2; ym=(ylimits(1)+ylimits(2))/2;
    ylim([ym-xsh,ym+xsh]);
  else
    ysh=(ylimits(2)-ylimits(1))/2; xm=(xlimits(1)+xlimits(2))/2;
    xlim([xm-ysh,xm+ysh]);
  end;
  set(event_obj,'ActionPostCallback', []);
  
function callback_sntoolbar_pan_on(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==0)
    vdata.data.sn(instance).state.actionmode=2;
    set(vdata.data.sn(instance).ui.toolbar_pointer,'State','off');
  end;
  if (vdata.data.sn(instance).state.actionmode==1)
    vdata.data.sn(instance).state.actionmode=2;
    set(vdata.data.sn(instance).ui.toolbar_zoom,'State','off');
  end;
  
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  pan on;
  
function callback_sntoolbar_pan_off(varargin)
  global vdata;
  instance=varargin{3};
  
  if (vdata.data.sn(instance).state.actionmode==2) %only if switching off pan mode
    vdata.data.sn(instance).state.actionmode=0;
    set(vdata.data.sn(instance).ui.toolbar_pointer,'State','on');
  end;
  set(0, 'currentfigure', vdata.data.sn(instance).fh);
  pan off;
  
  
function [] = callback_snsave(varargin)
  %Save a simple navigator image to a file
  global vdata;
  instance=varargin{3};
  
  if (min(size(vdata.data.sn(instance).filename>0)))
    targetname=[vdata.data.sn(instance).filename];
  else
    targetname=['simplenavigatorimage.mat'];
  end;
  [filename, pathname] = uiputfile({'*.mat';'*.*'},'Select simple navigator image file to save...',targetname);
  if (filename==0)
    %'Cancel' was pressed. Don't save.
    return;
  end;
  
  sndata=vdata.data.sn(instance).data;

  save([pathname filename],'sndata');
  vdata.data.sn(instance).filename=[pathname filename];
  set(vdata.data.sn(instance).fh,'name',['Simple Navigator Image ' vdata.data.sn(instance).filename]);
  
  
function [] = callback_snquit(varargin)
  global vdata;
  instance=varargin{3};
  
  res = questdlg('Close this simple navigator window?','Close Simple Navigator Window','Yes','No','Yes');
  if strcmp(res,'No') 
    return; 
  end
  
  %%%% CLEANUP
  if ishandle(vdata.data.sn(instance).fh) 
    delete(vdata.data.sn(instance).fh); 
  end
  vdata.data.sn(instance).open=0;
  vdata.data.sn(instance).fh=[];
  
  
function [] = callback_snresize(varargin)
  global vdata;
  instance=varargin{3};
  
  set(vdata.data.sn(instance).fh,'Units','pixels');
  pos = get(vdata.data.sn(instance).fh,'OuterPosition');
  set(vdata.data.sn(instance).ui.ax,'position',[20 20 pos(3)-50 pos(4)-100]);


function [] = callback_simplenavmousemove(varargin)
  global vdata;
  instance=varargin{3};
  
  p = get(vdata.data.sn(instance).ui.ax, 'currentpoint');
  if ((p(1)>=0) && (p(1)<size(vdata.data.sn(instance).data.image,1)) && (p(3)>=0) && (p(3)<size(vdata.data.sn(instance).data.image,2)))
    %Compute XYZ coordinates from image coordinates
  end;
  
function [] = callback_simplenavimageclick(varargin)
  global vdata;
  instance=varargin{3};
  
  p = get(vdata.data.sn(instance).ui.ax, 'currentpoint');
  if ((p(1)>=0) && (p(1)<size(vdata.data.sn(instance).data.image,2)) && (p(3)>=0) && (p(3)<size(vdata.data.sn(instance).data.image,1)))
    %Compute image coordinates
    sx=round(p(1,1));
    sy=round(p(1,2));
    sz=vdata.data.sn(instance).data.zmap(sy,sx);
    %Compute VAST coordinates from image coordinates
    sv=[sx sy sz 1]';
    tv=vdata.data.sn(instance).data.projectback*sv;
    vdata.vast.setviewcoordinates(tv(1),tv(2),round(tv(3)/vdata.data.sn(instance).data.stretchz));
  end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
function ret=getchildtreeids(data,parentlist)
  %Uses the 24-column-data matrix as analyzed from VAST color files
  %Gets a list of the IDs of the segment's children's tree (if it exists)

  ret=[];
  pal=parentlist(:);
  for p=1:1:size(pal,1)
    index=parentlist(p);
    if (data(index,15)>0)
      i=data(index,15);
      ret=[ret i getchildtreeids(data,i)]; %Add size of child tree
      while (data(i,17)>0) %add sizes of all nexts
        i=data(i,17);
        ret=[ret i getchildtreeids(data,i)];
      end;
    end;
  end;