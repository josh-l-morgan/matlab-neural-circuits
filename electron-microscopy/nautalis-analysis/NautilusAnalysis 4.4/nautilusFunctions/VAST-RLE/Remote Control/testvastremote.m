addpath './VASTControl';
vast=VASTControlClass();  %CAUTION: This contains javaaddpath which clears all globals. WTF, Matlab !

%Try to connect
res=vast.connect('127.0.0.1',22081,1000);
%res=vast.Connect('192.168.1.117',22081,1000);
if (res==0)
  warndlg('ERROR: Connecting to VAST Lite at 127.0.0.1 port 22081 failed.','Error connecting to VAST Lite');
  return;
end;
isconnected=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% info=vast.getinfo();
% 
% xs=info.datasizex; ys=info.datasizey; zs=info.datasizez;
% verts = [0 0 0;0 ys 0;xs ys 0;xs 0 0;0 0 zs;0 ys zs;xs ys zs;xs 0 zs];
% faces = [1 2 3 4;5 6 7 8;3 4 8 7;1 2 6 5;2 3 7 6;4 1 5 8];
% 
% fh=figure(1);
% h = patch('Faces',faces,'Vertices',verts,'FaceColor','b','EdgeColor','w','FaceAlpha',0.25);
% grid on;
% daspect([1/double(info.voxelsizex) 1/double(info.voxelsizey) 1/double(info.voxelsizez)]);
% 
% quit=0;
% lastvx=0;
% lastvy=0;
% lastvz=0;
% 
% while (quit==0)
%  info=vast.getinfo();
%   if ((info.currentviewx~=lastvx)||(info.currentviewy~=lastvy)||(info.currentviewz~=lastvz))
%     set(0, 'CurrentFigure', fh);
%     hold on;
%     plot3(info.currentviewx,info.currentviewy,info.currentviewz,'.');
%     hold off;
%     lastvx=info.currentviewx;
%     lastvy=info.currentviewy;
%     lastvz=info.currentviewz;
%     pause(0.1);
%   end;
% end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nrofsegments=vast.getnumberofsegments();
% 
% data=vast.getsegmentdata(16);
% name=vast.getsegmentname(16);
% 
% res=vast.setanchorpoint(16,1000,1000,1000);
% data2=vast.getsegmentdata(16);
% 
% res=vast.setsegmentname(16,'Bobbys Dendrite');
% name2=vast.getsegmentname(16);
% 
% res=vast.setsegmentcolor32(16,data.col1,data.col2);
% res=vast.setsegmentcolor8(16,255,0,0,5,108,0,0,0);

% [x,y,z,res] = vast.getviewcoordinates();
% [zoom,res] = vast.getviewzoom();
% res = vast.setviewcoordinates(3805, 9666, 1200);
% res = vast.setviewzoom(0);

%[segimage,res] = vast.getsegimage(0,4800,4900,8100,8200,1100,1110);

%[segdata,res] = vast.getallsegmentdata();
[segname,res] = vast.getallsegmentnames();

% [nr,res]=vast.getnroflayers();
% [layerinfo,res] = vast.getlayerinfo(1);

%Try to disconnect
res=vast.disconnect();
if (res==0)
  warndlg('ERROR: Disconnecting from VAST Lite failed.','Error disconnecting from VAST Lite');
  return;
end

  
