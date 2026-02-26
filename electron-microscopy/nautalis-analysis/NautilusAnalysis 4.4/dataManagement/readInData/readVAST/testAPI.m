





vast=VASTControlClass();

res=vast.connect('127.0.0.1',22081,1000);
vast.GETAPIVERSION


if (res==0)
warndlg('Connecting to VAST at 127.0.0.1, port 22081 failed.','Error');
else
vinfo=vast.getinfo();
disp(vinfo)
vast.disconnect();
end;