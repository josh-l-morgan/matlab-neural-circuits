

%% Grab(xoffset, yoffset, width (pixel) ,height, reduction,filename)
TPN = GetMyDir;
%TPN = 'C:\Images\joshm\testZeiss\';
WriteTo = [TPN 'testGrab.tif'];

sm.Get_ReturnTypeString(
sm.Set_PassedTypeSingle('DP_LINE_SCAN',1)

for i = 1:3
    
WriteTo = [TPN 'testGrab' num2str(i) '.tif'];
sm.Grab(0,0,1024,768,0,WriteTo)
sm.Get_ReturnTypeString('DP_HRRU_PHOTO_STATUS')
end