
%% Activate ActiveX
sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
sm.InitialiseRemoting

%% Command Format
%{
sm.GetMag
sm.Get_ReturnTypeSingle('AP_MAG')
sm.Get_ReturnTypeString('DP_COLUMN_TYPE')
sm.Set_PassedTypeSingle('AP_MAG',100)
%}

%% Test Stage Functions
oldX = sm.Get_ReturnTypeSingle('AP_STAGE_AT_X');
oldY = sm.Get_ReturnTypeSingle('AP_STAGE_AT_Y');


goX = oldX + .00010;
sm.Set_PassedTypeSingle('AP_STAGE_GOTO_X',goX)
goY = oldY + .00010;
sm.Set_PassedTypeSingle('AP_STAGE_GOTO_Y',goY)

%% Test Focus Functions
oldWD = sm.Get_ReturnTypeSingle('AP_WD');
goWD = oldWD + .0001;
sm.Set_PassedTypeSingle('AP_WD',goWD)

oldSX = sm.Get_ReturnTypeSingle('AP_STIG_X');
oldSY = sm.Get_ReturnTypeSingle('AP_STIG_Y');

sm.Execute('CMD_AUTO_FOCUS_FINE')
dpAuto = 'holding';
while ~strcmp(dpAuto,'Idle')  %wait to finish focusing
   pause(.01)
   dpAuto = sm.Get_ReturnTypeString('DP_AUTO_FUNCTION');   
end



autoWD = sm.Get_ReturnTypeSingle('AP_WD');

difWD = autoWD - oldWD
pause(3)
autoWD = sm.Get_ReturnTypeSingle('AP_WD');
difWD = autoWD - oldWD
pause(10)
autoWD = sm.Get_ReturnTypeSingle('AP_WD');
difWD = autoWD - oldWD
sm.Execute('CMD_AUTO_STIG')
dpAuto = sm.Get_ReturnTypeString('DP_AUTO_FUNCTION')


%%
dpAuto = sm.Get_ReturnTypeString('DP_AUTO_FUNCTION');


