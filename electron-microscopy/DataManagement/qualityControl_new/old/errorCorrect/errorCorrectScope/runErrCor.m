
%% Things to do
%{
focTarg Rotation
Check/Set scan rotation
Check inversion


%}


clear all

wif = GetMyWafer;

%% Analyze images
focSec = checkFocus(wif);
'ready for manual input'    
userEval = manFocus(wif,focSec)
retake = shapeOutput(wif,focSec)


if isfield(retake,'tiles')
'found retakes'

%% Activate ActiveX
sm = actxserver('VBComObjectWrapperForZeissAPI.KHZeissSEMWrapperComClass')
sm.InitialiseRemoting
sm.Set_PassedTypeSingle('AP_MAG',25);
sm.Fibics_Initialise();
sprintf(' Fibics Initializing, pausing 15 seconds...')
pause(15)

%%
retakeControl(sm,wif,retake);
end