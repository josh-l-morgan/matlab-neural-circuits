function runCellNav
clc
addpath(genpath([pwd '\cellNavigator']))
addpath(genpath([pwd '\functions']))
addpath(genpath([pwd '\func']))
addpath(genpath([pwd '\NautilusCN']))
addpath(genpath([pwd '\outsideFunctions']))
addpath(genpath([pwd '\Pictures']))
addpath(genpath([pwd '\bufo']))

makeCellNavGlob

global glob
makeMPNcnv


glob.cellNav.pwd = [pwd '\'];

cellNavGuide_App

disp('Connecting to figure...')
for i = 1:1000
    try
        set(glob.handles.figure1,'clipping','off')
        break
    end
    pause(.1)
end

'hi'

disp('Connected')















