function runCellNav





addpath(genpath([pwd '\cellNavigator']))
addpath(genpath([pwd '\functions']))
addpath(genpath([pwd '\NautilusCN']))
addpath(genpath([pwd '\outsideFunctions']))



makeCellNavGlob
global glob
makeMPNcnv

glob.cellNav.pwd = [pwd '\'];


cellNavGuide_37

for i = 1:1000
    try
        set(glob.handles.figure1,'clipping','off')
        break
    end
    pause(.1)
end

'hi'
















