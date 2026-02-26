

SPN = 'E:\IxQ_KarlsRetinaVG3_2019\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\swc\'
nrnFile = 'cid3323.npy';
nrnFile = 'nrn_cid3323.mat'

nrn = load([SPN nrnFile])

if 0
    pe = pyenv;
    if pe.Status == "NotLoaded"
        [~,exepath] = system("where python");
        pe = pyenv('Version',exepath);
    end
end

pe = pyenv('Version','C:\Users\jlmorgan\.conda\pkgs\python-3.8.12-h6244533_0\python.exe')
pyLibraryFolder2 = [pwd '\Python'];
insert(py.sys.path, int64(0),pyLibraryFolder2);

pe.Version
py.print("hello, Python?")

kwa = pyargs('cids',4)

cal = py.calendar.TextCalendar

formatmonth(cal,int32(2014),int32(12))

py.runNrnFromSwcSm

a = py.testFunction



