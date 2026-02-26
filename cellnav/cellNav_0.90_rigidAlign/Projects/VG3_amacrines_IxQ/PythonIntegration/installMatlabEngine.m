

pythonPath = 'C:\Users\jlmorgan\.conda\envs\neuron_new\';
addpath(genpath(pythonPath))
matRoot = matlabroot;
cd (fullfile(matlabroot,'extern','engines','python'))
system('python setup.py install')

%{

at windows system prompt
cd "matlabroot\extern\engines\python"
python setup.py install







%}








