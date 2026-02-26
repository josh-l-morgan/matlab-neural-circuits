function appoLocs=getAppos(group1,group2)
%connAppStruct=struct();
appRad=0.1;
appDensity=1;
global glob globFA
globFA.group1.idx=find(ismember(glob.cids,group1));
globFA.group2.idx=find(ismember(glob.cids,group2));
globFA.appDist=appRad;
globFA.outputGrid=appDensity;
%conAppStruct.appo=runFindAppositions();
appoLocs=runFindAppositions();


end
