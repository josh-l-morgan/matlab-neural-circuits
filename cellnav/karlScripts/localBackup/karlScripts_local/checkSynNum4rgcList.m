allVG3cidList=[2 3 4 5 10 11 13 14 20];
rgcList=[1110,1124,1180,2002,2003,3007,3008,3049,3223,1082,1124,1178 ...
    ,1187,1191,1209,1254,3171];
results=zeros(length(rgcList),2);
for p=1:length(rgcList)
    curResults=find(curTis.syn.edges(:,1)==rgcList(p)&ismember(curTis.syn.edges(:,2),allVG3cidList));
    results(p,:)=[rgcList(p) length(curResults)];
end