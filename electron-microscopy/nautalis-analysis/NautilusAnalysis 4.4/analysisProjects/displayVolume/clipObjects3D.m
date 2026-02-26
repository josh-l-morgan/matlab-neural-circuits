





%%

clip = [1000 1800 100 ;  1150 1951 300]; 


cellTarg = obI.cell.name == 4;
obIds = obI.cell.obIDs{cellTarg};

allSubs = double(uniqueSubs(cat(1,dsObj(obIds).subs)));
[subsIn subsOut] = clipSubs(allSubs,clip);
clf
renderCon(subsIn,[],[1 0 0],.5)
renderCon(subsOut,[],[1 1 1],.2)



%%







