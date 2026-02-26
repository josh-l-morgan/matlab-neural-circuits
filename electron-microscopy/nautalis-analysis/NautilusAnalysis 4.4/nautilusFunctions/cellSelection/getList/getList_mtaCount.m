function[checkID cellProp] = getList_mtaCount();

%%Returns a list of cells near the seed cell and the number of manually
%%identified multitubular arrays in their cell body

mtaCount = [108	1
511	0
177	0
115	1
400	0
201	0
114	2
120	0
267	2
506	0
901	0
210	0
905	0
273	0
109	3
144	0
271	3
903	0
908	1
151	2
906	2
925	1
129	4
508	2
902	1
907	1
217	0
203	0
232	1
224	1
241	2
221	0
230	1];

checkID = mtaCount(:,1);
cellProp = mtaCount(:,2);