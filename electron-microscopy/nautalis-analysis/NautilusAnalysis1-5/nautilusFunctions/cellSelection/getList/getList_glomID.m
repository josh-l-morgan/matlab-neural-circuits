function[axGlom] = getList_glomID(pickGlom)
%%First Column is axon names, second colum is Glomerulus ID 
%%id 1 = A, full reconstruction, 2= H full reconstruction crossover.
%%3,4,5,6 = B,G,K,N partial reconstructions
%%
axGlom = [
2001	1
2030	1
2017	1
2003	1
2004	1
2005	1
2006	1
2007	1
2008	1
2009	1
2010	1
2011	1
2012	3
2013	3
2014	3
2015	3
2016	3
2017	3
2018	3
2013	3
2019	3
2020	3
2032	2
2033	2
2034	2
2035	2
2038	5
2039	5
2040	5
2041	5
2037	6
2021	4
2022	4
2023	4
2024	4
2025	4
2045	4
2026	4
2027	4
2028	4
2029	4
];
%%

if exist('pickGlom','var')
    
    isPicked = zeros(size(axGlom,1),1);
    for i = 1:length(pickGlom)
       isPicked(axGlom(:,2)==pickGlom(i)) = 1; 
    end
    axGlom = axGlom(isPicked>0,1);
end



