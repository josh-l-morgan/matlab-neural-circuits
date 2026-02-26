%%Calculate relative contribution of proportion of sements with puncta
%%compared to density of puncta on segments. 
%%Must be run after God-whole
%%DD = NZ * Fill (where DD is full density, NZ is the density of non zero
%%sements and Fill is the proportion of nonzero segements


NZ7=mean(meanNonZeDDs(Age==7))
NZ30=mean(meanNonZeDDs(Age>30))

DD7=mean(DD(Age==7))
DD30=mean(DD(Age>30))

Fill=1-HistDC(:,1);
Fill7=mean(Fill(Age==7))
Fill30=mean(Fill(Age>30))

deltaDD=DD30-DD7
deltaFill=Fill30-Fill7
deltaNZ=NZ30-NZ7

ratDD=DD30/DD7
ratFill=Fill30/Fill7
ratNZ=NZ30/NZ7

PercentContributionByLocalDensity=(ratDD-ratFill)/(ratDD+ratFill-2)