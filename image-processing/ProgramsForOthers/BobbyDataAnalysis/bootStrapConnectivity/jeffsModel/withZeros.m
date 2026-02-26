

%{
 there are 753 potential interactions of which 183 that are zeros and the 
rest are synapses. So as one picks axons at random, if it is a zero axon 
it still gets a zero, the rest get synapses and are weighted so that the 
9 hit axon is 9 times more likely to be picked than a one hit or a zero 
hit (which are based on the number of one hit and zero hit axons in the 
bag you pick from). So with replacement keep picking until 753 axons have
been picked.
How does this sound?

synCount = [0:9]
synHist = [183 176 95 26 6 6 4 1 2 1]

%}

