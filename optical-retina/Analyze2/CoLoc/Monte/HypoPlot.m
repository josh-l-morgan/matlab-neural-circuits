


%Cs = probability of colocalization given synapse
Cr=.4 %Probability of colocalization givin random
Co=.5 %Observed number of colocalizations
Cs=[0:.01:1]

S=(Co-Cr)./(Cs-Cr)
plot(S)
ylim([0 1])