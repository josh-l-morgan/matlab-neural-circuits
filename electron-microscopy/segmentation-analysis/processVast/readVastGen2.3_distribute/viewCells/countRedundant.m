function[redun] = countRedundant(con);

%%Count redundant synapses in connectivity matrix

oneCon = con>0;
difCon = con-oneCon;
redun = sum(difCon>0);