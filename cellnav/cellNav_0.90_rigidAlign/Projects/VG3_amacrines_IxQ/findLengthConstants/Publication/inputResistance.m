



inputResistanceOhms = 1047 * 10^6; % from Grimes
SE = 477 * 10^6;
STD = SE * sqrt(10);

inputResistanceOhms = inputResistanceOhms /10


%% Calculate areas
um2cm = 10^-6/10^-2;
totalLengthCM = 1500 * um2cm;
cellDiameterCM = 10 * um2cm;
medianDiameterCM = 0.5 * um2cm;
dendCircCM = 2 * pi * medianDiameterCM/2;
dendAreaCM = dendCircCM * totalLengthCM;
cbAreaCM = 4 * pi * (cellDiameterCM/2)^2;
totAreaCM = dendAreaCM + cbAreaCM;


%% Calculate conductance
membraneConductanceTotSiemens = 1/inputResistanceOhms
specificMembraneConductanceSiements = membraneConductanceTotSiemens / totAreaCM;


%%
runCids = [ 2 3 ]

for v = 1:length(runCids)
    runCid = runCids(v)

    swcDir = 'Z:\Active\morganLab\karlsRetina\CellNavLibrary_IxQ\Volumes\Final\Analysis\swc\';
    nepFile = sprintf('%ssm2nrn_cid%d.mat',swcDir,runCid);
    loadNep = load(nepFile);

    eRads = loadNep.sm2nrn.nep.edgeRad;
    eLength = loadNep.sm2nrn.nep.props.edgeLength;

    eRadCM = eRads * um2cm;
    eLengthCM = eLength * um2cm;
    dendCircCM = 2 * pi * eRadCM/2;
    totAreaCM = sum(dendCircCM .* eLengthCM);

    membraneConductanceTotSiemens = 1/inputResistanceOhms;
    specificMembraneConductanceSiements = membraneConductanceTotSiemens / totAreaCM;
    disp(specificMembraneConductanceSiements)
end


