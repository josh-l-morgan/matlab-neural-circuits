function[sm] = updateSMwithNep(targID,nep)

%Update cell information after manual edit


%% load previous SM
load('MPN.mat')
fileName = sprintf('sm_cid%d.mat',targID);
load([WPN 'SMs\' fileName]);

if exist('nep','var')
    sm.nep = nep;
end

%% calculate distances on skel
%sm = getSkelForSM(sm,.5); % up sample edges so that there can be multiple nodes for each current edge
disp('Mapping skeleton distances')
sm = skel2skelDist(sm);
disp('Getting properties of skeletons')
sm = getSkelProps(sm);
disp('Mapping synapse to skeleton distances')
sm = syn2SkelSM(sm);
disp('saving')
save([WPN 'SMs\' sm.fileName],'sm','-v7.3')

%% generate swcS file
disp('Make SWC file')
sm = sm2swc(sm) %write swc file
save([WPN 'SMs\' sm.fileName],'sm','-v7.3')


disp('Save smx file')
smxName = sprintf('smx_cid%d.mat',sm.cid);
save([WPN 'SMs\' smxName],'-struct','sm','nep','syn',...
    'cell','skel2skel','syn2Skel','-v7.3');
disp('finished saving smx')




