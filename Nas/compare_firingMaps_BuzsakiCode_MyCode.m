






%% SHOW THE OUTPUTS OF FIRINGMAP1D FROM BUZSAKI AND MY CODE

%% Buzsaki
cd('C:\Users\McGill\Documents\GitHub\buzcode\tutorials\exampleDataStructs\20170505_396um_0um_merge\')

sessionInfo = bz_getSessionInfo;
lfpChan = sessionInfo.ca1;

% Load spikes,lfp and behavior
spikes = bz_GetSpikes('noprompts',true);
lfp = bz_GetLFP(lfpChan);
behavior = bz_LoadBehavior(pwd,'track');

% mapping spiking onto behavior
[firingMaps] = bz_firingMap1D(spikes,behavior,1,'savemat',false);



%% My code
my_behavior = testing_trials;
[firingMaps_mine] = bz_firingMap1D(spikes,my_behavior,1,'savemat',false);




%%

iCondition =1;
iNeuron = 96

% let's look at some ratemaps...
figure(1)
subplot(2,1,1)
imagesc(squeeze(firingMaps.rateMaps{iCondition}(iNeuron,:,:)))
ylabel('trial #')
xlabel('position')
title ({'Buzsaki Code';['Condition: ' num2str(iCondition)];['Neuron: ' num2str(iNeuron)]})



% let's look at some ratemaps...
figure(1)
subplot(2,1,2)
imagesc(squeeze(firingMaps_mine.rateMaps{iCondition}(iNeuron,:,:)))
ylabel('trial #')
xlabel('position')
title ({'My Code';['Condition: ' num2str(iCondition)];['Neuron: ' num2str(iNeuron)]})





a = squeeze(firingMaps.rateMaps_unsmooth{iCondition}(iNeuron,:,:));
b = squeeze(firingMaps_mine.rateMaps_unsmooth{iCondition}(iNeuron,:,:));

sum(sum(a-b))