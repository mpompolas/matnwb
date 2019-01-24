



%% Create a new sessionInfo with what is needed from the nwb file (just what is necessary for the functions to work - not an exhaustive list for now)

nwb_file = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\m120_25secs.nwb';
nwb2 = nwbRead(nwb_file);

%% Get the type of the recording (its key will be used to get info from the nwb file)

all_keys = keys(nwb2.acquisition);

for iKey = 1:length(all_keys)
    if ismember(all_keys{iKey}, {'ECoG','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
        iDataKey = iKey;
    end
end



[the_path,name_metadata,~] = fileparts(nwb2.acquisition.get(all_keys{iDataKey}).data.filename);
[~,name,~] = fileparts(nwb_file);

if name_metadata ~= name
    warning('The filename on the metadata doesnt correspond to the one that was loaded')
end
clear name


%%
%%%%% MAKE A CHECK HERE IF THE NWB IS THE RAW FILE OR THE CONVERTED TO LFP
% 



%% Create an LFP file from the raw nwb recording
new_path_for_files = [the_path filesep name_metadata];

if ~exist(new_path_for_files)
    mkdir(new_path_for_files)
end


% % If there is already an lfp file, this assumes that a sessionInfo has
% % already been computed.
% flfp = fullfile(the_path,name_metadata,[name_metadata,'.lfp']);
% if ~exist (flfp)


%%
sessionInfo = struct;

sessionInfo.nChannels      = nwb2.acquisition.get(all_keys{iDataKey}).data.dims(2);
sessionInfo.samples_NWB    = nwb2.acquisition.get(all_keys{iDataKey}).data.dims(1);
sessionInfo.nBits          = 16; % ASSUMING THAT NWB GIVES INT16 PRECISION
sessionInfo.rates.lfp      = 1250;    %1250 -  DEFAULT - CHECK THIS
sessionInfo.rates.wideband = nwb2.acquisition.get(all_keys{iDataKey}).starting_time_rate;
sessionInfo.rates.video    = 0;
sessionInfo.FileName       = name_metadata;% no extension - I DON'T USE THE FILENAME OF THE NWB HERE JUST IN CASE SOMEONE CHANGED IT. 
% sessionInfo.SampleTime     = 50; % 50 no idea
sessionInfo.nElecGps       = []; % 13
sessionInfo.ElecGp         = []; % 1x13 cell (struct with 1x12 cell inside)
% sessionInfo.HiPassFreq     = ;% probably the one from the LFP conversion
% sessionInfo.Date           = 
% sessionInfo.VoltageRange   = % 20
% sessionInfo.Amplification  = % 1000
% sessionInfo.Offset         = 0;
sessionInfo.lfpSampleRate  = 1250;    %1250 -  DEFAULT - CHECK THIS
% sessionInfo.AnatGrps       =
sessionInfo.spikeGroups.groups        = [];
sessionInfo.SpkGrps        = []; % I ADDED THIS FOR THE bz_EMGFromLFP
% sessionInfo.channels       =  % 1x128 % starts from 0
% sessionInfo.lfpChans       =
% sessionInfo.thetaChans     =
% sessionInfo.region         = % cell 1x128
% sessionInfo.depth          = 1;   % 3324 - single value
% sessionInfo.ca1            = 116; % 116 - single value
% sessionInfo.ca3            = [] % []
% sessionInfo.ls             = [];
% sessionInfo.animal         = 'MONKEY'% string
% sessionInfo.refChan        = 112; % single value



save([new_path_for_files filesep name_metadata '.sessionInfo.mat'], 'sessionInfo')



%% Get the .lfp file from the converted NWB

bz_LFPfromDat_NWB(new_path_for_files)
% bz_LFPfromDat(new_path_for_files)




%% Check that the lfp is loading from the bz_GetLFP
% bz_GetLFP needs to run within the folder that was created
cd (new_path_for_files)


% Load channels 1:3
lfp = bz_GetLFP(1);


% % Load channels 1:3
% lfp_interval = bz_GetLFP(1:3,'intervals',[10,20]);





%% Create the Ripple and EMG events

        basePath = new_path_for_files;
        [~, baseName] = fileparts(new_path_for_files);
        lfpChan = 1;

        % 2.   Ripple detection

        % Params (defaults listed below)
        ripthresh = [2 5]; %min std, max std
        ripdur = [30 100];
        ripfreq = [130 200];

        % Detect SPW-Rs or load structure with detection info if already detected
        if isempty(dir(fullfile(basePath, '*ripples.events.mat'))) == 1
            ripples = bz_FindRipples(basePath, lfpChan,'thresholds',ripthresh,'durations',ripdur,'show','on','saveMat',true);
        else
            ripples = bz_LoadEvents(basePath,'CA1Ripples');
        end

        % Calculate ripple stats
        [maps,data,stats] = bz_RippleStats(double(lfp.data),lfp.timestamps,ripples);

        % Sort by duration vs amplitude of ripple
        [~,dursort]=sort(data.duration,1,'descend');
        [~,ampsort]=sort(data.peakAmplitude,1,'descend');





%% 3.   Intro spectral processing: filter in ripple frequency and compare to output of bz_FindRipples

% Filter channel in ripple freq range
        ripfiltlfp = bz_Filter(lfp.data,'passband',ripfreq,'filter','fir1'); % lots of options in bz_Filter, read through it
        
%% 4.   Plots 
%% 4i.  Look: Plot some ripples
        figure()
        x=maps.ripples(dursort,:);
        for i=1:100
            subplot(10,10,i)
            plot(x(i,:))
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            axis off
        end
        
%% 4ii. Look: Detected ripples amplitude vs filtered signal (4 different views of the same thing)
        figure()
        subplot 221
        imagesc(maps.amplitude(ampsort,:))
        title('SPW-R Amplitude: sorted by amplitude')
        subplot 222
        imagesc(maps.amplitude(dursort,:))
        title('SPW-R Amplitude: sorted by duration')
        subplot 223
        imagesc(maps.ripples(ampsort,:))
        title('SPW-R Filtered Signal: sorted by amplitude')
        subplot 224
        imagesc(maps.ripples(dursort,:))
        title('SPW-R Filtered Signal: sorted by duration')
        
%% 4iii.Look: Scatterplot with marginals 
        figure()        
        scatterhist(log10(data.duration*1000),data.peakAmplitude,'kernel','on','Location','SouthWest',...
        'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
        'LineWidth',[2,2,2],'Nbins',[20 100], 'marker','.','markersize',10)
        box off
        LogScale('x',10)
        xlabel('logDuration (ms)'); ylabel('Amplitude (au)')
     
%% 5.  Make .evt file for use with neuroscope
        % This is useful for visual inspection of detected events in
        % neuroscope along with the corresponding LFP data. 
        % To load an event file, first open your data. Then, within neuroscope, go to 
        % File > Load Event File, and you will be prompted to select the .evt file.
        % Within neuroscope, a little tab will pop up beside the other two
        % in the top left area - click on that and then click on each of
        % the labels (e.g. start, peak, stop). When they are highlighted,
        % they will be displayed over your data.
       
        % Params
        eventtype = 'ripples'; %string for you that = event type you're saving
        numeventtypes = 3; % you have 3 diff types of events, start, peak, stop
        
        % Save as .evt for inspection
        ripbaseName = [baseName eventtype '.RO1.evt']; % you need the ROX bc neuroscope is buggy and uses this to parse files.
        
        % Below is for ripples specifically 
        % Populate events.time field
        lengthAll = numel(ripples.peaks)*numeventtypes;
        events.time = zeros(1,lengthAll);
        events.time(1:3:lengthAll) = ripples.timestamps(:,1);
        events.time(2:3:lengthAll) = ripples.peaks;
        events.time(3:3:lengthAll) = ripples.timestamps(:,2);
        
        % Populate events.description field
        events.description = cell(1,lengthAll);
        events.description(1:3:lengthAll) = {'start'};
        events.description(2:3:lengthAll) = {'peak'};
        events.description(3:3:lengthAll) = {'stop'};
        
        % Save .evt file for viewing in neuroscope - will save in your current directory
        SaveEvents(fullfile(basePath,ripbaseName),events) %Save and load into neuroscope along with corresponding LFP file




        %% Check that the bz_LoadEvents works


        basePath = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\m120_25secs';

        [ events,filename ] = bz_LoadEvents(basePath);



