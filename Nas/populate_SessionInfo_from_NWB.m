



%% Create a new sessionInfo with what is needed from the nwb file (just what is necessary for the functions to work - not an exhaustive list for now)

% nwb_file = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\m120_25secs.nwb';
nwb_file = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41.nwb';
nwb2 = nwbRead(nwb_file);

%% Get the type of the recording (its key will be used to get info from the nwb file)
% There will be 2 logical values:
% RawDataPresent
% LFPDataPresent

% A conversion to .lfp files will be prioritized on the raw data if both
% datatypes are present


[the_path, name,~] = fileparts(nwb_file);

% First check if there are raw data present

try
    all_raw_keys = keys(nwb2.acquisition);

    for iKey = 1:length(all_raw_keys)
        if ismember(all_raw_keys{iKey}, {'ECoG','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
            iRawDataKey = iKey;
            RawDataPresent = 1;
        else
            RawDataPresent = 0;
        end
    end
    % nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get('all_lfp').data
catch
    RawDataPresent = 0;
end




try
    % Check if the data is in LFP format
    all_lfp_keys = keys(nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries);

    for iKey = 1:length(all_lfp_keys)
        if ismember(all_lfp_keys{iKey}, {'all_lfp','bla bla bla'})   %%%%%%%% ADD MORE HERE, DON'T KNOW WHAT THE STANDARD FORMATS ARE
            iLFPDataKey = iKey;
            LFPDataPresent = 1;
            break % Once you find the data don't look for other keys/trouble
        else
            LFPDataPresent = 0;
        end
    end
catch
    LFPDataPresent = 0;
end


if ~RawDataPresent && ~LFPDataPresent
    error 'There is no data in this .nwb - Maybe check if the Keys are labeled correctly'
end


%% Create an LFP file from the raw nwb recording
new_path_for_files = [the_path filesep name];

if ~exist(new_path_for_files)
    mkdir(new_path_for_files)
end


% % If there is already an lfp file, this assumes that a sessionInfo has
% % already been computed.
% flfp = fullfile(the_path,name_metadata,[name_metadata,'.lfp']);
% if ~exist (flfp)


%%
sessionInfo = struct;

if RawDataPresent
    sessionInfo.nChannels      = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).data.dims(2);
    sessionInfo.samples_NWB    = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).data.dims(1);
    sessionInfo.rates.wideband = nwb2.acquisition.get(all_raw_keys{iRawDataKey}).starting_time_rate;
    sessionInfo.rates.lfp      = 1250;    %1250 -  DEFAULT - CHECK THIS
    sessionInfo.lfpSampleRate  = 1250;    %1250 -  DEFAULT - CHECK THIS % tH lfpSampleRate bypasses 

elseif LFPDataPresent
    sessionInfo.nChannels      = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.dims(2);
    sessionInfo.samples_NWB    = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.dims(1);
    sessionInfo.rates.wideband = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;
    
    sessionInfo.rates.lfp      = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;    %I assign the LFP sampling rate that was already used. Not sure yet if a different value than 1250 Hz will cause problems
    sessionInfo.lfpSampleRate  = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).starting_time_rate;    %I assign the LFP sampling rate that was already used. Not sure yet if a different value than 1250 Hz will cause problems
    
    if sessionInfo.rates.wideband <1250
        warning 'Something weird will happen. Not thoroughly tested yet'
    end
end





sessionInfo.nBits          = 16; % ASSUMING THAT NWB GIVES INT16 PRECISION
sessionInfo.rates.video    = 0;
sessionInfo.FileName       = name;% no extension - I DON'T USE THE FILENAME OF THE NWB HERE JUST IN CASE SOMEONE CHANGED IT. 
% sessionInfo.SampleTime     = 50; % 50 no idea
sessionInfo.nElecGps       = []; % 13
sessionInfo.ElecGp         = []; % 1x13 cell (struct with 1x12 cell inside)
% sessionInfo.HiPassFreq     = ;% probably the one from the LFP conversion
% sessionInfo.Date           = 
% sessionInfo.VoltageRange   = % 20
% sessionInfo.Amplification  = % 1000
% sessionInfo.Offset         = 0;
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

sessionInfo.useRaw           = RawDataPresent; % I add this to choose which data to convert to .lfp files

save([new_path_for_files filesep name '.sessionInfo.mat'], 'sessionInfo')



%% Get the .lfp file from the converted NWB



bz_LFPfromDat_NWB(new_path_for_files)
% bz_LFPfromDat(new_path_for_files)




%% Check that the lfp is loading from the bz_GetLFP
% bz_GetLFP needs to run within the folder that was created


            new_path_for_files = 'C:\Users\McGill\Documents\GitHub\matnwb\Nas\YutaMouse41' ;



cd (new_path_for_files)


% Load channel 1 % Starts from 0 in bz_GetLFP
ichannel = 0;

lfp = bz_GetLFP(ichannel);
figure(1);
plot(lfp.data)

% % Load channels 1:3
% lfp_interval = bz_GetLFP(1:3,'intervals',[10,20]);


% compare to nwb
data_nwb_channel = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.load([ichannel+1,1],[ichannel+1, Inf]);

figure(2); plot(data_nwb_channel)

%% Create the Ripple and EMG events


        % lfp must be from only 1 channel for this to work!
        % Make sure not multiple channels are loaded on the previous
        % section from bz_GetLFP



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


        
        
        
        
        
  %% Work on the spikes      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

