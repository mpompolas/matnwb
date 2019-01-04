



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
sessionInfo.Filename       = name_metadata;% no extension - I DON'T USE THE FILENAME OF THE NWB HERE JUST IN CASE SOMEONE CHANGED IT. 
% sessionInfo.SampleTime     = 50; % 50 no idea
sessionInfo.nElecGps       = []; % 13
sessionInfo.ElecGp         = []; % 1x13 cell (struct with 1x12 cell inside)
% sessionInfo.HiPassFreq     = ;% probably the one from the LFP conversion
% sessionInfo.Date           = 
% sessionInfo.VoltageRange   = % 20
% sessionInfo.Amplification  = % 1000
% sessionInfo.Offset         = 0;
% sessionInfo.lfpSampleRate  = 
% sessionInfo.AnatGrps       =
sessionInfo.spikeGroups.groups        = [];
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
lfp = bz_GetLFP(1:3);


% Load channels 1:3
lfp_interval = bz_GetLFP(1:3,'intervals',[10,20]);




















