



%% Example conversion to NWB of the A1407_190416_113437 dataset

% Folder that contains all the files - This is the only input needed
folder_path = 'F:\Adrien';

%% Start Adding fields to NWB

% Get info from the rhd file
rhd = Intan2NWB.GetRHDInfo(folder_path);

% Add general info to the NWB file
nwb = Intan2NWB.GeneralInfo(rhd);

% Add electrode info
nwb = Intan2NWB.getElectrodeInfo(rhd, nwb);

% Add electrophysiological channels
nwb = Intan2NWB.getElectrophysiology(rhd, nwb);

% Add auxiliary channels
nwb = Intan2NWB.getAuxiliary(rhd,nwb);

% Add analogIn channel
nwb = Intan2NWB.getAnalogIn(rhd,nwb);


%% Export to nwb
nwbExport(nwb, 'PoorMouse_converted.nwb')





