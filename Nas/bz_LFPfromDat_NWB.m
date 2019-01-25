function bz_LFPfromDat_NWB(basepath,varargin)
% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
%
% note that sincFilter was altered to accomodate GPU filtering
%
%INPUTS
%   basePath    path where the recording files are located
%               where basePath is a folder of the form: 
%                   whateverPath/baseName/
%
%               Assumes presence of the following files:
%                   basePath/baseName.dat
%                   -or-
%                   basePath/amplifier.dat
%
%                   (optional parameters files)
%                   basePath/baseName.xml
%                   basePath/baseName.sessionInfo.mat
%
%               If basePath not specified, tries the current working directory.
%
%   (options)
%       'outFs'         (default: 1250) downsampled frequency of the .lfp
%                       output file. if no user input and not specified in
%                       the xml, use default
%       'lopass'        (default: 450) low pass filter frequency 
%       'noPrompts'     (default: true) prevents prompts about
%                       saving/adding metadata
%
%
%OUTPUT
%   Creates file:   basePath/baseName.lfp
%
%   If no sessionInfo.mat file previously exists, creates one with 
%   the information from the .xml file, with the .lfp sampling frequency 
%   and the lowpass filter used.
%
%
%Dependency: iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
%
%SMckenzie, BWatson, DLevenstein 2018


% Added support for NWB: Konstantinos Nasiotis 2018

%% Input handling
if ~exist('basepath','var')
    basepath = pwd;
end
basename = bz_BasenameFromBasepath(basepath);

defaultoutFS = 1250; %used later

p = inputParser;
addParameter(p,'noPrompts',true,@islogical);
addParameter(p,'outFs',[],@isnumeric);
addParameter(p,'lopass',450,@isnumeric);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
outFs = p.Results.outFs;
lopass = p.Results.lopass;

import iosr.dsp.*

useGPU = false;
% I disable the GPU since I can't find the gather_try function

% try
%     if gpuDeviceCount>0
%         useGPU = true;
%     end
% end
sizeInBytes = 2; % Assumes int16 precision

%% files check
fxml = fullfile(basepath, [basename '.xml']);
fsessioninfo = fullfile(basepath,[basename,'.sessionInfo.mat']);
fdat = fullfile(basepath,[basename,'.dat']);
flfp = fullfile(basepath,[basename,'.lfp']);

fNWB = [basepath '.nwb']; % The .nwb is not inside the new folder


%If there's already a .lfp file, make sure the user wants to overwrite it
if exist(flfp,'file')
    overwrite = input([basename,'.lfp already exists. Overwrite? [Y/N]']);
    switch overwrite
        case {'y','Y'}
            delete(flfp)
        case {'n','N'}
            return
        otherwise
            error('Y or N please...')
    end
end
% % % % % 
% % % % % %Check the dat
% % % % % if ~exist(fdat,'file')
% % % % %     fdat = fullfile(basepath,'amplifier.dat'); %Try amplifier.dat
% % % % %     if ~exist(fdat,'file')
% % % % %         error('Dat file does not exist')
% % % % %     end
% % % % %     
% % % % % end
% % % % % fInfo = dir(fullfile(basepath, [basename '.dat']));

%Get the metadata
if ~exist(fxml,'file') && ~exist(fsessioninfo,'file')
    warning('No xml or sessionInfo file, using defaults and creating a minimal sessionInfo file')
else
    %Get everything from the xml/sessionInfo
    sessionInfo = bz_getSessionInfo(basepath,'noPrompts',noPrompts);
    inFs = sessionInfo.rates.wideband;
    nbChan = sessionInfo.nChannels;
    
    %set output sampling rate from xml, user input    
    if ~isempty(outFs)          %If user input - priority (keep from above)
        outFs = outFs;          %redundant, clearly.
    elseif isfield(sessionInfo,'lfpSampleRate') %If not, use from xml
        outFs = sessionInfo.lfpSampleRate;    
    else                                        %If not in xml, use default
        outFs = defaultoutFS;
    end
end
sessionInfo.lfpSampleRate = outFs;
sessionInfo.LFPLoPassFreq = lopass;
save(fsessioninfo,'sessionInfo');  %Save the sessioninfo with the parameters used


if lopass> outFs/2
    warning('low pass cutoff beyond Nyquist')
end
 

ratio =lopass/(inFs/2) ;
sampleRatio = (inFs/outFs);

%% Set Chunk and buffer size at even multiple of sampleRatio


% % % % Consider using:
% % % % [~,sys] = memory;
% % % % chunksize = sys.SystemMemory.Available * 0.75; % Use 75% of the available memory for demultiplexing. Will speed things up
% % % % clear sys
chunksize = 1e5; % depends on the system... could be bigger I guess





if mod(chunksize,sampleRatio)~=0
    chunksize = chunksize + sampleRatio-mod(chunksize,sampleRatio);
end

%ntbuff should be even multiple of sampleRatio
ntbuff = 525;  %default filter size in iosr toolbox
if mod(ntbuff,sampleRatio)~=0
    ntbuff = ntbuff + sampleRatio-mod(ntbuff,sampleRatio);
end

% nBytes = fInfo.bytes;
nBytes = nbChan * sizeInBytes * sessionInfo.samples_NWB;

nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize));

samplesInEachChunk = floor(sessionInfo.samples_NWB/nbChunks);


%% Check what type of data is installed in the nwb file
% Finding the key/type of data is annoying, check with Ben.

nwb2 = nwbRead(fNWB);

% Check 

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









if sessionInfo.useRaw



    %% GET LFP FROM DAT!
    fidI = fopen(fdat, 'r');
    fprintf('Extraction of LFP begun \n')
    fidout = fopen(flfp, 'a');


    for ibatch = 1:nbChunks

        if mod(ibatch,10)==0
            if ibatch~=10
                fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/nbChunks)), ' percent complete'])]))
            end
            fprintf('%d percent complete', round(100*ibatch/nbChunks));
        end

        h=waitbar(ibatch/(nbChunks+1));


        if ibatch>1
    %         fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');

    %         dat = fread(fidI,nbChan*(chunksize+2*ntbuff),'int16');
    %         dat = reshape(dat,[nbChan (chunksize+2*ntbuff)]);
                dat_nwb = nwb2.acquisition.get(all_raw_keys(iRawDataKey)).data.load([1,((ibatch-1)*chunksize-ntbuff)+1],[nbChan, ((ibatch-1)*chunksize-ntbuff)+1 + chunksize+2*ntbuff]) * 10^(6); % I assume the data is in int16/?V - convert to Volts
        else
    %         dat = fread(fidI,nbChan*(chunksize+ntbuff),'int16');
    %         dat = reshape(dat,[nbChan (chunksize+ntbuff)]);
            dat_nwb = nwb2.acquisition.get(all_raw_keys(iRawDataKey)).data.load([1,1],[nbChan, chunksize+ntbuff]) * 10^(6); % I assume the data is in int16/?V - convert to Volts
        end


    %     DATA = nan(size(dat,1),chunksize/sampleRatio);
        DATA = nan(size(dat_nwb,1),chunksize/sampleRatio);
        for ii = 1:size(dat_nwb,1)

            d = double(dat_nwb(ii,:));
            if useGPU
                d = gpuArray(d);
            end

            tmp=  iosr.dsp.sincFilter(d,ratio);
            if useGPU
                if ibatch==1
                    DATA(ii,:) = gather_try(int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff))));
                else
                    DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff))));
                end

            else
                if ibatch==1
                    DATA(ii,:) = int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff)));
                else
                    DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
                end

            end

        end

        fwrite(fidout,DATA(:),'int16'); 
    end


    % I assume there is no remainder for now. The way that nwb loads data now
    % doesn't need bounds. If the requested value is out of bounds (more
    % samples than those that are actually there, it still loads up to the
    % end).

    % When the way that things load changes, this should be revisited.


    remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunksize;
    if ~isempty(remainder)

    %     fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
    %     dat = fread(fidI,nbChan*(remainder+ntbuff),'int16');
    %     dat = reshape(dat,[nbChan (remainder+ntbuff)]);
        dat_nwb = nwb2.acquisition.get(all_raw_keys(iRawDataKey)).data.load([1,((ibatch-1)*chunksize-ntbuff)+1],[nbChan, ((ibatch-1)*chunksize-ntbuff)+1 + remainder+ntbuff]) * 10^(6); % I assume the data is in int16/?V - convert to Volts

        DATA = nan(size(dat_nwb,1),floor(remainder/sampleRatio));
        for ii = 1:size(dat_nwb,1)
            d = double(dat_nwb(ii,:));
            if useGPU
                d = gpuArray(d);
            end

            tmp=  iosr.dsp.sincFilter(d,ratio);

            if useGPU

                DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end))));
            else
                DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end)));
            end
        end

        fwrite(fidout,DATA(:),'int16');
    end

    close(h);

    % fclose(fidI);
    fclose(fidout);
    
    
    
else % Convert lfp in nwb to .lfp
    
    fidI = fopen(fdat, 'r');
    fprintf('Extraction of LFP begun \n')
    fidout = fopen(flfp, 'a');


    for ibatch = 1:nbChunks

        if mod(ibatch,10)==0
            if ibatch~=10
                fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/nbChunks)), ' percent complete'])]))
            end
            fprintf('%d percent complete', round(100*ibatch/nbChunks));
        end

        h=waitbar(ibatch/(nbChunks+1));
    
        % I inverse the matrix here. I assume that 
        DATA = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.load([1,((ibatch-1)*chunksize)+1],[nbChan, (ibatch*chunksize)]);
    
        fwrite(fidout,DATA(:),'int16'); 
    end


    % I assume there is no remainder for now. The way that nwb loads data now
    % doesn't need bounds. If the requested value is out of bounds (more
    % samples than those that are actually there, it still loads up to the
    % end).

    % When the way that things load changes, this should be revisited.


    remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunksize;
    if ~isempty(remainder)        
        DATA = nwb2.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries.get(all_lfp_keys{iLFPDataKey}).data.load([1, ibatch*chunksize+1],[nbChan, Inf])';    
        fwrite(fidout,DATA(:),'int16');
    end

    close(h);

    % fclose(fidI);
    fclose(fidout);
end
    
    
    
    

disp(' ........baseName.lfp file created! Huzzah!')
end

