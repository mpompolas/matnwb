% function convertNWBtoDat()

%                                 nwb_file = ['C:/Users/McGill/Desktop/YutaMouse41150903.nwb'];
                                nwb_file = ['C:\Users\McGill\Documents\GitHub\matnwb\Nas\raw_data_no_header_m120_200secs.nwb'];
nwb2 = nwbRead(nwb_file);


nwb2.acquisition; % This gives the sets that exist in the NWB file


% info = nwb2.acquisition.get('all_lfp').data;
info = nwb2.acquisition.get('ECoG').data;
nChannels = info.dims(2);


% Select how big the chunk will be
ram = 1e9; % 1 GB
max_samples = ram / 8 / nChannels;



%% Convert the acquisition system file to an int16 without a header.
converted_raw_File = ['C:\Users\McGill\Documents\GitHub\matnwb\Nas\' nwb_file '.dat'];

fid = fopen(converted_raw_File, 'a');

isegment = 1;
nsegment_max = 0;

while nsegment_max < info.dims(1)
    nsegment_min = (isegment-1) * max_samples+1;
    nsegment_max = isegment * max_samples - 1;
    if nsegment_max > info.dims(1)
        nsegment_max = info.dims(1);
    end
%     F = in_fread(sFile, ChannelMat, [], [nsegment_min,nsegment_max], [], []);
    F = nwb2.acquisition.get('ECoG').data.load([1,nsegment_min ], [nChannels, nsegment_max]);
    fwrite(fid, F,'int16');
    
    isegment = isegment + 1;
end
fclose(fid);


    
% end
