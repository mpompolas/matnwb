%% nwbFiles - You may download from the browser instead - I have it here so we can all use the same folders
path_to_download_files = 'C:\Users\knasi\Desktop\nwb\Files';
fileURLs = {'https://girder.dandiarchive.org/api/v1/item/5fa438fa271e91aad754e4e9/download'; % sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys 2.2.5
            'https://girder.dandiarchive.org/api/v1/item/5eda806c99f25d97bd279816/download'; % sub-703279277_ses-719161530_probe-729445648_ecephys.nwb  2.2.2
            };

fileNames = {'sub-YutaMouse41_ses-YutaMouse41-150829_behavior+ecephys.nwb';
             'sub-703279277_ses-719161530_probe-729445648_ecephys.nwb'
             };
     
mkdir(path_to_download_files)
%Download files
for i = 1:length(fileURLs) 
    fullFilePath = fullfile(path_to_download_files, fileNames{i});
    if ~isfile(fullFilePath)
        outFileNames{i} = websave(fullFilePath, fileURLs{i});
    else
        outFileNames{i} = fullFilePath;
    end
end
    

%% MatNWB - download and extract versions of MatNWB
path_to_download_MatNWB = 'C:\Users\knasi\Desktop\nwb\MatNWB';



matNWBURLs = {'https://github.com/NeurodataWithoutBorders/matnwb/archive/v2.2.5.1.zip';
              'https://github.com/NeurodataWithoutBorders/matnwb/archive/v2.2.5.0.zip';
              'https://github.com/NeurodataWithoutBorders/matnwb/archive/v2.2.4.0.zip';
              'https://github.com/NeurodataWithoutBorders/matnwb/archive/0.2.3.zip';
              'https://github.com/NeurodataWithoutBorders/matnwb/archive/v0.2.2.zip'};

MatNWBNames = {'2.2.5.1';
               '2.2.5.0';
               '2.2.4.0';
               '0.2.3'; % 2.2.2
               '0.2.2'}; % 2.2.1

%Download
mkdir(path_to_download_MatNWB)
for i= 1:length(matNWBURLs)
    cd(path_to_download_MatNWB)
    fullMatNWBPath = fullfile(path_to_download_MatNWB,MatNWBNames{i});
    outFileNameMatNWB{i} = websave(fullMatNWBPath, matNWBURLs{i});
    
    % Unzip downloaded MatNWB
    unzip(outFileNameMatNWB{i})
    % Final folder that has all the matNWB files
    finalMatNWBFolder{i} = fullfile(path_to_download_MatNWB, ['matnwb-' MatNWBNames{i}]);
end





%% Check the version of the NWB file and load the appropriate matnwb
for iFile = 1:length(outFileNames)
    
    % GenerateCore within 2.2.5.0 so you have access to util.getSchemaVersion
    cd(finalMatNWBFolder{2})
    generateCore() % This fails here for 2.2.5.1 - HOW GO AROUND THIS?
	fileSchemaVersion = util.getSchemaVersion(outFileNames{iFile});
    rmpath(genpath(finalMatNWBFolder{2})); % Remove from path what was used to get the util.getSchemaVersion
    
    % Add what is needed for changing version here and then call nwbRead
    if strcmp(fileSchemaVersion, '2.2.1')
        associated_folder_string = '0.2.2';
    elseif strcmp(fileSchemaVersion, '2.2.2')
        associated_folder_string = '0.2.3';
    elseif strcmp(fileSchemaVersion, '2.2.4')
        associated_folder_string = '2.2.4.0';
    elseif strcmp(fileSchemaVersion, '2.2.5')
        associated_folder_string = '2.2.5.0';
    else
        error('version is not supported by Brainstorm')
    end
    
    % Index of cached schema that correspond to the fileVersion
    iSchema = find(ismember(MatNWBNames, associated_folder_string));
    % Enter the folder that corresponds to that version
    cd(finalMatNWBFolder{iSchema})
    generateCore()

    
    %%%%%%%%%%%% ADD HERE %%%%%%%%%%%%%%
   
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nwb2 = nwbRead(outFileNames{iFile});
    disp('Success')
%     rmpath(genpath(finalMatNWBFolder{iSchema})); % Remove from path what was used to load the file

end

    
    
    
    
    