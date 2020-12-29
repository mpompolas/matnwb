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
    

%% Download matNWB 2.2.5.0
path_to_download_MatNWB = 'C:\Users\knasi\Desktop\nwb\MatNWB';

matNWBURLs = {'https://github.com/NeurodataWithoutBorders/matnwb/archive/v2.2.5.0.zip'};
MatNWBNames = {'2.2.5.0'};   

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

%% Schema versions

path_to_download_schemas = 'C:\Users\knasi\Desktop\nwb\schemas';



schemaURLs = {'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.5.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.4.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.3.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.2.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.1.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.2.0.zip';
              'https://github.com/NeurodataWithoutBorders/nwb-schema/archive/2.1.0.zip'};

schemaNames = {'2.2.5';
               '2.2.4';
               '2.2.3';
               '2.2.2';
               '2.2.1';
               '2.2.0';
               '2.1.0'};

%Download
mkdir(path_to_download_schemas)
for i= 1:length(schemaURLs)
    cd(path_to_download_schemas)
    fullschemaPath = fullfile(path_to_download_schemas,schemaNames{i});
    outFileNameSchema{i} = websave(fullschemaPath, schemaURLs{i});
    
    % Unzip downloaded MatNWB
    unzip(outFileNameSchema{i})
    % Final folder that has all the matNWB files
    finalSchemaFolder{i} = fullfile(path_to_download_schemas, ['nwb-schema-' schemaNames{i}]);
end


%% Check the version of the NWB file and load the appropriate matnwb
matNWBpath = finalMatNWBFolder{1}; % MatNWB 2.2.5.0
cd(matNWBpath)
generateCore()
addpath(matNWBpath)


for iFile = 1:length(outFileNames)
    
	fileSchemaVersion = util.getSchemaVersion(outFileNames{iFile});    
    
    % Index of cached schema that correspond to the fileVersion
    iSchema = find(ismember(schemaNames, fileSchemaVersion));
    
    if isempty(iSchema)
        error('This version is not supported by Brainstorm')
    end
    
    % Entering and adding schema folders to the path
    cd(finalSchemaFolder{iSchema})
    addpath(genpath(finalSchemaFolder{iSchema}))

    % Read and then remove the schema version from the path
    nwb2 = nwbRead(outFileNames{iFile});
    disp('Success')
    rmpath(genpath(finalSchemaFolder{iSchema})); % Remove from path what was used to load the file

end

    
    
    
    
    