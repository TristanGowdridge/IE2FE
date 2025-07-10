clear;clc;close all

% Script Notes
% This script is used to convert the output from the FE models in PEAR into
% the feature data in the PBSHM schema
% Required functions:

% Revision Record
% 01, 11/06/2025, COH, Draft

%Inputs: - Output txt Data File from Lusas 

%Outputs: PBSHM schema complaint feature data files 

%% <<< Inputs >>>

% data_folder = 'C:\Users\conno\Documents\Postdoc - ROSEHIPS\Projects\08_PBSHM_Dataset\BeamSlab\SHMII Dataset BeamSlab\fe_outputs';
% output_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\features';


%% <<<  >>>

tic
% fileList = dir('1-1-*.txt');
fileList = dir('bridge-1-1*.txt');
fileList = {fileList.name};

% fe_fre_all = {};

% count = 1;
for fileRef = 1:length(fileList)
% for fileRef = 1
    FileID = fopen(fileList{1, fileRef});

    % For testing
    % FileID = fopen('bridge-1-1-994_v10_displacement_results.txt');

    % import acceleration data
    Data = textscan(FileID,'%s','delimiter','\t');
    
    % Close the file
    fclose(FileID);

    % extract file creation time and current time
    fileCreationDateID = find(~cellfun(@isempty, regexp(Data{1, 1}, 'Created', 'once')))+1;
    fileCreationDate = Data{1, 1}{fileCreationDateID, 1};
    fileCreationDate = datetime(fileCreationDate,'InputFormat','dd/MM/yy HH:mm');
    CreationTimeStamp = uint64(posixtime(fileCreationDate)*1e+9);
    CurrentTimeStamp = uint64(posixtime(datetime)*1e+9);

    % extract structure name
    StructureName = extractBefore(fileList{1, fileRef},'_');      

    % change the name REMOVE if filename is correct
    StructureName = insertAfter(StructureName,'bridge-1-1','-1-0');

    % Define the pattern: start with one or more digits, a colon, then one or more characters
    pattern = '^\d+:.+$';

    % Apply regexp to each cell; 'once' returns only the first match (or empty if no match)
    LoadCaseID = find(~cellfun(@isempty, regexp(Data{1, 1}, pattern, 'once')));
    LoadCaseID(1) = [];
    LoadCaseID = [LoadCaseID; length(Data{1, 1})+1];

    if contains(fileList{1, fileRef},'ModeShapes','IgnoreCase',true)
        dataType = 'eigenMode';        
        output_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\beam-and-slab\features\eigenMode';
    else
        error('Data type not supported')
    end

    for LoadCasesRef = 1:length(LoadCaseID)-1
        LoadCaseTemp = {Data{1, 1}{LoadCaseID(LoadCasesRef)+1:LoadCaseID(LoadCasesRef+1)-1, 1}};
        % LoadCaseTemp = [0 LoadCaseTemp];
        if mod(LoadCasesRef,2) ~= 0
            LoadCaseTemp = reshape(LoadCaseTemp ,[12, length(LoadCaseTemp)/12])';
        else
            LoadCaseTemp = reshape(LoadCaseTemp ,[9, length(LoadCaseTemp)/9])';  
            continue
        end
        TableNames = LoadCaseTemp(1,:);
        TableNames{1} = 'TableRef';
        valuesUnits = extractBetween(TableNames{1, 6},'[',']');

        LoadCaseTemp = cellfun(@str2double, LoadCaseTemp(2:end,:));
        LoadCaseTemp = array2table(LoadCaseTemp,"VariableNames",TableNames);
        if matches(dataType,'reaction')
             LoadCaseTemp(end,:) = [];
        end
       
        % create json file for current loadcase 'LoadCaseTemp'
        LoadCaseName = lower(Data{1, 1}{LoadCaseID(LoadCasesRef), 1});

        if contains(LoadCaseName,'mode','IgnoreCase',true)
            typeRoot = 'eigen-mode';
            typeHeader = struct('name','eigenMode');
            variantID = 'acceleration';
            ModeShapeID = extractBetween(LoadCaseName,'mode ',' frequency');
            scenarioID = '3';         
            featureName = [typeRoot '-whole-structure-mode-shape' '-' ModeShapeID{1, 1}];
        
        end

        % <<< Header Inputs >>>

        version = "1.3.1";
        population = 'pear-bridge-beam-and-slab-1';         

        % <<< write json file header >>>      

        StructureNameFileName = insertAfter(StructureName,'bridge-1-1-1-0',['-' scenarioID]);
        % fid = fopen([output_folder '\'  StructureName '_' extractAfter(LoadCaseName,':') '_' dataType '.json'],'w');
        fid = fopen([output_folder '\'   StructureNameFileName '-'  featureName '.json'],'w');
        fprintf(fid,'{');
        fprintf(fid,'"version": "%s",', version);
        fprintf(fid,'"name": "%s",', StructureName);
        fprintf(fid,'"population": "%s",', population);
        % fprintf(fid,'"timestamp": %d,', CreationTimeStamp);
        % fprintf(fid,'"features": [');

        % source
        timeStamps = struct('generated', CreationTimeStamp,'stored',CurrentTimeStamp);

        sourceSelection = {struct('name',StructureName)};

        softwareParameters = struct(...
            "verticalDirection", "Z",...
            "analysisCategory", "3D",...
            "modelUnits", "kN,m,t,s,C",...
            "useLineDistance", 'true',...
            "lineMeshSpacing", [2],...
            "lineMeshDistance", [0.5],...
            "surfaceMeshSpacing", [0, 0],...
            "volumeMeshSpacing", [2, 2, 2],...
            "performDynamicAnalysis", 'false',...
            "performStaticAnalysis", 'true',...
            "nModalFeatures", 10,...
            "saveDataDisplacement",'true',...            
            "saveDataReaction", 'true',...
            "saveDataForceMomentBeam", 'false',...
            "saveDataForceMomentShell", 'false',...
            "saveDataLoading", 'false',...
            "saveDataModal", 'false');     

        software1 = {struct('name','IE2FE',...
            'version','0.1',...
            'source','https://github.com/TristanGowdridge/IE2FE',...
            'parameters',softwareParameters)};

        if contains(LoadCaseName,'mode','IgnoreCase',true)            
            environmentDescription = sprintf('Mode Shape %s of the structure, natural frequency in Hz',ModeShapeID{1, 1});
            freID = str2double(extractAfter(LoadCaseName,'frequency = '));
            environmentValue = freID;
            environmentKey1 = 'naturalFrequency';

            environment = struct(environmentKey1,struct('description',environmentDescription,'value',environmentValue));
        end

        % environment = struct(environmentKey1,struct('description',environmentDescription,'value',environmentValue));

        source = struct('timestamp',timeStamps,...
            'selection',{sourceSelection},...
            'software',{software1},...
            'environment',environment);

        fprintf(fid,['"source":' jsonencode(source) ',']);
        
        % data type

        % if contains(LoadCaseName,'mode','IgnoreCase',true)
        %     typeHeader = struct('name','eigenMode');
        %     variantID = 'acceleration';
        % end


        % typeHeader = struct('name','spatial','type',struct('name',dataType));
        % dataHeader = struct('name',[dataType '-' LoadCaseName],'type',typeHeader);
        % fprintf(fid,[jsonencode(dataHeader) ',']);
        % s = jsonencode(source,PrettyPrint=true);

        % fprintf(fid,'"variant": "static",');

        % coordinates
        coordinates_indices = struct('start',0,'end',length(LoadCaseTemp.X)-1);
        x_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.X(1:5));        
        y_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.Y(1:5));
        z_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.Z(1:5));
       
        translational = struct('translational',struct('unit','m','x',x_coordinates,'y',y_coordinates,'z',z_coordinates));

        % coordinates =  struct('coordinates',struct('global',translational));
        % fprintf(fid,[jsonencode(coordinates) ',']);

        % values
        % x_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(1:5,6)));
        % y_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(1:5,7)));
        % z_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(1:5,8)));

        % including the imaginary part

        % these values need to be all in the same precision e.g real:
        % 0.56e-6 imaginary: 0e0,
        x_values_objs = [LoadCaseTemp(1:5,6) table(round(zeros(5,1),4,'significant'))];
        x_values_objs = renamevars(x_values_objs,["DX[m]", "Var1"],["real", "imaginary"]);       
        x_values = struct('indices',coordinates_indices, 'vector',x_values_objs);

        y_values_objs = [LoadCaseTemp(1:5,7) table(zeros(5,1))];
        y_values_objs = renamevars(y_values_objs,["DY[m]", "Var1"],["real", "imaginary"]);
        y_values = struct('indices',coordinates_indices, 'vector',y_values_objs);

        z_values_objs = [LoadCaseTemp(1:5,8) table(zeros(5,1))];
        z_values_objs = renamevars(z_values_objs,["DZ[m]", "Var1"],["real", "imaginary"]);
        z_values = struct('indices',coordinates_indices, 'vector',z_values_objs);

        modeShape = struct('coordinates',struct('global',translational),...
            'values',struct('unit',valuesUnits{1, 1},'x',x_values,'y',y_values,'z',z_values));

        % values = struct('values',struct('unit',valuesUnits{1, 1},'x',x_values,'y',y_values,'z',z_values));
        % fprintf(fid,[jsonencode(values)]);

        % fprintf(fid,'} ] }');
        % featureName = lower(extractAfter(LoadCaseName,':'));
        features = struct('name',featureName,...
            'type',typeHeader,...
            'variant',variantID,...
            'naturalFrequency', struct('unit','Hz','value',freID),...
            'modeShape', modeShape);

        fprintf(fid,['"features": [' jsonencode(features) '] }']);

        fclose(fid);
    end

    if mod(fileRef,100) == 0
        sprintf('Complete:%d',fileRef)
    end
end
fclose('all');
% format default
toc


