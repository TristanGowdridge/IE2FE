clear;clc;close all

% Script Notes
% This script is used to convert the output from the FE models in PEAR into
% the feature data in the PBSHM schema
% Required functions:N/A

% Revision Record
% 01, 11/06/2025, COH, Draft
% 02, 09/07/2025, COH, Added reading of SI models to get loading
% locations

%Inputs: - Output txt Data File from Lusas
% si models

%Outputs: PBSHM schema complaint feature data files

%% <<< Inputs >>>

% data_folder = 'C:\Users\conno\Documents\Postdoc - ROSEHIPS\Projects\08_PBSHM_Dataset\BeamSlab\SHMII Dataset BeamSlab\fe_outputs';
% output_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\features';
% load('C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\BeamSlabDataSetSHMII.mat')
si_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\si-models';

%% <<<  >>>

tic
% fileList = dir('1-1-*.txt');
fileList = dir('bridge-1-1*.txt');
fileList = {fileList.name};

% fe_fre_all = {};

% count = 1;
for fileRef = 1:length(fileList)
% for fileRef = 1:2
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

    % read in si model data to get position of load
    siFileName = [si_folder '\' StructureName '.json'];
    siData = fileread(siFileName);
    siData = jsondecode(siData); % Using the jsondecode function to parse JSON from string

    span_lengths = siData.models.structuralInformation.parameters.span.length.values;
    bridge_width = (siData.models.structuralInformation.parameters.beams.quantity-1)*...
        siData.models.structuralInformation.parameters.beams.center_to_center.value;

    x_temp = cumsum(span_lengths)';
    x_temp = [0 x_temp];

    centre_span_location = [];
    mid_y = bridge_width/2;
    for i = 1:length(span_lengths)
        mid_x = mean([x_temp(i) x_temp(i+1)]);
        centre_span_location(i,1:3) = [mid_x mid_y 0];
    end

    % Define the pattern: start with one or more digits, a colon, then one or more characters
    pattern = '^\d+:.+$';

    % Apply regexp to each cell; 'once' returns only the first match (or empty if no match)
    LoadCaseID = find(~cellfun(@isempty, regexp(Data{1, 1}, pattern, 'once')));
    LoadCaseID(1) = [];
    LoadCaseID = [LoadCaseID; length(Data{1, 1})+1];

    if contains(fileList{1, fileRef},'displacement','IgnoreCase',true)
        dataType = 'displacement';
        output_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\features\dead-load-mid-span-displacement';
    elseif contains(fileList{1, fileRef},'reaction','IgnoreCase',true)
        dataType = 'reaction';
        LoadCaseID = LoadCaseID(1:2);
        output_folder = 'C:\Users\conno\Documents\ROSEHIPS\PEAR UoS July 25\features\dead-load-reaction';
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

        if contains(LoadCaseName,'deadload','IgnoreCase',true)
            typeRoot = 'spatial';
            typeHeader = struct('name',typeRoot,'type',struct('name',dataType));
            variantID = 'static';
            if matches(dataType,'displacement')
                scenarioID = '1';
            elseif matches(dataType,'reaction')
                scenarioID = '2';
            end
            featureName = [typeRoot '-' dataType '-whole-structure-dead-load'];
        elseif contains(LoadCaseName,'Span','IgnoreCase',true)
            spanID = extractBetween(LoadCaseName,'span','_');
            typeRoot = 'spatial';
            typeHeader = struct('name',typeRoot,'type',struct('name',dataType));
            variantID = 'static';
            scenarioID = '3';
            featureName = [typeRoot '-' dataType '-whole-structure' '-' spanID{1, 1} '-midspan-point-load'];
        end

        % <<< Header Inputs >>>

        version = "1.3.0";
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

        if contains(LoadCaseName,'deadload','IgnoreCase',true)
            environmentDescription = 'Structure is loaded under own weight in z axis in m/s^2' ;
            environmentValue = 9.81;
            environmentKey1 = 'load';

            environment = struct(environmentKey1,struct('description',environmentDescription,'value',environmentValue));
        elseif contains(LoadCaseName,'Span','IgnoreCase',true)

            environmentDescription = sprintf('Point load at midspan of span %s in -z axis in kN',spanID{1, 1});
            environmentValue = 400;
            environmentKey1 = 'loadValue';
            environmentKeyStruct1 = struct('description',environmentDescription,'value',environmentValue);

            environmentDescription = 'X coordinates of load position in meters';
            environmentValue = centre_span_location(str2double(spanID{1, 1}),1);
            environmentKey2 = 'loadLocationX';
            environmentKeyStruct2 = struct('description',environmentDescription,'value',environmentValue);

            environmentDescription = 'Y coordinates of load position in meters';
            environmentValue = centre_span_location(str2double(spanID{1, 1}),2);
            environmentKey3 = 'loadLocationY';
            environmentKeyStruct3 = struct('description',environmentDescription,'value',environmentValue);

            environmentDescription = 'Z coordinates of load position in meters';
            environmentValue = centre_span_location(str2double(spanID{1, 1}),3);
            environmentKey4 = 'loadLocationZ';
            environmentKeyStruct4 = struct('description',environmentDescription,'value',environmentValue);

            environment = struct(environmentKey1,environmentKeyStruct1,...
                environmentKey2,environmentKeyStruct2,...
                environmentKey3,environmentKeyStruct3,...
                environmentKey4,environmentKeyStruct4);
        end

        % environment = struct(environmentKey1,struct('description',environmentDescription,'value',environmentValue));

        source = struct('timestamp',timeStamps,...
            'selection',{sourceSelection},...
            'software',{software1},...
            'environment',environment);

        fprintf(fid,['"source":' jsonencode(source) ',']);

        % data type

        % if contains(LoadCaseName,'deadload','IgnoreCase',true)
        %     typeRoot = 'spatial';
        %     typeHeader = struct('name',typeRoot,'type',struct('name',dataType));
        %     variantID = 'static';
        % elseif contains(LoadCaseName,'Span','IgnoreCase',true)
        %     typeRoot = 'spatial';
        %     typeHeader = struct('name',typeRoot,'type',struct('name',dataType));
        %     variantID = 'static';
        % end


        % typeHeader = struct('name','spatial','type',struct('name',dataType));
        % dataHeader = struct('name',[dataType '-' LoadCaseName],'type',typeHeader);
        % fprintf(fid,[jsonencode(dataHeader) ',']);
        % s = jsonencode(source,PrettyPrint=true);

        % fprintf(fid,'"variant": "static",');

        % coordinates
        coordinates_indices = struct('start',0,'end',length(LoadCaseTemp.X)-1);
        x_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.X);
        y_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.Y);
        z_coordinates = struct('indices',coordinates_indices, 'vector',LoadCaseTemp.Z);


        translational = struct('translational',struct('unit','m','x',x_coordinates,'y',y_coordinates,'z',z_coordinates));

        % coordinates =  struct('coordinates',struct('global',translational));
        % fprintf(fid,[jsonencode(coordinates) ',']);

        % values
        x_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(:,6)));
        y_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(:,7)));
        z_values = struct('indices',coordinates_indices, 'vector',table2array(LoadCaseTemp(:,8)));

        % values = struct('values',struct('unit',valuesUnits{1, 1},'x',x_values,'y',y_values,'z',z_values));
        % fprintf(fid,[jsonencode(values)]);

        % fprintf(fid,'} ] }');
        % featureName = lower(extractAfter(LoadCaseName,':'));
        % name ref spatial-displacement-wholestructure-1-midspan-point-load

        features = struct('name',featureName,...
            'type',typeHeader,...
            'variant',variantID,...
            'coordinates',struct('global',translational),...
            'values',struct('unit',valuesUnits{1, 1},'x',x_values,'y',y_values,'z',z_values));

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


