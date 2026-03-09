%% main script for calibration

%% cleaning
close all;
clear all;
clc;

%% Adding project to path
editorPath = matlab.desktop.editor.getActiveFilename;
scriptFolder = fileparts(editorPath);
targetFolder = fileparts(scriptFolder);

addpath(genpath(targetFolder));

%% define file paths
cnaPath = "C:\Users\emanuel.lange\CNA2025.2_ca7\CellNetAnalyzer"; % download CNA and define path to it

modelStruct.Mitocore_Preliminary = ".\autopacmen_output\Mitocore_Preliminary\Mitocore_Preliminary_uncalibrated.xml";
modelStruct.Mitocore_MitoMammal = ".\autopacmen_output\Mitocore_MitoMammal\Mitocore_MitoMammal_uncalibrated.xml";
modelStruct.Mitocore_aligned_to_Human1 = ".\autopacmen_output\Mitocore_aligned_to_Human1\Mitocore_aligned_to_Human1_uncalibrated.xml";

candidateReactionsPath =  ".\autopacmen_output\kcat_candidates.json";
scenariosPath = ".\model_data_input\platelet_scenario.json";
calibratedParamsPath = ".\autopacmen_output\calibrated_factors.json";

%% start cna
addpath(cnaPath);

startcna

%% read sceanrio struct
scenarioStruct = fileread(scenariosPath);
scenarioStruct = jsondecode(scenarioStruct);

%% get candidate reactions
candidatesString = fileread(candidateReactionsPath);
candidatesStruct = jsondecode(candidatesString);

%% perform calibration for each
global currentModelName

allCalibratedParams = struct([]);

modelNames = fieldnames(modelStruct);

for modelIndex = 1 : numel(modelNames)
   
   modelNameChar = string(modelNames(modelIndex)); 
   candidateReactions = string(fieldnames(candidatesStruct.(modelNameChar).reactions));
   
   disp("====================================");
   disp(modelNameChar);
   
   currentModelName = modelNameChar;
   
   scenarioMatrix = {"platelet", "NA"};
   
   % read max change factor
   maxChangeFactor = candidatesStruct.(modelNameChar).max_change_factor;
      
   % char bakes a string, otherwise CNA complains..
   modelPath = char(modelStruct.(modelNameChar));
   
   calibratedParams = calibrate_05(modelPath, candidateReactions, scenarioStruct, scenarioMatrix, maxChangeFactor);
      
   allCalibratedParams(1).(modelNameChar) = calibratedParams;
   
   % visualize predictions
    [modelUncalibrated, ~] = CNAsbmlModel2MFNetwork(modelPath);

    calibratedModel = applyKcatsToModel(modelUncalibrated, calibratedParams.reactions);
    performFBAForScenarios(calibratedModel, scenarioStruct);
end

%% write calibrated parameters to json
jsonString = jsonencode(allCalibratedParams);
fid = fopen(calibratedParamsPath,'w');
fprintf(fid, '%s', jsonString);
fclose(fid);