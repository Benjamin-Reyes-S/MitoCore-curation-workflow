function [] = performFBAForScenarios(model, scenarios)

    scenarioNames = string(fieldnames(scenarios));
    fbaResults = zeros(1, numel(scenarioNames));
    measurements = zeros(1, numel(scenarioNames));
        
    for scenarioIndex = 1:numel(scenarioNames)
       
        scenarioName = scenarioNames(scenarioIndex);
        scenarioData = scenarios.(scenarioName);
        
        scenario_model = p_apply_scenario_on_model(model, scenarioData);
        
        [~, success, ~, maximization_result] = CNAoptimizeFlux(scenario_model, [], [], 0, -1);
        
        fbaResults(scenarioIndex) = -maximization_result;
        measurements(scenarioIndex) = scenarioData.target.value;
                
    end
        
    comparisonTable = table(fbaResults', measurements');
    comparisonTable.Properties.VariableNames = {'fba' 'seahorse'};
    
    disp('FBA results and measurements:')
    disp(comparisonTable);
%     
%     maxValue = max([fbaResults measurements]);
%         
%     figure();
%     plot(measurements, fbaResults, 'bx', measurements, measurements);
%     xlim([0 maxValue]);
%     ylim([0 maxValue]);
end

