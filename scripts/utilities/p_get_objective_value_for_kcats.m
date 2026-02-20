% Copyright 2019 PSB
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Addition:
% Author: Emanuel Lange
% See explanation in smoment_kcat_optimization.m
% Minimization scenario was not used, so it was removed

function [f] = p_get_objective_value_for_kcats(x)
    % Returns f, the objective value (which shall be minimized) with the given
    % current state of optimized kcat values.
    
    % Arguments:
    % x ~ vector of kcat factors

    % Get given global variables set in smoment_kcat_optimization.m
    global g_cna_model;
    global g_reactions_to_change;
    global g_scenarios;
    global g_scenarios_matrix;
    global g_start_kcats;
    global use_measured_fluxes;

    objective_value = 0;

    % Get changed model
    changed_model = g_cna_model;

    % Changed kcat factors
    kcat_factors = x(1:end);
        
    % original reaction names are stored in the fieldnames of reactions to change
    origReactionNames = convertCharsToStrings(fieldnames(g_reactions_to_change));
    
    % apply kcat factors
    for i = 1:numel(origReactionNames)
        
        % get all coupled split reactions
        origReactionName = origReactionNames(i);
        coupledReactions = g_reactions_to_change.(origReactionName);
        
        % for every original reaction iterate all split reactions and
        % multiply the stoichiometry of pseudo mets with the new change
        % factor
        for j = 1:numel(coupledReactions)
           
            reactionName = coupledReactions(j);
            reaction_index = strmatch(reactionName, changed_model.reacID, 'exact');
        
            new_kcat_factor = kcat_factors(i);
            pseudo_metabolite_indices = get_all_pseudo_metabolite_indices(reactionName, changed_model.stoichMat, changed_model.reacID, changed_model.specID);
            
            for k = 1:numel(pseudo_metabolite_indices)
                metabolite_index = pseudo_metabolite_indices(k);
                changed_model.stoichMat(metabolite_index, reaction_index) = g_start_kcats.(origReactionName).(reactionName)(k) * 1/new_kcat_factor;
            end
            
        end
    end

    % Get values for each scenario line
    number_scenario_lines = size(g_scenarios_matrix, 1);
    fba_error = false;
    sum_scenario_weights = 0;
    for i = 1:number_scenario_lines
        line = g_scenarios_matrix(i,:);
        % Get scenario names \o/
        maximization_scenario_name = line(1);
        maximization_scenario_name = maximization_scenario_name{1};
        % Set maximization variables (used if subsequent minimization is
        % used :-)
        maximization_result = NaN;

        % Apply scenario on model
        scenario = g_scenarios.(maximization_scenario_name);
        scenario_model = p_apply_scenario_on_model(changed_model, scenario);

        % Set model objective :-)
        maximization_target = scenario.target.reaction;
        maximization_target_id = strmatch(maximization_target, scenario_model.reacID, 'exact');
        scenario_model.objFunc(:) = 0;
        scenario_model.objFunc(maximization_target_id) = -1;

        % Perform FBA :D
        [optFlux, success, ~, maximization_result] = CNAoptimizeFlux(scenario_model, [], [], 0, -1);
        maximization_result = maximization_result * -1;
        if (success == 0) || (isnan(maximization_result))
            disp("FBA error (maximization) D:");
            fba_error = true;
        end

        % Calculate objective value
        target_value = scenario.target.value;
        % objective_value_addition = abs( abs(maximization_result) - abs(target_value) );
        objective_value_addition = abs( maximization_result / target_value );
        if (objective_value_addition < 1.0)
            objective_value_addition = 1 / objective_value_addition;
        end
        objective_value_addition = objective_value_addition - 1;

        if strmatch('weight', fieldnames(scenario), 'exact')
            objective_value_addition = objective_value_addition * scenario.weight;
            sum_scenario_weights = sum_scenario_weights + scenario.weight;
        else
            sum_scenario_weights = sum_scenario_weights + 1;
        end

        objective_value = objective_value + abs(objective_value_addition);

    % Return objective value :D
    objective_value = objective_value / sum_scenario_weights;
    if ~fba_error
        f = objective_value;
        fprintf('FBA value: %.3f \n', maximization_result)
        fprintf('fmincon objective: %.3f \n', objective_value)
        disp('______________________')
    else
        f = 100000;
    end
    % sprintf('%.12f', objective_value)
end