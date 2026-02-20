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
% autoPACMEN was originally only calibrating the stoichiometric
% coefficients M/kcat of the global pool metabolite (M_prot_pool) in
% reactions and did not consider pools for individual proteins and their 
% stoichiometric coefficients (1/kcat). Instead of changing the
% stoichiometric coefficients directly, I introduced a factor
% (start_kcat_factors) that scales the stoichiometric coefficients of
% pseudometabolites in the objective functions. Effectively, this factor is
% calibrated. The advantage here is that all stoichiometric coefficients of
% pseudo mets are scaled at the same time and there is no need to keep 
% track of them seperately(they are all based on the kcat).

function [best_kcat_factors, start_kcat_factors] = smoment_kcat_optimization(cna_model, reactions_to_change, scenarios, scenarios_matrix, max_change_factor)
    % AutoPACMEN kcat model calibrator.
    %
    % Calibrates the stoichiometric coefficients (M/kcat and 1/kcat)
    % so that the model output fits best with the given scenarios.
    % For a more detailed explanation on scenarios and this
    % function's context within AutoPACMEN, see AutoPACMEN's manual.
    %
    % Usage example:
    % [best_kcats, start_kcats] = smoment_kcat_optimization(cna_model,...
    %                                                       reactions_to_change,...
    %                                                       scenarios,...
    %                                                       scenarios_matrix,...
    %                                                       max_change_factor)
    %
    % Arguments:
    % cna_model ~ A CellNetAnalyzer model, e.g. loaded with
    %             CNAsbmlModel2MFNetwork
    % reactions_to_change ~ A struct containing original reaction 
    %                   name + direction and for every original reaction,
    %                   the split reactions.
    % scenarios ~ A struct from a scenarios JSON loaded with jsondecode.
    %             See the manual for more on this.
    % scenarios_matrix ~ The scenario dependency matrix. See the manual for
    %                    more on this.
    % max_change_factor ~ Number; The kcats will be restriced to be
    %                     minimally max_change_factor times lower or
    %                     max_change_factor times higher than their start
    %                     value.
    %
    % Outputs:
    % best_kcats ~ A list of optimized MW/kcat values in the order of
    %              the reactions given in reactions_to_change.
    % start_kcats ~ The unoptimized MW/kcat values in the order of
    %               the reactions given in reactions_to_change. Useful in
    %               order to compare them with best_kcats.
    
    % Set global variables used in fmincon
    global g_cna_model;
    g_cna_model = cna_model;
    global g_reactions_to_change
    g_reactions_to_change = reactions_to_change;
    global g_scenarios;
    g_scenarios = scenarios;
    global g_scenarios_matrix;
    g_scenarios_matrix = scenarios_matrix;
    global g_start_kcats;

    origReactionNames = convertCharsToStrings(fieldnames(g_reactions_to_change));
    
    % preparation:
    % for every original reaction, collect stoichiometric coefficients of
    % parallel reactions
    for i = 1:numel(origReactionNames)
 
        origReactionName = origReactionNames(i);
        
        % get all coupled reactions for the reaction candidates
        splitReactions = g_reactions_to_change.(origReactionName);
        
        % for every split reaction, collect stoichiometric coefficients
        for splitReactionIndex = 1:numel(splitReactions)
                        
            splitReactionName = splitReactions(splitReactionIndex);
            reaction_index = strmatch(splitReactionName, cna_model.reacID, 'exact');
            pseudo_metabolite_indices = get_all_pseudo_metabolite_indices(splitReactionName, cna_model.stoichMat, cna_model.reacID, cna_model.specID);
            
            % Stores stoichiometric coefficients for every split reaction,
            % belonging to one original reaction
            g_start_kcats.(origReactionName).(splitReactionName) = zeros(size(pseudo_metabolite_indices));

            for j = 1:numel(pseudo_metabolite_indices)

                metabolite_index = pseudo_metabolite_indices(j);
                
                g_start_kcats.(origReactionName).(splitReactionName)(j) = cna_model.stoichMat(metabolite_index, reaction_index);

            end
        end
    end

    % every original reaction receives one kcat factor representing the
    % kcat.
    start_kcat_factors = ones(size(origReactionNames));
    lower_bounds = start_kcat_factors / max_change_factor;
    upper_bounds = start_kcat_factors * max_change_factor;

    % Run optimization :D
    options = optimoptions('fmincon','Algorithm','sqp');
       
    best_kcat_factors = fmincon(@p_get_objective_value_for_kcats, start_kcat_factors * 10, [], [], [], [], lower_bounds, upper_bounds, []);
    
    % Delete global variables
    clear g_cna_model;
    clear g_reactions_to_change;
    clear g_scenarios;
    clear g_scenarios_matrix;
end
