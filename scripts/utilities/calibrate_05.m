function [outputStruct] = calibrate41(modelPath, candidateReactions, scenarioStruct, scenariosMatrix, changeFactor)
    % performs calibration for a model and one scenario. It optimizes kcats
    % first and then determines a protein pool by linesearch.
    
    % Arguments:
    % modelPath ~           path to an unclibrated smoment model
    % candidateReactions ~  string array encoding the reactions whose kcats should
    %                       be optimized.
    % scenarioStruct ~      autoPACMEN scenarios
    % scenariosMatrix ~     scenario dependency matrix
    % changeFactor ~        maximal factor the kcats are allowed to be changed    

    %% read model
    [modelUncalibrated, ~] = CNAsbmlModel2MFNetwork(modelPath);

    % collect calibrated values
    outputStruct = struct([]);
    
    % autoPacmen splits reactions with OR rule in their GPR. Candidate
    % reactions can be one of these split reactions. However, the kcat
    % applies to the complete reaction. Therefore, the kcat should be
    % changed for all splits. Coupled reaction is a nested struct containing all
    % splits in one direction for an original reaction
    % (struct.originalReactionName.splitReactionName).
    
    coupledReactions = getCoupledReactionsFromCandidates(candidateReactions, modelUncalibrated);

    %% Step 1: calibrate individual kcats
    poolReactionIndex = getIndexFromChar(modelUncalibrated.reacID, 'R_ER_pool_TG_');
    globalPoolLimit = modelUncalibrated.reacMax(poolReactionIndex);
    
%     modelUncalibrated.reacMax(poolReactionIndex) = 1000;
        
    originalReactionNames = convertCharsToStrings(fieldnames(coupledReactions));
                
    % calibrate kcats
    [best_kcat_factors, start_kcat_factors] = smoment_kcat_optimization(modelUncalibrated, coupledReactions, scenarioStruct, scenariosMatrix, changeFactor);
    
    global currentModelName;
        
    disp('current model:' + currentModelName);
    
    % update all candidate reactions
    for i = 1:numel(originalReactionNames)
        
        % get all split reactions for original reaction
        originalReactionName = originalReactionNames(i);
        splitReactions = coupledReactions.(originalReactionName);
        
        % for every split reaction write kcats to output dict and update in
        % model
        for j = 1:numel(splitReactions)
           
            splitReactionName = splitReactions(j);
            
            disp('Applying calibration to ' + splitReactionName);
            
            % write factor to output struct
            outputStruct(1).reactions.(splitReactionName) = best_kcat_factors(i);

            reactionIndex = getIndexFromChar(modelUncalibrated.reacID, splitReactionName);

            pseudoMetaboliteIndices = get_all_pseudo_metabolite_indices(splitReactionName, modelUncalibrated.stoichMat, modelUncalibrated.reacID, modelUncalibrated.specID);

            % for every psudo metabolite, apply kcat factor
            disp('new stoichiometries:');
            
            for k = 1:numel(pseudoMetaboliteIndices)

                metaboliteIndex = pseudoMetaboliteIndices(k);
                metaboliteName = modelUncalibrated.specID(metaboliteIndex,:);
                                
                stoichiometry = modelUncalibrated.stoichMat(metaboliteIndex, reactionIndex);
                newStoichiometry = stoichiometry * 1/best_kcat_factors(i);
                
                disp(string(metaboliteName) + ' ' + newStoichiometry);

                modelUncalibrated.stoichMat(metaboliteIndex, reactionIndex) = newStoichiometry;
            end
        end
    end
    
    [~, success, ~, maximization_result] = CNAoptimizeFlux(modelUncalibrated, [], [], 0, -1);
    disp('FBA result before prot pool calibration:')
    disp(-maximization_result)
    
    %% Step 2: optimize global protein pool by linesearch
%     modelUncalibrated.reacMax(poolReactionIndex) = globalPoolLimit;
    
    [~, success, ~, maximization_result] = CNAoptimizeFlux(modelUncalibrated, [], [], 0, -1);
    
    protPool = [0.0:0.01:1];
    outputs = zeros(size(protPool));
    
    global g_cna_model;
    g_cna_model = modelUncalibrated;
    global g_scenarios;
    g_scenarios = scenarioStruct;
    global g_scenarios_matrix;
    g_scenarios_matrix = scenariosMatrix;
    
    for i=1:length(protPool)

        outputs(i) = p_get_objective_value_for_prot_pool(protPool(i));

    end

    minimumVal = min(outputs(2:length(outputs)));
    minimumIndex = outputs(2:length(outputs)) == minimumVal;
    prot_pool_skip_first = protPool(2:length(outputs));
    bestPool = prot_pool_skip_first(minimumIndex);
    bestPool = bestPool(1);
    
    outputStruct(1).bestPool = bestPool;
    
    fprintf('best prot pool: %.3f \n', bestPool);
  
    % Note:
    % Calculating the objective value is based on setting the
    % stoichiometry of the prot_pool pseudometabolite in its delivery
    % reaction, which is 1 by default. However, this delivery is actually
    % restricted by its upper constraint. As the stoichiometry just
    % scales the "efficiency of global protein production", its also
    % possible to scale the upper constraint with the new stoichiometry
    % to achieve the same limitiation of the global protein pool.
    % Effectively, this makes no difference. The stoichiometry could
    % be seen as the utilized fraction of the measured protein content,
    % making the model slightly more informative than just updating the
    % upper constraint.
    % modelUncalibrated.reacMax(poolReactionIndex) = bestPool * globalPoolLimit;
      
    poolMetaboliteIndex = getIndexFromChar(modelUncalibrated.specID, 'M_prot_pool');
    modelUncalibrated.stoichMat(poolMetaboliteIndex, poolReactionIndex) = bestPool;
    
    [~, success, ~, maximization_result] = CNAoptimizeFlux(modelUncalibrated, [], [], 0, -1);
    
    outputStruct(1).newObjectiveVal = -maximization_result;

    disp('Final FBA result: ')
    disp(-maximization_result)
end

