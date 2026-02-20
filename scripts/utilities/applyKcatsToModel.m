function [alteredModel] = applyKcatsToModel(model, calbibratedKcatFactors)
    
    alteredModel = model;

    reactionNames = fieldnames(calbibratedKcatFactors);

    % apply kcat factors
    for i = 1:numel(reactionNames)
           
        reactionName = reactionNames(i);
        reactionName = string(reactionName);
        reaction_index = strmatch(reactionName, alteredModel.reacID, 'exact');

        new_kcat_factor = calbibratedKcatFactors.(reactionName);
        pseudo_metabolite_indices = get_all_pseudo_metabolite_indices(reactionName, alteredModel.stoichMat, alteredModel.reacID, alteredModel.specID);

        for k = 1:numel(pseudo_metabolite_indices)
            metabolite_index = pseudo_metabolite_indices(k);
            alteredModel.stoichMat(metabolite_index, reaction_index) = alteredModel.stoichMat(metabolite_index, reaction_index) * 1/new_kcat_factor;
        end
            
    end 

end

