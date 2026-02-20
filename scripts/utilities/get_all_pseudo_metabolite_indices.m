function [pseudo_met_indices] = get_all_pseudo_metabolite_indices(reactionID ,stoichMatrix, reactions, metabolites)
% Returns the indices of all psudometabolites in the provided reaction

% Arguments:
% reactionID ~      reaction identifier of reaction to look up
% stoichMatrix ~    stoichiometric matrix from cna model
% reactions ~       vector of reaction ids from cna model
% metabolites ~     vector of metabolite ids from cna model

% Outputs:
% pseudo_met_indices ~  Vector of pseudo met indices

    reactionIndex = find(strcmp(cellstr(reactionID), reactions) == 1);
    
    metabolite_column = stoichMatrix(:, reactionIndex);
    metabolite_column(metabolite_column ~= 0) = 1;
    metabolite_indices = metabolite_column == 1;
    metaboliteIDs = cellstr(metabolites(metabolite_indices, :));
    
    pseudo_mets = metaboliteIDs(startsWith(metaboliteIDs, "M_ENZYME_"));
    pseudo_mets = [pseudo_mets; "M_prot_pool"];
    
    pseudo_met_indices = zeros(size(pseudo_mets));
    
    for i = 1:numel(pseudo_met_indices)
        pseudo_met_indices(i) = find(strcmp(pseudo_mets(i, :), metabolites) == 1);
    end 
        
end

