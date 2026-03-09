function [f] = p_get_penalty_for_mesured_reactions(model, flux_vector, scenario_struct)
    
    total_penalty = 0;

    if isfield(scenario,'secondary_targets')
        fildnames = fieldnames(scenario.secondary_targets);
        for secondary_target_index = 1:numel(fildnames)
            fieldname = string(fildnames(secondary_target_index));
            target_obj = scenario.secondary_targets.(fieldname);

            reaction_fields = fieldnames(target_obj.reactions);
            
            target_penalty = 0;
            
            % in case of multiple reactions, calculate linear combination
            % of fluxes and factors
            for reaction_index = 1:numel(reaction_fields)
                reaction_name = string(reaction_fields(reaction_index));
                factor = target_obj.reactions.(reaction_name);
                
                reaction_index = strmatch(reaction_name, model.reacID, 'exact');
                reaction_flux = flux_vector(reaction_index) * factor;
                
                target_penalty = target_penalty + reaction_flux;
            end
            
            % determine penalty
            target_penalty = abs(target_penalty/target_obj.value);
            if (target_penalty < 1.0)
                target_penalty = 1 / target_penalty;
            end
            target_penalty = target_penalty - 1;
            
            total_penalty = total_penalty + target_penalty;
        end
        f = total_penalty * 0.8; % multiply by weight for measurements
    else
        f = 0;
    end

end
