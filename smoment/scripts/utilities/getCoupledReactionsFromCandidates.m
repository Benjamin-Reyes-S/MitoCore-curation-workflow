function [coupledReactions] = getCoupledReactionsFromCandidates(canidateReactions, model)
% AutoPACMEN splits reactions with isoenzymes into "parallel" reactions.
% For each reaction for calibration, this function determines these
% parallel reactions.

% Argumnts:
% canidateReactions ~   string array of reaction names for calibration
% model ~               CNA model

    coupledReactions = struct;
    
    for rIndex = 1 : numel(canidateReactions)
        
        % reactionName = strcat("R_", canidateReactions(rIndex));
        reactionName = canidateReactions(rIndex);
        direction = "";
        idSplit = split(reactionName, '_TG_');

        % if it is a reversible reaction, get the direction suffix
        if numel(idSplit) == 2
            direction = strcat("_", idSplit(2));
        end

        originalReactionId = split(idSplit(1),'_GPRSPLIT_');
        originalReactionId = originalReactionId(1);
        
        splitReactions = getSplitReactions(reactionName, model);
        
        coupledReactions.(strcat(originalReactionId, direction)) = splitReactions;
        
    end
end

