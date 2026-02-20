function [splitReactions] = getSplitReactions(reactionName, model)

    splitReactions = reactionName;
    
    direction = "";
    
    idSplit = split(reactionName, '_TG_');
    
    % if it is a reversible reaction, get the direction suffix
    if numel(idSplit) == 2
        direction = idSplit(2);
    end
    
    originalReactionId = split(idSplit(1),'_GPRSPLIT_');
    originalReactionId = originalReactionId(1);
    
    for reactionIndex = 1:numel(model.reacID(:,1))
       reactionId = strip(convertCharsToStrings(model.reacID(reactionIndex,:)));
       
       origReactionForCompare = split(reactionId,'_GPRSPLIT_');
       origReactionForCompare = origReactionForCompare(1);
       
       if strcmp(origReactionForCompare, originalReactionId) && endsWith(reactionId, direction)
           if ~any(contains(splitReactions, reactionId))
                splitReactions = [splitReactions; reactionId];
           end
       end 
    end
end

