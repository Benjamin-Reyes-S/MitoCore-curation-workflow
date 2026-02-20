function [result] = getIndexFromChar(charArray, query)
    result = find(matches(strip(string(charArray)), query));
end

