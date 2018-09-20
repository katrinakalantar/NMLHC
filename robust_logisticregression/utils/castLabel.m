% Description 
%        Function to cast between different type of 
%        label representation i.e. from {-1,1} to {0,1} to {1,2}
% Author
%        Jakramate Bootkrajang
% Input  
%        y = label vector
%        t = target representation choose from {-1,0,2}
% Output  
%        y = casted label vector
%==========================================================================

function y = castLabel(y, t)

if (length(y) == 1)
    error('All value of y required to recognise current format');
end

if (contains(y,-1))
    % {-1,1} input 
    if (t == -1)
        y = y; % do nothing, included for clarity
    elseif (t == 0)
        y = (y + 1) ./ 2;
    elseif (t == 2)
        y = (y + 3) ./ 2;        
    end
elseif (contains(y,0))
    % {0,1} input
    if (t == -1)
        y = y .* 2 - 1; 
    elseif (t == 0)
        y = y; % do nothing, included for clarity
    elseif (t == 2)        
        y = y + 1;
    end
elseif (contains(y,2))
    % {1,2} input    
    if (t == -1)
        y = y .* 2 - 3;
    elseif (t == 0)
        y = y - 1;
    elseif (t == 2)
        y = y; % do nothing, included for clarity
    end
end
    
            
    