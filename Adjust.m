function [dish, labels] = Adjust(agentNum, dish, labels, NewNum, x1, x2)
dish(x1,x2) = NewNum ; 
    if NewNum == 0 
        labels(x1,x2) = 0 ; 
    else 
        labels(x1,x2) = agentNum; 
    end 
end 