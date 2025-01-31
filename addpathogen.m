function [agent,dish,dish_labels,bt3,bt1,bt2,overwrite] = addpathogen(seed,tot,dish, dish_labels,agent,paramsBT3,bt1,bt2)
[x,y] = find(dish==0) ; 
numstart = (length(agent))+1; 
dummycount = 1 ; 
bt3 = 0 ; 

    cols = ceil(seed/tot) ; 
    remend = rem(seed,tot)+1 ; 
    count = numstart ; 

    for count2 = tot:-1:tot-cols
        for count1 = 1:1:tot
            if count2 == tot-cols && count1 == remend
                break
            end 
      
            overwrite = dish_labels(count1,count2) ; 
            if overwrite ~= 0
                agent{overwrite}.alive = 19 ;
                agent{overwrite}.position = [ ] ;

                if isa(agent{overwrite}, 'BT1')
                    bt1 = bt1 - 1;
                end
                if isa(agent{overwrite}, 'BT3')
                    bt3 = bt3 - 1;
                end

                if isa(agent{overwrite}, 'BT2')
                    bt2 = bt2 - 1 ;
                end

                pos = [count1,count2] ;
                [dish,dish_labels] = Adjust(count, dish,dish_labels,3,pos(1),pos(2)) ;
                last_rep = ceil(rand*paramsBT3.rep_rate) ;
                food = ceil(rand*2) ;
                agent{count} = BT3(paramsBT3.rep_rate, last_rep,food, pos,1) ;
                bt3 = bt3 + 1 ;
                count = count + 1 ; % next agent number
            end 
        end 
     end 
% end 
end
