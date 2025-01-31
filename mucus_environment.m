function [updated_dish,butConc,human_proteo,bac_proteo] = mucus_environment(old_dish, human_proteo, bac_proteo, tot, butConc, mucBox) 

checking = [ 0, 9, 9.5] ; 
        
Refresh = [ 9.5 0 9.5 ; %replacement code for proteolyzed(?) mucus 
            0 0 0 ;
            9.5 0 9.5] ; 

updated_dish = old_dish ; 

for count1 = 2:1:tot-1
    for count2 = 2:1:tot-1
        mucRepl = Refresh ;
        if old_dish(count1, count2) == 9 && ...
                human_proteo(count1,count2) >= 75 && ...
                bac_proteo(count1,count2) >= 75 
            rF = count1-1 ; 
            rL = count1+1 ; 
            cF = count2-1 ; 
            cL = count2+1 ; 
            mucCount = 1 ; 
            for count4 = cF:1:cL % columns
                for count3 = rF:1:rL %rows
                    if ~ismember(old_dish(count3,count4),checking) % only want to change mucus/open places
                        mucRepl(mucCount) = old_dish(count3,count4) ;
                    end 
                    mucCount = mucCount+1 ; 
                end 
            end 
            human_proteo(count1,count2) = 0 ; 
            bac_proteo(count1,count2) = 0 ; 
            butConc(count1,count2) = butConc(count1,count2) + mucBox ; 
            updated_dish(rF:rL, cF:cL) = mucRepl ;
          
        end
    end
end            