function [check, newcoors] = Cell_neighbors(x1, x2, X, plate, tot, outputamt)
%% INPUTS

% x1, x2 = center position 
% X  = number you are seeking 
% plate = entire dish  
% tot = size of dish 
% outputamt = max number of coordinates you want outputed 

poss_coors = [ (x1+1) x2 ; ...% 1 = [(cpos(1)+1) cpos(2)] ; % look right
    (x1-1) x2 ; ... % 2 = [(cpos(1)-1) cpos(2)] ; % look left
    x1 (x2+1) ; ... % 3 = [cpos(1) (cpos(2)+1)] ; % look down
    x1 (x2-1) ; ... % 4 = [cpos(1) (cpos(2)-1)] ; % look up
    (x1-1) (x2+1) ; ... % 5 = [(cpos(1)-1) (cpos(2)+1)] % diag down left
    (x1+1) (x2+1) ; ... % 6 = [(cpos(1)+1) (cpos(2)+1)] % diag down right
    (x1+1) (x2-1) ; ... % 7 = [(cpos(1)+1) (cpos(2)-1)] % diag up right
    (x1-1) (x2-1) ] ; % 8 = [(cpos(1)-1) (cpos(2)-1)] % diag up left

  %[ agent{agnts+agnts_new}, agnts_new ] = neighbors(cpos(1),cpos(2), dish, curr) ; 
  % 

choices = [ ] ; 

 if x1 ~= tot && plate((x1+1), x2) == X % look right
     choices = [choices 1] ; 
 end 

 if x1 ~= 1 && plate((x1-1), x2) == X % look left
     choices = [choices 2] ; 
 end 

 if x2 ~= tot && plate(x1, (x2+1)) == X % look down
    choices = [choices 3] ; 
 end 

 if x2 ~= 1 && plate(x1, (x2-1)) == X % look up
    choices = [choices 4] ; 
 end 

 if x1 ~= 1 && x2 ~= tot &&  plate((x1-1), (x2+1)) == X % look down left
     choices = [choices 5] ; 
 end 

 if x1 ~= tot && x2 ~= tot &&  plate((x1+1), (x2+1)) == X % look down right
     choices = [choices 6] ; 
 end 

 if x1 ~= tot && x2 ~= 1 &&  plate((x1+1), (x2-1)) == X % look up right
     choices = [choices 7] ; 
 end 

 if x1 ~= 1 && x2 ~= 1 &&  plate((x1-1), (x2-1)) == X % look up left
     choices = [choices 8] ; 
 end 

 blankcheckmatrix = [ ]  ; 
 check = isequal(choices, blankcheckmatrix) ; % negative check to make sure choices matrix is not still empty
 
 if check == 0 % False if choices is NOT empty
     if size(choices,2) == 1  % randsample doesn't work with one scalar input
         choice = choices ; 
         newcoors = poss_coors(choice, :)  ;
     elseif size(choices,2) >= outputamt % want at max the desired output amount
         choice = randsample(choices,outputamt) ; % randsample returns a vector
         newcoors = poss_coors(choice, :)  ; 
     else
         number = size(choices,2) ; % if less than the output amount, want all of them
         choice = randsample(choices,number) ;
         newcoors = poss_coors(choice, :) ; 
     end 
 else
     newcoors = [ ]  ; 
 end 

end 