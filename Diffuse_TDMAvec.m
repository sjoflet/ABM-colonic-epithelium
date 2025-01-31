function [C_new] = Diffuse_TDMAvec(C_0, Dcoeff, total, t_fin, del, BCtypes,BCalpha) 

% ** Input Information **
% 
% C_0 = initial concentrations for b vector
% Dcoeff = diffusion coefficient matrix
% total = size of matrix
% del_x = delta x
% del_y = delta y
% t_fin = desired final time in seconds (assuming t_0 to be 0)
% del_t = delta t 
% BCs = matrix with Boundary Condition information
    % 1 = Epithelial
    % 2 = Top
    % 3 = Lumen
    % 4 = Bottom

% -----------------------------------------------------------------------
% n_r = (total - 1) ; % defines n-step by number of points (n0 is the initial condition and thus not included)
% n_c = (total - 1); %                       "" 

C_1 = zeros(total,total) ; % matrix that will update with new values every iteration

% * Alternating Direction Implicit Method Loop 
% del_t = del(1) ; 
% del_r = del(2) ; 
% del_c = del(3) ; 

dt = del(1) ;

for t = dt:dt:t_fin  % iterations in time

    % -----------------------
    % ** Boundary Conditions **
       % 4x2 matrix with first column true/false of Dirichlet and second as
       % values 

       %   Row 1 = Boundary 1 ( epithelial, j = 1 ) % No flux - tight junctions  
       %   Row 2 = Boundary 2 ( top, i = 1 ) % No Flux - symmetry
       %   Row 3 = Boundary 3 ( lumen, j = total ) % constant Concentration
       %   Row 4 = Boundary 4 ( bottom, i = total ) % No flux - symmetry
    
   [u,d,l,bvec] = resetTDMA(total) ; 
   % look at boundaries for r-sweep

   Bnum = 1 ; % Epithelial  
   [u,d,l,bvec] = BoundaryCalcs(Bnum, Dcoeff, total, del, C_0, u,d,l, bvec, BCtypes,BCalpha) ;
   C_1(:,1) = TDMA(u,d,l,bvec) ;

   [u,d,l,bvec] = resetTDMA(total) ; 

   Bnum = 3 ; % Lumen
   [u,d,l,bvec] = BoundaryCalcs(Bnum, Dcoeff, total, del, C_0, u,d,l, bvec, BCtypes,BCalpha) ; 
   C_1(:,total) = TDMA(u,d,l,bvec) ;

   % % inner r-sweep 
   for countc = 2:1:total-1 % r-sweep without c=1 and c = total
       
      [u,d,l,bvec] = resetTDMA(total) ; % reset

      [u(1),d(1),l(1),bvec(1)] = Edges(1, countc, total, Dcoeff(1,countc), C_0, del, BCtypes) ; 
      [u(total),d(total),l(total),bvec(total)] = Edges(total, countc, total, Dcoeff(total,countc), C_0, del, BCtypes) ; 

      d1 = Dcoeff(2:total-1,countc)/2*del(1)/(del(2)^2) ; 
      d2 = Dcoeff(2:total-1,countc)/2*del(1)/(del(3)^2) ; 
      a = - d1 ; 
      b = (1+2*d1) ; 

      u(2:total-1) = a ;
      d(2:total-1) = b ; 
      l(2:total-1) = a ; 
      bvec(2:total-1) = d2.*C_0(2:total-1,countc-1) + (1-2*d2).*C_0(2:total-1,countc) + d2.*C_0(2:total-1,countc+1) ;

      C_1(:,countc) = TDMA(u,d,l,bvec) ; % solve for new values after going through every r for that c
   end

   C_0 = C_1 ; % update to n+1/2
   
 
   %look at boundaries for c -sweep - have to solve after r-sweep to
   %obtain values n +1/2 for r=2 and r = total-1
   
   %reset
   [u,d,l,bvec] = resetTDMA(total) ; 
   
   Bnum = 2 ; % top
   [u,d,l,bvec] = BoundaryCalcs(Bnum, Dcoeff, total, del, C_0, u,d,l, bvec, BCtypes,BCalpha) ; 
   C_1(1,:) = TDMA(u,d,l,bvec) ; 
    
   %reset
   [u,d,l,bvec] = resetTDMA(total) ; 

   Bnum = 4 ; % bottom
   [u,d,l,bvec] = BoundaryCalcs(Bnum, Dcoeff, total, del, C_0, u,d,l, bvec, BCtypes,BCalpha) ; 
   C_1(total,:) = TDMA(u,d,l,bvec) ;   % updated to C n+1/2 values
  
    for countr = 2:1:total-1 % c-sweep without r = 1, r = total
       %reset
       [u,d,l,bvec] = resetTDMA(total) ; 

       [u(1),d(1),l(1),bvec(1)] = Edges(countr, 1, total, Dcoeff(countr,1), C_0, del, BCtypes) ;
       [u(total),d(total),l(total),bvec(total)] = Edges(countr, total, total, Dcoeff(countr,total), C_0, del, BCtypes) ;

       d1 = Dcoeff(countr,2:total-1)/2*del(1)/(del(2)^2) ;
       d2 = Dcoeff(countr,2:total-1)/2*del(1)/(del(3)^2) ;
       a = - d2 ;
       b = (1+2*d2) ;

       u(2:total-1) = a ;
       d(2:total-1) = b ;
       l(2:total-1) = a ;
       bvec(2:total-1) = d1.*C_0(countr-1,(2:total-1)) + (1-2*d1).*C_0(countr,(2:total-1)) + d1.*C_0(countr+1,(2:total-1)) ;

       C_1(countr,:) = TDMA(u,d,l,bvec);  % solve for new values after going through every r for that c
    end
C_0 = C_1 ; % update to n+1
end  

C_new = C_1 ; 
end 

function [u,d,l,bvec] = BoundaryCalcs(Bnum, D, total, del, C_0, u,d,l, bvec, BCtypes,BCalpha) 

    if Bnum == 1 
        c = 1 ; 
% BC1: Epithelial Layer (:,1) - all rows at c = 1 
    % Dirichlet     
    % C_1(:,1) = C_1(:,1) ; 
        if BCtypes(Bnum) == 1 % Dirichlet
            d(:,1) = 1 ; % identity vector
            bvec(:,1) = C_0(:,c) ; % values unchanged 
            return
        else 
            % Neumann 

            d2temp = D(2:total-1,c)/2*del(1)/(del(3)^2) ;
            atemp = - d2temp ;
            btemp = (1+2*d2temp) ;

            u(2:total-1) = atemp ;
            d(2:total-1) = btemp ;
            l(2:total-1) = atemp ;
            bvec(2:total-1) = 2*d2temp.*C_0(2:total-1, c+1)   + (1-2*d2temp).*C_0(2:total-1,c) - 2*d2temp*BCalpha(1)*del(3) ;
           
            clear d2temp atemp btemp
           
            % -------------------------------------------------------------  
            r = 1 ;  
            if BCtypes(2) == 1 
                d(r) = 1 ; % identity
                bvec(r) = C_0(r,c) ; % unchanged
            else
                [a, b, ~, d2] = abcdDef(del, D(r,c), 1) ;
                d(r) = b ; 
                u(r) = 2*a ; 
                bvec(r) = 2*d2*C_0(r, c+1) + (1-2*d2)*C_0(r,c) - 2*d2*BCalpha(1)*del(3) + 2*a*BCalpha(2)*del(2) ;
            end

            r = total ; 
            if BCtypes(4) == 1 
                d(r) = 1 ; 
                bvec(r) = C_0(r,c) ; 
            else 
                [a, b, ~, d2] = abcdDef(del, D(r,c),1) ;
                d(r) = b ; 
                l(r) = 2*a ; 
                bvec(r) = 2*d2*C_0(r, c+1) + (1-2*d2)*C_0(r,c) - 2*d2*BCalpha(1)*del(3) - 2*a*BCalpha(2)*del(2) ;
            end

        end 
        return % don't need to look at other Bnums 
    end 

    % BC2: Top (1, :) - all columns at r = 1 
    % Dirichlet

    if Bnum == 2 % Top side: r = 1 
        r = 1 ; 

        if BCtypes(Bnum) == 1 % Dirichlet
            d(:,1) = 1 ; % identity vector
            bvec(:,1) = C_0(r,:) ; % values unchanged 
            return
        else 
            % Neummann
            % b(j) = 2*d1*C_1(i+1,j) + (1-2*d1)*C(i,j) - 2*alpha*d1*delx

            d1temp = D(r,2:total-1)/2*del(1)/(del(2)^2) ;
            atemp = - d1temp ;
            btemp = (1+2*d1temp) ;

            u(2:total-1) = atemp ;
            d(2:total-1) = btemp ;
            l(2:total-1) = atemp ;
            bvec(2:total-1) = 2*d1temp.*C_0(r+1, 2:total-1)   + (1-2*d1temp).*C_0(r,2:total-1) - 2*d1temp*BCalpha(2)*del(2) ;
           
            clear d1temp atemp btemp
                
            
    % -------------------------------------------------------------  
            c = 1 ;  
            if BCtypes(1) == 1 
                d(c) = 1 ; 
                bvec(c) = C_0(r,c) ; 
            else 
                [a, b, d1, ~] = abcdDef(del, D(r,c),2) ;
                d(c) = b ; 
                u(c) = 2*a ; 
                bvec(c) = 2*d1*C_0(r+1, c)   + (1-2*d1)*C_0(r,c) - 2*d1*BCalpha(2)*del(2) + 2*a*BCalpha(1)*del(3) ;
            end

            c = total ; 
            if BCtypes(3) == 1 
                d(c) = 1 ; 
                bvec(c) = C_0(r,c) ; 
            else 
                [a, b, d1, ~] = abcdDef(del, D(r,c),2) ;
                d(c) = b ; 
                l(c) = 2*a ; 
                bvec(c) = 2*d1*C_0(r+1,c)   + (1-2*d1)*C_0(r,c) - 2*d1*BCalpha(2)*del(2) - 2*a*BCalpha(3)*del(3) ;
            end

        end 
        
        return % don't need to look at other Bnums 
    end 

% BC3: Lumen (:,total) - all rows at j = total
    % Dirichlet
    % C_1(:,total) = C_1(:,total) 
    % continue

    % Neumann
    % b(:) = 2*d2*C_1(i,j-1) + (1-2*d2)*C_1(i,j) + 2*d2*alpha*dely 

    % b(1) = b(1) + 2*alpha*a1*delx
    % A(1, 1:2) = [b1 (a1+c1) ]
    %
    % b(total) = b(total) - 2*alpha*c1*delx
    % A(total, total-1:total) = [ (a1+c1) b1] 
    if Bnum == 3 
        c = total ; 
        if BCtypes(Bnum) == 1 % Dirichlet
            d(:,1) = 1 ; % identity vector
            bvec(:,1) = C_0(:,c) ; % values unchanged 
            return
        else 
            
            d2temp = D(2:total-1,c)/2*del(1)/(del(3)^2) ;
            atemp = - d2temp ;
            btemp = (1+2*d2temp) ;

            u(2:total-1) = atemp ;
            d(2:total-1) = btemp ;
            l(2:total-1) = atemp ;
            bvec(2:total-1) = 2*d2temp.*C_0(2:total-1, c-1)   + (1-2*d2temp).*C_0(2:total-1,c) + 2*d2temp*BCalpha(3)*del(3) ;
           
            clear d2temp atemp btemp 
            
    % -------------------------------------------------------------  
            r = 1 ;  
            if BCtypes(2) == 1 
                d(r) = 1 ; 
                bvec(r) = C_0(r,c) ; 
            else 
                [a, b, ~, d2] = abcdDef(del, D(r,c), 1) ;
                d(r) = b ; 
                u(r) = 2*a ; 
                bvec(r) = 2*d2*C_0(r, c-1) + (1-2*d2)*C_0(r,c) - 2*d2*BCalpha(3)*del(3) + 2*a*BCalpha(2)*del(2) ;
            end

            r = total ; 
            if BCtypes(4) == 1 
                d(r) = 1 ; 
                bvec(r) = C_0(r,c) ; 
            else 
                [a, b, ~, d2] = abcdDef(del, D(r,c),1) ;
                d(r) = b ; 
                l(r) = 2*a ; 
                bvec(r) = 2*d2*C_0(r, c-1) + (1-2*d2)*C_0(r,c) - 2*d2*BCalpha(3)*del(3) - 2*a*BCalpha(4)*del(2) ;
            end

        end 
        return % don't need to look at other Bnums 
    end 

% BC4:  Bottom (total, :) -> all j for i = total

    % Dirichlet
    % C_1(total,:) = C_1(total,:) 
    % continue

    % Neumann 
    % b(:) = 2*d1*C_1(i+1,j) + (1-2*d1)*C_1(i,j) + 2*a1*alpha*delx

    % b(1) = b(1) + 2*alpha*a2*dely
    % A(1,1:2) = [b2 (a2+c2)]
    % 
    % b(total) = b(total) - 2*alpha*c2*dely
    % A(total, total-1:total] = [(a2+c2) b2]

    if Bnum == 4 % Bottom side: r = total 
        r = total ; 

        if BCtypes(Bnum)== 1 % Dirichlet
            d(:,1) = 1 ; % identity vector
            bvec(:,1) = C_0(r,:) ; % values unchanged 
            return
        else 
    
            d1temp = D(r,2:total-1)/2*del(1)/(del(2)^2) ;
            atemp = - d1temp ;
            btemp = (1+2*d1temp) ;

            u(2:total-1) = atemp ;
            d(2:total-1) = btemp ;
            l(2:total-1) = atemp ;
            bvec(2:total-1) = 2*d1temp.*C_0(r-1,2:total-1)   + (1-2*d1temp).*C_0(r,2:total-1) + 2*d1temp*BCalpha(4)*del(2) ;
           
            clear d1temp atemp btemp

    % -------------------------------------------------------------        
            c = 1 ;  
            if BCtypes(1) == 1 
                d(c) = 1 ; 
                bvec(c) = C_0(1,c) ; 
            else 
                [a, b, d1, ~] = abcdDef(del, D(r,c),2) ; 
                d(c) = b ; 
                u(c) = 2*a ; 
                bvec(c) = 2*d1*C_0(r-1, c) + (1-2*d1)*C_0(r,c) - 2*d1*BCalpha(4)*del(2) + 2*a*BCalpha(1)*del(3) ;
            end

            c = total ; 
            if BCtypes(3) == 1 
                d(c) = 1 ; 
                bvec(c) = C_0(1,c) ; 
            else 
                [a, b, d1, ~] = abcdDef(del, D(r,c),2) ;
                d(c) = b ; 
                l(c) = 2*a ; 
                bvec(c) = 2*d1*C_0(r-1,c) + (1-2*d1)*C_0(r,c) - 2*d1*BCalpha(4)*del(2) - 2*a*BCalpha(3)*del(3) ; 
            end

        end 
        return % don't need to look at other Bnums 
    end 
 
end 

function [u,d,l,bvec] = Edges(r,c,total,D,C_0,del, BCtypes) 
u = 0 ;
d = 0 ; 
l = 0 ; 
bvec = 0 ; 
   if c == 1 % Epithelial layer, c-sweep edges for r= 2 to r= total-1 (r=1 and r =total are already solved
       if BCtypes(1) == 1 %Dirichlet
           d = 1 ; 
           bvec = C_0(r,c) ; 
       else 
           [a, b, d1, ~] = abcdDef(del, D, 2); 
           bvec = d1*C_0(r-1,c) + (1-2*d1)*C_0(r,c) + d1*C_0(r+1,c)  ;
           d = b ; 
           u = 2*a ; 
       end 
   end

   if c == total % Lumen c-sweep 
       if BCtypes(3) == 1 %Dirichlet
           d = 1 ; 
           bvec = C_0(r,c) ; 
       else 
           [a, b, d1, ~] = abcdDef(del, D, 2); 
           bvec = d1*C_0(r-1,c) + (1-2*d1)*C_0(r,c) + d1*C_0(r+1,c) ;
           d = b ; 
           l = 2*a ; 
       end 
   end

   if r == 1 % top, r -sweep 
       if BCtypes(2) == 1 % Dirichlet
           d = 1 ; 
           bvec = C_0(r,c) ; 
       else % Neumann
           [a, b, ~, d2] = abcdDef(del, D, 1);
           bvec = d2*C_0(r,c-1) + (1-2*d2)*C_0(r,c) + d2*C_0(r,c+1)  ;
           d = b ; 
           u = 2*a ; 
       end
   end 
   
   if r == total % bottom, r -sweep 
       if BCtypes(4) == 1 % Dirichlet
           d = 1 ; 
           bvec = C_0(r,c) ; 
       else 
           [a, b, ~, d2] = abcdDef(del, D, 1);
           bvec = d2*C_0(r,c-1) + (1-2*d2)*C_0(r,c) + d2*C_0(r,c+1)  ;
           d = b ; 
           l = 2*a ; 
       end
   end 

end

function [a,b, d1, d2] = abcdDef(del, D, which) 
    if which == 1 % want a1 b1
        d1 = (D/2)  * del(1)/(del(2)^2) ; 
        d2 = (D/2)  * del(1)/(del(3)^2) ; 
        a = -d1 ; 
        b = (1+2*d1) ;
    elseif which == 2 % want a2,b2
        d1 = (D/2)  * del(1)/(del(2)^2) ; 
        d2 = (D/2)  * del(1)/(del(3)^2) ; 
        a = -d2 ; 
        b = (1+2*d2) ;
    end 
end

function [Cnew] = TDMA(u,d,l,b)  
    total = length(b) ;  
    lambda = zeros(total,1) ; 
    rho = zeros(total,1) ; 
    X = zeros(total,1) ; 

    % For n = 1 
    lambda(1) = u(1)/d(1) ; 
    rho(1) = b(1)/d(1) ; 

    % For all other n 
    for n = 2:1:total
        lambda(n) = u(n) / (d(n)-l(n)*lambda(n-1)) ; 
        rho(n) = (b(n)-l(n)*rho(n-1)) / (d(n)-l(n)*lambda(n-1)) ; 
    end 

    % For n = total
    X(total) = rho(total) ; 

    % For all other n 
    for n = total-1:-1:1
        X(n) = rho(n) - lambda(n)*X(n+1) ; 
    end

    Cnew = X ; 
end

function [u,d,l,b] = resetTDMA(total) % function to clear TDMA column vectors easily based on dish size
u = zeros(total,1) ; % upper diagonal
d = zeros(total,1) ; % diagonal
l = zeros(total,1) ; % lower diagonal
b = zeros(total,1) ; % b vector
end 