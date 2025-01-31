function [metabolites,totalmets] = MetaboliteStruct(sz,fiberfeed) 

% Useful Pre-sets

blank_field = ones(sz,sz) ; 

D_sugarsAQ = 6.7e2 ; % um^2/s
% D_sugarsBio = 1.675e2 ; % um^2/s
D_oxygen = 2e-2 ; % um^2/s
D_proteo = 1e-2  ; % um^2/s  
D_antibio = 1e-3  ; % um^2/s
D_fibers = D_sugarsAQ*5 ; % um^2/s
D_H2aq = 4.5e-5 ; % Aqueous H2 diffusion at T = 298 K
D_nitrateAQ = 1.7e3 ; % Aqueous NO3- diffision at 303 K
D_IgA = 1e3 ; % um^2/s
D_cytokines = 1e-3; % um^2/s

metabolites = struct('type',[],'alpha', [],'conc',[], 'Dcoeff',[],'Dbaseline', [], 'DdM', [], 'DlM', [] ) ; 

%% Metabolite Key
% 1 = Galactose 
% 2 = Galactitol
% 3 = Butyrate
% 4 = Bacteria Proteolytic Factors
% 5 = Human Proteolytic Factors
% 6 = Antibiotics
% 7 = Oxygen
% 8 = Antipathogens
% 9 = PAMPs
% 10 = Hydrogen
% 11 = Dietary Fiber
% 12 = Nitrate
% 13 = IgA
% 14 = MAMPs

totalmets = 14 ; 

%% Boundary Condition Key
    % 1 = Epithelial
    % 2 = Top of dish
    % 3 = Lumen
    % 4 = Bottom of dish 

    % 1 = Dirichlet: constant value
    % 2 = Neumann: Constant flux over boundaries

    % t = type
    % a = alpha

    % * Carbon Soources isolated within system - all Neumann BCs * 
    
%% 1: Galactose 

metabolites(1).type(1:4) = 2 ; % Neumann
metabolites(1).alpha(1:4) = 0 ; % No Flux

metabolites(1).conc = 10*blank_field ; % Initial Concentration Matrix

metabolites(1).Dcoeff = D_sugarsAQ*blank_field ; 
metabolites(1).trapDm = .25 ; 
metabolites(1).trapLm = .125 ; 
metabolites(1).sat = 10 ; 

metabolites(1).degradeMAX = 0 ; 
metabolites(1).degradeMIN = 0 ; 


%% 2: Galactitol 

metabolites(2).type(1:4) = 2 ; % Neumann
metabolites(2).alpha(1:4) = 0 ; % No Flux

metabolites(2).conc = 20*blank_field ; % Initial Concentration Matrix

metabolites(2).Dcoeff = D_sugarsAQ*blank_field ; 
metabolites(2).trapDm = .25 ; 
metabolites(2).trapLm = .125 ;
metabolites(2).sat = 10 ; 

metabolites(2).degradeMAX = 0 ; 
metabolites(2).degradeMIN = 0 ;

%% 3: Butyrate
metabolites(3).type(1:4) = 2 ; % Neumann
metabolites(3).alpha(1:4) = 0 ; % No Flux

metabolites(3).conc = .10*blank_field ; % Initial Concentration Matrix

metabolites(3).Dcoeff = D_sugarsAQ*blank_field ; 
metabolites(3).trapDm = .25 ; 
metabolites(3).trapLm = .125 ; 
metabolites(3).sat = 0 ; 

metabolites(3).degradeMAX = 0 ; 
metabolites(3).degradeMIN = 0 ; 

%% 4: Bacteria Proteolytic Factors

metabolites(4).type([1,3]) = 1 ; % Dirichlet (stay 0 at epithelial layer and always come from lumen)
metabolites(4).type([2,4]) = 2 ;  % Neumann
metabolites(4).alpha(1:4) = 0 ; % No Flux

metabolites(4).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(4).conc(:,sz) = 100 ; % Constant Lumen value

metabolites(4).Dcoeff = D_proteo*blank_field ; 
metabolites(4).trapDm = 0 ; 
metabolites(4).trapLm = 0 ; 
metabolites(4).sat = 0 ; 

metabolites(4).degradeMAX = 1 ;
metabolites(4).degradeMIN = .5 ;

%% 5: Human Proteolytic Factors  

metabolites(5).type(1) = 1 ; % Dirichlet constant epithelial release (value can change during overall iterations but constant during diffusion time)
metabolites(5).type(2:4) = 2 ;  % Neumann
metabolites(5).alpha(1:4) = 0 ; % No Flux

metabolites(5).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(5).conc(:,1) = 20 ; % Initial Epithelial value

metabolites(5).Dcoeff = D_proteo*blank_field ; 
metabolites(5).trapDm = 0 ; 
metabolites(5).trapLm = 0 ; 
metabolites(5).sat = 0 ; 

metabolites(5).degradeMAX = 1 ; 
metabolites(5).degradeMIN = .5 ; 

%% 6: Antibiotic - against E coli

metabolites(6).type([1,3]) = 1 ; % Dirichlet (constant release at epithelial layer and always come from lumen)
metabolites(6).type([2,4]) = 2 ;  % Neumann
metabolites(6).alpha(1:4) = 0 ; % No Flux

metabolites(6).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(6).conc(:,1) = 0.4e-12 ; % Starting Epithelial value

metabolites(6).Dcoeff = D_antibio*blank_field ; 
metabolites(6).trapDm = .5 ; 
metabolites(6).trapLm = .25 ;
metabolites(6).sat = 0 ; 

metabolites(6).degradeMAX = 0.75 ; 
metabolites(6).degradeMIN = 0.25 ; 

%% 7: Oxygen

metabolites(7).type(1) = 1 ; % Dirichlet constant epithelial release (value can change during overall iterations but constant during diffusion time)
metabolites(7).type(2:4) = 2 ;  % Neumann
metabolites(7).alpha(1:4) = 0 ; % No Flux

metabolites(7).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(7).conc(:,1:round(sz*(1/3))) = 32*(91)*10^(-15) ;  % ug/um^3

metabolites(7).Dcoeff = D_oxygen*blank_field ; 
metabolites(7).trapDm = .5 ; 
metabolites(7).trapLm = .25 ; 
metabolites(7).sat = 0 ; 

metabolites(7).degradeMAX = 0 ; 
metabolites(7).degradeMIN = 0 ; 

%% 8: Anti-Pathogen - against pathogen

metabolites(8).type([1,3]) = 1 ; % Dirichlet (constant release at epithelial layer and always come from lumen)
metabolites(8).type([2,4]) = 2 ;  % Neumann
metabolites(8).alpha(1:4) = 10 ; % No Flux

metabolites(8).conc = 0*blank_field ; % None to start, should be triggered by pathogens

metabolites(8).Dcoeff = D_antibio*blank_field ; 
metabolites(8).trapDm = .5 ; 
metabolites(8).trapLm = .25 ; 
metabolites(8).sat = 3 ; 

metabolites(8).degradeMAX = 1 ; 
metabolites(8).degradeMIN = 0 ; 

%% 9: PAMP - Pathogen cytokine

metabolites(9).type(1:4) = 2 ; % Neumann
metabolites(9).alpha(1:4) = 0 ; % No Flux

metabolites(9).conc = 0*blank_field ; % Initial Concentration Matrix

metabolites(9).Dcoeff = D_antibio*blank_field ; 
metabolites(9).trapDm = .5 ; 
metabolites(9).trapLm = .25 ; 
metabolites(9).sat = 15 ; 

metabolites(9).degradeMAX = 0.75 ;
metabolites(9).degradeMIN = 0.25 ; 

%% 10: Hydrogen

metabolites(10).type(1:4) = 2 ; % Neumann
metabolites(10).alpha(1:4) = 0 ; % No Flux

metabolites(10).conc = 0*blank_field ; % Initial Concentration Matrix

metabolites(10).Dcoeff = D_H2aq*blank_field; 
metabolites(10).trapDm = 0 ; 
metabolites(10).trapLm = 0 ; 
metabolites(10).sat = 0 ; 

metabolites(10).degradeMAX = 0 ; 
metabolites(10).degradeMIN = 0 ; 

%% 11: Dietary Fiber

metabolites(11).type(1:4) = 2 ; % Neumann
metabolites(11).alpha(1:4) = 0 ; % No Flux

metabolites(11).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(11).conc(:,sz) = fiberfeed ; 

metabolites(11).Dcoeff = D_fibers*blank_field ; 
metabolites(11).trapDm = .75 ; 
metabolites(11).trapLm = .5 ; 
metabolites(11).sat = 10 ; 

metabolites(11).degradeMAX = 0 ; 
metabolites(11).degradeMIN = 0 ; 

%% 12: Nitrate

metabolites(12).type(1) = 1 ; % Dirichlet constant epithelial release (value can change during overall iterations but constant during diffusion time)
metabolites(12).type(2:4) = 2 ;  % Neumann
metabolites(12).alpha(1:4) = 0 ; % No Flux

metabolites(12).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(12).conc(:,sz) = fiberfeed ; 

metabolites(12).Dcoeff = D_nitrateAQ*blank_field ; 
metabolites(12).trapDm = 0 ; 
metabolites(12).trapLm = 0 ; 
metabolites(12).sat = 0 ; 

metabolites(12).degradeMAX = 0.25 ; 
metabolites(12).degradeMIN = 0 ; 

%% IgA
metabolites(13).type(1) = 1 ; % Dirichlet constant epithelial release (value can change during overall iterations but constant during diffusion time)
metabolites(13).type(2:4) = 2 ;  % Neumann
metabolites(13).alpha(1:4) = 0 ; % No Flux

metabolites(13).conc = 50*blank_field ; % Initial Concentration Matrix
metabolites(13).conc(:,sz) = fiberfeed ; 

metabolites(13).Dcoeff = D_IgA*blank_field ; 
metabolites(13).trapDm = .90 ; 
metabolites(13).trapLm = .6 ; 
metabolites(13).sat = 40 ;

metabolites(13).degradeMAX = 0.25 ; 
metabolites(13).degradeMIN = 0 ; 

%% MAMPs
metabolites(14).type(1:4) = 2 ; % Neumann
metabolites(14).alpha(1:4) = 0 ; % No Flux

metabolites(14).conc = 0*blank_field ; % Initial Concentration Matrix
metabolites(14).conc(:,sz) = fiberfeed ; 

metabolites(14).Dcoeff = D_cytokines*blank_field ; 
metabolites(14).trapDm = .9 ; 
metabolites(14).trapLm = .5 ; 
metabolites(14).sat = 12 ;  

metabolites(14).degradeMAX = 0 ; 
metabolites(14).degradeMIN = 0 ; 

end 