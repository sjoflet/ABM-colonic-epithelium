% Agent-based model of the human colon to investigate mechanisms of pathogen colonization resistance
% 
% Samantha Fletcher, Carly Ching, Mark Suprenant, Darash Desai, and Muhammad Zaman. 
% Boston University 2025

%% Number values - Model Key

% * Dish * 
% 0 = blank 
% 1 = BT1
% 2 = BT2 
% 3 = BT3 (pathogen)
% 9 = dense mucus
% 9.5 = loose mucus
% 10 = goblet cells 

% * Alive object feature * 
% 0 - dead
% 1 - alive
% * Causes of elimination from model * 
% 15 - cleared
% 17 - starved 
% 18 - killed by antibiotics
% 19 - killed by pathogen
% 20 - Too much oxygen

%% Initialize

folder = "SimulationFolder" ; % folder to save figures and workspace
parentfoldername = "ParentFolder" ; 
mkdir(parentfoldername, folder) % directs to ParentFolder within curent directory
%% 
for countinput = 3:1:6  % for loop to run simulation repeatedly

clearvars -except DEADPATH countinput folder parentfoldername

% Create folder
saveit = 1 ; % 1 = yes or 0 = no to keep folder making code

testtype = 2 ; % 1 = Dietary fiber varying, 2 = Pathogen Varying


% ** Inputs **     
final_graph = 30 ; % figure number for final graph
dish_graph = countinput+final_graph ; % figure number for dish screenshot graphs
plotcolornum = countinput+3 ; % establishes plot color   

        
% general
nsteps = 50 ; % iterations per simulation
inter =  10; % interval to publish dish screenshot on Dish Graph Figure
tot = 100; % Size of dish

% DEFINE FIGURE SIZE
figsNeeded = nsteps/inter+1 ; % Number of intervals plus #1
figrows = ceil(figsNeeded/2) ; % Two figures per row  
iterationrate =  1 ; % min per iteration - changes parameters 
feedhours = 7 ; % add dietary fiber every X hours
fiberFeedRate  = feedhours*60*iterationrate ; % hours*min/hour*iterationsperminute

if testtype == 1     
    fiberOps = [ 0 250 500 750 1000 1250 1500 1750 2000] ; % element number should correspond to countinputs max
    fiberFeedAmt = fiberOps(countinput) ; % selects value based on countinput
elseif testtype == 2 
    fiberFeedAmt = 1000 ; % sets constant amount for testing other variable
end 

% bacterial cells 
bt1 = round(0.04*tot^2) ; % scale initial cells to total size of dish
bt2 = round(0.04*tot^2) ; %                    " " 
BT1_ini = bt1 ; % starting amount of BT1                  
BT2_ini = bt2 ;  % starting amount of BT2
cellstart = bt1+bt2 ; % starting amount of total cells 
flag_PrePath = 0 ; % flag for time prior to pathogen seeding
AllPathDead = 0 ; % Tick Number of complete pathogen resistance

%% Pathogen inputs   

pathseed =  300; % number of seeding tick    

if testtype == 1 
    pathogen_amt = 1000 ; % Set pathogen seeding value for dietary fiber varying simulations

elseif testtype ==2 
    pathops = 1000.*[.5 1 2 3 4 5 ] ; % varying pathogen amounts 
    pathogen_amt = pathops(countinput) ;
end 

bt3 = 0 ; % placeholder 



% Get Parameters - tuned by iteration rate [min]
[paramsBT1, paramsBT2, paramsGC, paramsBT3] = parametersBT1BT2(iterationrate) ; 

%% Figure Labels

if testtype == 1
    legendtitles = "Dietary Fiber Seeded" ;
    Finalgraphtitle = "Dietary Fiber Varying with " +num2str(pathogen_amt) + " Pathogens";
    figname = "DishGraph_" + num2str(fiberFeedAmt) + ".fig" ;
    linelabels =  num2str(fiberFeedAmt) ; 
    DishgraphTitles = "DF = " +  num2str(fiberFeedAmt) ;
    WStitle = "DF" + num2str(fiberFeedAmt) + "_workspace" ;
end

if testtype == 2 
    legendtitles = "Pathogens Seeded" ;
    Finalgraphtitle = "Pathogens Varying with " + num2str(fiberFeedAmt) +" DF" ;
    figname = "DishGraph_" + num2str(pathogen_amt) + ".fig" ;
    linelabels = num2str(pathogen_amt) ; 
    DishgraphTitles = "Pathogens = " + num2str(pathogen_amt) ; 
    WStitle = "DF" + num2str(pathogen_amt) + "_workspace" ;    
end 

if saveit == 1 
workspacetitle = parentfoldername + "/" +folder+ "/" + WStitle  ;
end 

% Establish Geometry
dish = 9*ones(tot,tot) ; % impenetrable mucus = 9 -> start with it everywhere 
dish_labels = zeros(tot,tot) ; % Matrix for agent #s, no agents yet 

% ** Metabolites ** 

[metabolites, mets] = MetaboliteStruct(tot,fiberFeedAmt) ; 
del = [1 ; 1 ; 1] ; % [dt = 1s, dx = 1 um, dy = 1 um] 
tfinal = iterationrate*60 ; % mins to secs

% Bkron = ones(5,5) ; 

%% ** Seed Cells ** 

% bacteria
BT1locN = round(1+(tot-1)*rand(bt1,1)) ; % can be any row
BT2locN = round(1+(tot-1)*rand(bt2,1)) ; % "

BT1locM = round(round(tot*.25)+(tot-round(tot*.25))*rand(bt1,1)) ; % can be anywhere except first 50 impenetrable layers (columns)
BT2locM = round(round(tot*.75)+(tot-round(tot*.75))*rand(bt2,1)) ; % outer layer starts at 75% to prevent immediate oxygen death

BT1loc = [BT1locN BT1locM] ; % pair randomly generated rows and columns 
BT2loc = [BT2locN BT2locM ] ; % pair randomly generated rows and columns

for count = 1:bt1 
    pos = BT1loc(count, :) ; 
    while dish(pos(1), pos(2)) ~= 9
    pos(1) = round(1+(tot-1)*rand(1,1)) ; 
    pos(2) = round(round(tot*.25)+(tot-round(tot*.25))*rand(1,1)) ; 
    end
    [dish, dish_labels] = Adjust(count,dish,dish_labels, 1, pos(1),pos(2)) ;  % override MUC2 at points and 1 == BT1 and place agent numbers at their location
    last_rep = ceil(rand*paramsBT1.rep_rate) ;
    food = ceil(rand*2) ; 
    agent{count} = BT1(paramsBT1.rep_rate,last_rep,pos,food,1,1) ; % alive = 1 for yes
end 

count1 = 1 ; 
for count2 = (bt1+1):(bt1+bt2)
    pos = BT2loc(count1,:) ; 
    while dish(pos(1), pos(2)) ~= 9
       pos(1) = round(1+(tot-1)*rand(1,1)) ; 
       pos(2) = round(round(tot*.60)+(tot-round(tot*.60))*rand(1,1)) ; 
    end
    [dish, dish_labels] = Adjust(count2, dish,dish_labels, 2, pos(1),pos(2)) ;  % override MUC2 at points and 2 == Anaerobic bac and place agent numbers at their location
    last_rep = ceil(rand*paramsBT2.rep_rate) ; 
    food = ceil(rand*2) ; 
    agent{count2} = BT2(paramsBT2.rep_rate, last_rep, pos,food, 1) ; %alive = 1 for yes 
    count1 = count1 + 1 ;   
end

% Create goblet cells 
N = 1:1:tot ; 
dummycol = 10*ones(tot,1) ; 
muc2_productionrate = paramsGC.mucProd ; 

count1 = 1 ; 
for count3 = (bt1+bt2+1):(bt1+bt2+tot)
    pos = [count1 1] ; 
    randlastProd = ceil(rand*muc2_productionrate) ;
    randButTot = ceil(rand*paramsGC.butMin) ; 
    agent{count3} = Goblet(muc2_productionrate,randlastProd, pos, randButTot, 1) ; % alive = 1 for yes
    count1 = count1 + 1 ;
end 

dish(:,1) = dummycol ; % add goblet cell column to dish column 1 

%  ** Data Storage ** 
% Record total population values after each iteration
BT1_TOTPOP = zeros(nsteps,1) ; 
BT2_TOTPOP = zeros(nsteps,1) ; 
BT3_TOTPOP = zeros(nsteps,1) ; 
DISH = struct("dishvals", cell(1,nsteps)) ; 
DenseMucusThickness = zeros(nsteps,1) ; 
LooseMucusThickness = zeros(nsteps,1) ; 
METABOLITES(1:mets,nsteps) = struct('conc', zeros(tot,tot)) ; 

% Record individual time of replication or mucus production
BT1Reps = [] ; 
BT2Reps = [] ; 
BT3Reps = [] ; 
GC_MUCProd= [] ; 

% ** Plotting Features ** 

Labels = 1:1:tot ; 
cut_size = tot/2 ; 
MyLabels = string(Labels) ; 
MyLabels(mod(Labels,cut_size) ~= 0) = " " ; 


% Why are the Cells Dying ? 

BT1_die = struct("antibio", 0, "starve",0, "cleared",0) ; 
BT2_die = struct("oxygen", 0, "starve",0, "cleared",0) ; 
BT3_die = struct("antipath",0, "starve",0, "cleared",0);

%% Run

plotcount = 1 ; % start plot count for dish graph subplots

for tick = 1:nsteps

    % ** Add fiber ** 
    if rem(tick,fiberFeedRate) == 0 
        metabolites(11).conc(:,tot) = metabolites(11).conc(:,tot) + fiberFeedAmt ; 
    end 

    % ** Add pathogens 
    if tick >= pathseed && pathogen_amt > 0
        [agent,dish,dish_labels,bt3,bt1,bt2,overwrite] = addpathogen(pathogen_amt,tot,dish, dish_labels,agent,paramsBT3,bt1, bt2) ;
        pathogen_amt = 0 ;
        flag_PrePath = 1 ; % flag the end of Pre-pathogen time
    end 

    prev_agnts = length(agent) ; % find number of agents 
    agnts_new = 0 ; % reset new agent count tracker for replication and seeding 
    dummycol = 1:prev_agnts ; 
    randomized_agents = dummycol(randperm(length(dummycol))) ; % list of agent numbers in random order to not always start with the same one (mix-up opportunities for resources) 

    %% ** Metabolite Trapping and Diffusion **

    % find locations of dense and loose mucus for metabolite trappig
    dish_Dm = zeros(tot,tot) ; 
    dish_Lm = zeros(tot,tot) ; 
    dish_blank = zeros(tot,tot) ; 

    dummyfindDM = find(dish==9) ; 
    dummyfindLM = find(dish==9.5) ; 

    dish_Dm(dummyfindDM) = 1 ;
    dish_Lm(dummyfindLM) = 1 ; 

    % Trap Metabolites and Diffuse Excess
    for metcount = 1:1:mets
        clear traptemp trapamt diffConc tempConc degradeConc
        randomDegrade = 1 - ( metabolites(metcount).degradeMIN + (metabolites(metcount).degradeMAX - metabolites(metcount).degradeMIN)*rand );
        degradeConc = randomDegrade*metabolites(metcount).conc ;
        if metabolites(metcount).sat ~= 0
            trapamt = dish_Dm*metabolites(metcount).trapDm + dish_Lm*metabolites(metcount).trapLm ;
            traptemp = trapamt.*degradeConc ;
            if max(traptemp) > metabolites(metcount).sat 
                dummyfindSat = find(traptemp > metabolites(metcount).sat) ; 
                traptemp(dummyfindSat) = metabolites(metcount).sat ;
            end 
            tempConc = degradeConc - traptemp ;
            [diffConc] = Diffuse_TDMAvec(tempConc, metabolites(metcount).Dcoeff, tot, tfinal, del, metabolites(metcount).type, metabolites(metcount).alpha) ;
            metabolites(metcount).conc = diffConc + traptemp ;
        else 
            [metabolites(metcount).conc] = Diffuse_TDMAvec(degradeConc, metabolites(metcount).Dcoeff, tot, tfinal, del, metabolites(metcount).type, metabolites(metcount).alpha) ;
        end 
          
    end 

    %% ** Start agent iterations ** 
    for num = 1:prev_agnts
        cn = randomized_agents(num) ; % take next random agent number from generated list
        curr = agent{cn} ; % current agent
        [dish,metabolites(3).conc,metabolites(5).conc, metabolites(4).conc] = mucus_environment(dish, metabolites(5).conc, metabolites(4).conc, tot, metabolites(3).conc, paramsGC.MucusBox) ; 

   %% BT1
        if isa(curr, 'BT1')
            if curr.alive == 1

                cpos = curr.position ; 
                cbiomass = curr.biomass ; 
                crep_rate = curr.rep_rate ; 
                clast_r = curr.last_rep ; 

                % Antibiotics
                if metabolites(6).conc(cpos(1),cpos(2)) > paramsBT1.antibio
                    [dish, dish_labels] = Adjust(cn, dish, dish_labels, 0, cpos(1), cpos(2)) ; 
                    bt1 = bt1 - 1 ;  
                    BT1_die.antibio = BT1_die.antibio + 1 ; 
                    metabolites(6).conc(cpos(1),cpos(2)) = 0 ; 
                    agent{cn} = BT1(crep_rate,clast_r,[], cbiomass, 0, 18) ; % not anchored and not alive
                    continue
                end 

                % Oxygen
                if metabolites(7).conc(cpos(1),cpos(2)) > paramsBT1.O2 && metabolites(7).conc(cpos(1),cpos(2)) < paramsBT1.O2consume 
                    metabolites(7).conc(cpos(1),cpos(2)) = 0 ; % subtract amount consumed by bacteria
                    if crep_rate > paramsBT1.repMin 
                        crep_rate = crep_rate - paramsBT1.O2Rep ; % faster when aerobic resp
                    end 
                elseif  metabolites(7).conc(cpos(1),cpos(2)) >= paramsBT1.O2consume 
                    metabolites(7).conc(cpos(1),cpos(2)) = metabolites(7).conc(cpos(1),cpos(2)) - paramsBT1.O2consume  ; % subtract max
                    if crep_rate > paramsBT1.repMin 
                        crep_rate = crep_rate - paramsBT1.O2Rep ; % faster when aerobic resp
                    end 
                else 
                    crep_rate = crep_rate  + paramsBT1.O2Rep ; % slower when anaerobic resp
                end 

                % Eat 
                if metabolites(2).conc(cpos(1),cpos(2)) > 0 
                    [temp_galtirate] = MM(metabolites(2).conc(cpos(1),cpos(2)),paramsBT1.galtiVmax,paramsBT1.galtiKm) ; % Michaelis-Menten calculation for uptake rate at the time
                    metabolites(2).conc(cpos(1),cpos(2)) = metabolites(2).conc(cpos(1),cpos(2)) - temp_galtirate ; % lower instant location concentration
                    cbiomass = cbiomass + temp_galtirate ;
                    clear temp_galtirate
                    if crep_rate > paramsBT1.repMin 
                        crep_rate = crep_rate - paramsBT1.galtiRep ; % faster when fed -> allows recovery! 
                    end 
                elseif metabolites(1).conc(cpos(1),cpos(2)) > 0 
                    [temp_galacrate] = MM(metabolites(1).conc(cpos(1),cpos(2)),paramsBT1.galacVmax,paramsBT1.galacKm) ; 
                    metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) - temp_galacrate ; 
                    cbiomass = cbiomass + temp_galacrate ; 
                    clear temp_galacrate
                    if crep_rate < paramsBT1.galacRepMax 
                        crep_rate = crep_rate + paramsBT1.galacRep ; 
                    end 
                else 
                    cbiomass = cbiomass - iterationrate*paramsBT1.biomassLost ;
                    crep_rate = crep_rate + paramsBT1.starveRep ; % replication slows at no food 
                end 

                if cbiomass <= paramsBT1.starve % dies 
                    [dish, dish_labels] = Adjust(cn, dish, dish_labels, 0, cpos(1), cpos(2)) ; 
                    bt1 = bt1 - 1 ; 
                    BT1_die.starve = BT1_die.starve + 1 ; 
                    agent{cn} = BT1(crep_rate,clast_r,[], cbiomass, 0, [17 tick]) ; % not anchored and not alive
                    continue 
                end 

                % Bacteria Proteolytic Factors
                metabolites(4).conc(cpos(1),cpos(2)) = metabolites(4).conc(cpos(1),cpos(2)) + paramsBT1.proteo ;
                
                % MAMPs
                metabolites(14).conc(cpos(1),cpos(2)) = metabolites(14).conc(cpos(1),cpos(2)) + paramsBT1.MAMPs ;
     
                % Iga
                metabolites(13).conc(cpos(1),cpos(2)) = 0 ; % IgA Binds Cells 
                [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9,dish,tot, 3) ;
                if check == 0
                    if size(M_spot,1) == 2 && randsample(paramsBT1.degradationOdds,1) == 1 % 1 in degradation odds chance of eating mucus 
                        metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox ; % breakdown mucus and add galactose 
                        dish(M_spot(1,1), M_spot(1,2)) = 0 ; % remove mucus from dish
                     elseif size(M_spot,1) == 3
                         metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox ; % breakdown mucus and add galactose 
                         dish(M_spot(1,1), M_spot(1,2)) = 0 ; % remove mucus from dish
                         dish(M_spot(2,1), M_spot(2,2)) = 9.5 ; % add loose mucus to other spot
                    end 
                end 

                [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9.5,dish,tot, 1) ;
                if check == 1 % No anchoring
                    % First: Eat
                    if metabolites(2).conc(cpos(1),cpos(2)) > 0 % eat again to make up for no anchoring 
                        [temp_galtirate] = MM(metabolites(2).conc(cpos(1),cpos(2)),paramsBT1.galtiVmax,paramsBT1.galtiKm); % Michaelis-Menten calculation for consumption rate at the time
                        metabolites(2).conc(cpos(1),cpos(2)) = metabolites(2).conc(cpos(1),cpos(2)) - temp_galtirate ; % lower instant location concentration
                        clear temp_galtirate
                    elseif metabolites(1).conc(cpos(1),cpos(2)) > 0
                        [temp_galacrate] = MM(metabolites(1).conc(cpos(1),cpos(2)),paramsBT1.galacVmax,paramsBT1.galacKm) ;
                        metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) - temp_galacrate ;
                        cbiomass = cbiomass + temp_galacrate ;
                        clear temp_galacrate
                        if crep_rate < paramsBT1.galacRepMax 
                            crep_rate = crep_rate + paramsBT1.galacRep ;
                        end
                    else 
                        if crep_rate < paramsBT1.AnchorRepMax
                            crep_rate = crep_rate - paramsBT1.AnchorRep ; % slow down replication
                        end
                    end
                    
                    % Second: Move and release PAMPs
                    metabolites(9).conc(cpos(1),cpos(2)) = metabolites(9).conc(cpos(1),cpos(2)) + paramsBT1.PAMPs ; % release PAMPs when using flagella

                    [check,opencoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,1) ; % no anchoring so want to move into an empty space
                    if check == 1 % no empty spaces to move into
                        anchor = 0 ; 
                    else 
                        dish(cpos(1), cpos(2)) = 0 ; %get out of old space
                        dish_labels(cpos(1),cpos(2)) = 0 ; %get out of old space
                        agent{cn}.position = [opencoors(1) opencoors(2) ] ; %move into random nearby open space
                        cpos = agent{cn}.position ; %rewrite cpos
                        [dish, dish_labels] = Adjust(cn,dish,dish_labels, 1, cpos(1),cpos(2)) ; 
                        [check,~] = Cell_neighbors(cpos(1),cpos(2),9.5,dish,tot, 2) ; %look again
                        if check == 1 
                            if crep_rate < paramsBT1.AnchorRepMax
                                crep_rate = crep_rate - paramsBT1.AnchorRep ; % slow down replication from moving
                            end 
                            cbiomass = cbiomass - paramsBT1.noAnchorFood ; % loose biomass from moving
                            anchor = 0 ; 
                        else 
                            anchor = 1 ; 
                        end
                    end 
                else
                     anchor = 1 ; 
                end 

                % release antipathogens and release more antipathogens if there are pathogens nearby
                metabolites(8).conc(cpos(1),cpos(2)) = metabolites(8).conc(cpos(1),cpos(2)) + paramsBT1.antipath*(1+metabolites(9).conc(cpos(1),cpos(2))) ;
                metabolites(9).conc(cpos(1),cpos(2)) =  metabolites(9).conc(cpos(1),cpos(2))*.5 ;

                % Replicate if its been enough time 
                if tick-clast_r >= crep_rate
                    [isnew, newcoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,1) ; % Looking for empty space to go into
                    if isnew == 0 % space to replicate! 
                        actrep = tick-clast_r ; 
                        BT1Reps = vertcat(BT1Reps, actrep) ; 
                        crep_rate = crep_rate - paramsBT1.parentpoleRep ; % parental poles for OG cell ; 
                        clast_r = tick ; % rewrite last rep
                        agnts_new = agnts_new + 1 ; 
                        split = cbiomass/2 ; 
                        cbiomass = split ; 
                        posnew = [newcoors(1) newcoors(2)] ; 
                        bt1 = bt1 + 1 ;
                        [check, M_spot] = Cell_neighbors(newcoors(1),newcoors(2),9,dish,tot,1) ; 
                        if check == 1 
                            anchor = 0 ; 
                        else 
                            anchor = 1 ; 
                        end 
                        agent{prev_agnts+agnts_new} = BT1(paramsBT1.rep_rate,tick, posnew, cbiomass, anchor, 1) ; 
                        [dish, dish_labels] = Adjust(prev_agnts+agnts_new,dish,dish_labels, 1, newcoors(1),newcoors(2)) ; 
                    end 
                end

                % End Matters 

                agent{cn} = BT1(crep_rate, clast_r, cpos, cbiomass, anchor, 1) ;

            else 
                continue
            end
 
            %% Anaerobic BT2

        elseif isa(curr, 'BT2')
            if curr.alive == 1
                cpos = curr. position ; 
                cbiomass = curr.biomass ; 
                crep_rate = curr.rep_rate ; 
                clast_r = curr.last_rep ; 

                % Oxygen
                if metabolites(7).conc(cpos(1),cpos(2)) > paramsBT2.O2 % death condition
                    [dish, dish_labels] = Adjust(cn,dish,dish_labels, 0, cpos(1),cpos(2)) ; 
                    bt2 = bt2 - 1 ; 
                    metabolites(7).conc(cpos(1),cpos(2)) = 0 ; 
                    BT2_die.oxygen = BT2_die.oxygen + 1; 
                    agent{cn} = BT2(crep_rate,clast_r,[], cbiomass, 20) ; % not anchored and not alive
                    continue
                end 
                
                % Eat and Produce
                
                % Degrade for Food  

                if metabolites(11).conc(cpos(1),cpos(2)) > paramsBT2.fiberMin
                    [fiberconsumedtemp] = MM(metabolites(11).conc(cpos(1),cpos(2)),paramsBT2.dfVmax,paramsBT2.dfKm) ; 
                    metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + fiberconsumedtemp*paramsBT2.fiberdegrade ;
                    metabolites(11).conc(cpos(1),cpos(2)) = metabolites(11).conc(cpos(1),cpos(2)) - fiberconsumedtemp ; 
                    clear fiberconsumedtemp
                else
                    % breakdown mucus
                    [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9.5,dish,tot,8) ;
                    if check == 0 % mucus is near 
                        possibilities = size(M_spot,1) ; 
                        mucusdegrad = min([possibilities, randsample(paramsBT2.degradationOdds,1)]) ; 
                        for countmd = 1:1:mucusdegrad
                            metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox ; % galactose released from mucus
                            dish(M_spot(countmd,1),M_spot(countmd,2)) = 0 ; % no more mucus at spot
                        end 
                    end

                    [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9,dish,tot,8) ;
                    if check == 0 % mucus is near 
                        possibilities = size(M_spot,1) ; 
                        mucusdegrad = min([possibilities, randsample(paramsBT2.degradationOdds,1)]) ; 
                        for countmd = 1:1:mucusdegrad
                        metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox*2 ; % galactose released from mucus
                        dish(M_spot(countmd,1),M_spot(countmd,2)) = 0 ; % no more mucus at spot
                        end 
                    end
                end

                % Eat 

                if metabolites(1).conc(cpos(1),cpos(2)) > 0 
                    [temp_galacrate] = MM(metabolites(1).conc(cpos(1),cpos(2)),paramsBT2.galacVmax,paramsBT2.galacKm) ; % Michaelis-Menten calculation for consumption rate at the time
                    metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) - temp_galacrate ; % lower instant location concentration
                    metabolites(3).conc(cpos(1),cpos(2)) = metabolites(3).conc(cpos(1),cpos(2)) + paramsBT2.but*temp_galacrate ; % yield of butyrate from galactose 
                    metabolites(2).conc(cpos(1),cpos(2)) = metabolites(2).conc(cpos(1),cpos(2)) + paramsBT2.galti*temp_galacrate ; % yield of galactitol from galactose 
                    cbiomass = cbiomass + temp_galacrate ;
                    if crep_rate > paramsBT2.repMin 
                        crep_rate = crep_rate - paramsBT2.galacRep ; % faster when fed -> allows recovery! 
                    end 
                else 
                    cbiomass = cbiomass - iterationrate*paramsBT2.biomassLost ;
                    crep_rate = crep_rate + paramsBT2.galacRep ; % replication slows at no food 
                end 

                if cbiomass <= paramsBT2.starve % dies 
                    [dish, dish_labels] = Adjust(cn,dish,dish_labels, 0, cpos(1),cpos(2)) ; 
                    bt2 = bt2 - 1 ; 
                    BT2_die.starve = BT2_die.starve + 1; 
                    agent{cn} = BT2(crep_rate,clast_r,[], cbiomass, 0) ; 
                    continue % go to next agent 
                end 

                % move at random
                moveit = randi([1 5],5) ; % 1 in 5 chance of random movement
                if moveit == 3 
                    [check,opencoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,1) ; 
                    if check == 0
                        dish(cpos(1), cpos(2)) = 0 ; %get out of old space
                        dish_labels(cpos(1),cpos(2)) = 0 ; %get out of old space
                        agent{cn}.position = [opencoors(1) opencoors(2) ] ; %move into random nearby open space
                        cpos = agent{cn}.position ; %rewrite cpos
                        [dish, dish_labels] = Adjust(cn,dish,dish_labels, 2, cpos(1),cpos(2)) ; 
                    end 
                end

                % Bacteria Proteolytic Factors
                metabolites(4).conc(cpos(1),cpos(2)) = metabolites(4).conc(cpos(1),cpos(2)) + paramsBT2.proteo ; 

                % MAMPs
                metabolites(14).conc(cpos(1),cpos(2)) = metabolites(14).conc(cpos(1),cpos(2)) + paramsBT2.MAMPs ;
     
                % IgA
                metabolites(13).conc(cpos(1),cpos(2)) = 0 ; % IgA Binds Cells 

                % Antibacterials 
                [check,opencoors] = Cell_neighbors(cpos(1),cpos(2),1,dish,tot,8) ; % Looking for surrounding BT1
                if check == 0 
                    numB1 = size(opencoors,1) ; 
                elseif check == 1 
                    numB1 = 0 ; 
                end 

                metabolites(6).conc(cpos(1),cpos(2)) = metabolites(6).conc(cpos(1),cpos(2)) + (1+numB1*10)*paramsBT2.antibio ;

                % Antipathogens
                metabolites(8).conc(cpos(1),cpos(2)) = metabolites(8).conc(cpos(1),cpos(2)) + paramsBT2.antipath*(1+metabolites(9).conc(cpos(1),cpos(2))) ;
                metabolites(9).conc(cpos(1),cpos(2)) =  metabolites(9).conc(cpos(1),cpos(2))*.5 ;

               
                % Replicate if its been enough time 
                if tick-clast_r >= crep_rate
                    [isnew, newcoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,1) ; % Looking for empty space to go into
                    if isnew == 0 % space to replicate!
                        actrep = tick-clast_r ; 
                        BT2Reps = vertcat(BT2Reps, actrep) ; 
                        crep_rate = crep_rate - paramsBT2.parentpoleRep ; % parental poles for OG cell ; 
                        clast_r = tick ; % rewrite last rep
                        agnts_new = agnts_new+1 ; 
                        split = cbiomass/2 ; 
                        cbiomass = split ; 
                        posnew = [newcoors(1) newcoors(2)] ; 
                        bt2 = bt2 + 1 ;
                        agent{prev_agnts+agnts_new} = BT2(paramsBT2.rep_rate,tick, posnew, cbiomass, 1) ; 
                        [dish, dish_labels] = Adjust(prev_agnts+agnts_new,dish,dish_labels, 2, newcoors(1),newcoors(2)) ;
                    else % make room 
                        metabolites(6).conc(cpos(1),cpos(2)) = metabolites(6).conc(cpos(1),cpos(2)) + paramsBT2.antibio*10 ;
                    end 
                end
                % End Matters 
                agent{cn} = BT2(crep_rate, clast_r, cpos, cbiomass, 1) ; 
            else 
                continue
            end

            %% Goblet Cells 
        elseif isa(curr, 'Goblet')
            if curr.alive == 1
                 cpos = curr.position ; 
                 m_rate = curr.muc2_prodrate  ;
                 lastRep = curr.muc2_lastProd ; 
                 cButyTot = curr.butyrateTot ; 

                 % Consume Butyrate 

                 if metabolites(3).conc(cpos(1),cpos(2)) > 0 
                     [temp_butyrate] = MM(metabolites(3).conc(cpos(1),cpos(2)), paramsGC.butVmax, paramsGC.butKm) ; 
                     metabolites(3).conc(cpos(1),cpos(2)) = metabolites(3).conc(cpos(1),cpos(2)) - temp_butyrate ; 
                     cButyTot = cButyTot+temp_butyrate ; 
                     clear temp_butyrate
                     if cButyTot > paramsGC.butMin 
                         metabolites(7).conc(cpos(1),cpos(2)) = paramsGC.O2 ;  % reset O2 with proper butyrate
                         cButyTot = cButyTot-paramsGC.butMin ; 
                         if m_rate > paramsGC.butScale+paramsGC.mucMin 
                            m_rate = m_rate - paramsGC.butScale ; % butyrate encourages mucus production
                         end
                     end 
                     agent{cn}.butyrateTot = cButyTot ; 
                 else 
                     metabolites(7).conc(cpos(1),cpos(2)) = 2*metabolites(7).conc(cpos(1),cpos(2)) ; % Double O2 released 
                     m_rate = m_rate + paramsGC.noBut;
                 end 

                 agent{cn}.muc2_prodrate = m_rate ; 
                
                 % IgA and Inflammation
                if metabolites(13).conc(cpos(1),cpos(2)) > paramsGC.IgAMin % scale inflammation based on previous 
                    % recovery
                    if metabolites(7).conc(cpos(1),cpos(2)) > paramsGC.O2 && metabolites(6).conc(cpos(1),cpos(2)) > paramsGC.antibioReg && ...
                            metabolites(12).conc(cpos(1),cpos(2)) > paramsGC.NitReg 
                        metabolites(7).conc(cpos(1),cpos(2)) = metabolites(7).conc(cpos(1),cpos(2)) - paramsGC.O2rec ;
                        metabolites(9).conc(cpos(1),cpos(2)) = 2*metabolites(8).conc(cpos(1),cpos(2)) + paramsGC.antipath*(1+metabolites(9).conc(cpos(1),cpos(2))); % overreaction during recovery

                     % maintainence 
                    else 
                        metabolites(7).conc(cpos(1),cpos(2)) = paramsGC.O2 ; 
                        metabolites(6).conc(cpos(1),cpos(2)) = paramsGC.antibioReg ; 
                        metabolites(8).conc(cpos(1),cpos(2)) = metabolites(8).conc(cpos(1),cpos(2)) + paramsGC.antipath*(1+metabolites(9).conc(cpos(1),cpos(2))) ;  % scale anti-pathogens by PAMPs 
                    end 
                else % Inflammation
                    metabolites(7).conc(cpos(1),cpos(2)) = 10^3* metabolites(7).conc(cpos(1),cpos(2)) ; % significant oxygen increase
                    metabolites(7).conc(cpos(1),cpos(2)) = 10 * metabolites(6).conc(cpos(1),cpos(2)) ; % significant antibiotic increase
                    metabolites(7).conc(cpos(1),cpos(2)) = 10 * metabolites(8).conc(cpos(1),cpos(2)) ; % significant antipathogen increase
                    metabolites(7).conc(cpos(1),cpos(2)) = 10 * metabolites(12).conc(cpos(1),cpos(2)) ; % significant nitrate increase
                end 

                % PAMP-based immune response always active 
                metabolites(8).conc(cpos(1),cpos(2)) = metabolites(8).conc(cpos(1),cpos(2)) + paramsGC.antipath*(1+metabolites(9).conc(cpos(1),cpos(2))) ;  % scale anti-pathogens by PAMPs 

                % Increase IgA by MAMPs
                metabolites(13).conc(cpos(1),cpos(2)) = metabolites(13).conc(cpos(1),cpos(2)) + paramsGC.IgAscale*(1+metabolites(14).conc(cpos(1),cpos(2))) ; 
                if  metabolites(13).conc(cpos(1),cpos(2)) > paramsGC.IgAMax 
                    metabolites(13).conc(cpos(1),cpos(2)) = paramsGC.IgAMax ; 
                end 

                 if metabolites(9).conc(cpos(1),cpos(2)) >= paramsGC.epiMin
                     O2multi = floor(metabolites(9).conc(cpos(1),cpos(2))/paramsGC.epiMin) ; 
                     metabolites(7).conc(cpos(1),cpos(2)) = metabolites(7).conc(cpos(1),cpos(2)).*(1+O2multi) ; % leaky oxygen 
                     metabolites(9).conc(cpos(1),cpos(2)) = 0 ; %reset PAMPS 
                 end 

                 % Human Proteolytic Factors
                 metabolites(5).conc(cpos(1),cpos(2)) = metabolites(5).conc(cpos(1),cpos(2)) + paramsGC.proteo ; 

                 % curr.muc2_prodrate = m_rate ; 
                 mucamt = ceil(round(rand*10)) ; 
                 

                %%  secrete mucus 
                 if tick-lastRep >= m_rate % it has been long enough to produce mucus
                     actGCrep = tick-lastRep ; 
                     GC_MUCProd = vertcat(GC_MUCProd,actGCrep ) ; 
                     agent{cn}.muc2_lastProd = tick ; % rewrite last tick so that you can have a current value
                     for count = 1:1:mucamt % mucus amount: want cells to be able to produce >1 mucus
                         clear dummyvex dummysize extract_dish extract_labels implant_labels
                         extract_dish = dish(cpos(1), 2:tot) ; % take out row of interest except for goblet cell itself
                         dummyvex(1,1) = 9 ; % start with a mucus
                       
                         for count2 = 1:1:tot-2  % look at every spot (goblet cell isnt there so there is no "tot" spot, don't need tot-1 because dummy vex shouldnt be longer than tot-1
                             dummyvex(1,count2+1) = extract_dish(1,count2) ; % move spot one over in dummy vex
                         end 

                         dummysize = size(dummyvex,2) ;
                         lastcol = dummysize+1 ; %first column is goblet cell
                         dish(cpos(1),2:lastcol) = dummyvex ; % dish adjusted and mucus moved

                         % fix labels and adjust agents
                         extract_labels = dish_labels(cpos(1),2:lastcol)  ; 

                         for count2 = 1:1:dummysize

                            if extract_labels(1,count2) ~= 0 
                                agntnum = extract_labels(1,count2) ; 
                                 Prepos = agent{agntnum}.position ; 
                                 if Prepos(2) == tot % last cell that is cleared
                                     newCol = [ ] ; 
                                     agent{agntnum}.alive = tick ; % clearance #

                                     if isa(agent{agntnum}, 'BT1')
                                        bt1 = bt1 - 1; 
                                        BT1_die.cleared = BT1_die.cleared + 1; 
                                     end 

                                     if isa(agent{agntnum}, 'BT3')
                                        bt3 = bt3 - 1; 
                                        BT3_die.cleared = BT3_die.cleared + 1; 
                                     end 

                                     if isa(agent{agntnum}, 'BT2') 
                                        bt2 = bt2 - 1 ; 
                                        BT2_die.cleared = BT2_die.cleared + 1; 
                                     end 

                                 else 
                                     newCol = [ Prepos(1) (Prepos(2)+1)] ; 
                                 end 
                                 agent{agntnum}.position = newCol ; % rewrite location in object
                            end
                         end 
                         extract_labels(1,dummysize) = 0 ; %override last place
                         implant_labels = circshift(extract_labels,1) ; % move everything one space
                         dish_labels(cpos(1), 2:lastcol) = implant_labels  ; % labels moved
                      end 
                 end 
            else
                continue
            end  

        %% Pathogenic BT3
        elseif isa(curr, 'BT3')
            if curr.alive == 1

                crep_rate = curr.rep_rate ; 
                clast_r = curr.last_rep ;
                cpos = curr.position ; 
                cbiomass = curr.biomass ; 

                % Antipathogens
                if metabolites(8).conc(cpos(1),cpos(2)) > paramsBT3.antipath % could introduce randomization! might be too finite  
                    [dish, dish_labels] = Adjust(cn, dish, dish_labels, 0, cpos(1), cpos(2)) ; 
                    bt3 = bt3 - 1 ;  
                    BT3_die.antipath = BT3_die.antipath + 1 ; 
                    metabolites(8).conc(cpos(1),cpos(2)) = 0 ; %reset space
                    agent{cn} = BT3(crep_rate,clast_r, cbiomass,[], 18) ; % killed by antipathogen
                    continue
                end 

                % Oxygen
                if metabolites(7).conc(cpos(1),cpos(2)) > metabolites(7).conc(cpos(1),cpos(2))*paramsBT3.O2consume % don't want to go into negatives
                    metabolites(7).conc(cpos(1),cpos(2)) = metabolites(7).conc(cpos(1),cpos(2)) - metabolites(7).conc(cpos(1),cpos(2))*paramsBT3.O2consume ; % subtract amount consumed by bacteria
                    if crep_rate > paramsBT3.repMin 
                        crep_rate = crep_rate - paramsBT3.O2Rep ; %faster when aerobic resp
                    end 
                end 

                % Bacteria Proteolytic Factors
                metabolites(4).conc(cpos(1),cpos(2)) = metabolites(4).conc(cpos(1),cpos(2)) + paramsBT3.proteo ; 

                % PAMPs
                metabolites(9).conc(cpos(1),cpos(2)) = metabolites(9).conc(cpos(1),cpos(2)) + paramsBT3.PAMPs ; 


                % IgA
                metabolites(13).conc(cpos(1),cpos(2)) = 0 ; % IgA Binds Cells 
                
                % Eat 
                if metabolites(2).conc(cpos(1),cpos(2)) > 0 
                    temp_galtirate = MM(metabolites(2).conc(cpos(1),cpos(2)),paramsBT3.galtiVmax,paramsBT3.galtiKm) ; % Michaelis-Menten calculation for consumption rate at the time
                    metabolites(2).conc(cpos(1),cpos(2)) = metabolites(2).conc(cpos(1),cpos(2)) - temp_galtirate ; % lower instant location concentration
                    cbiomass = cbiomass + temp_galtirate ;
                    if crep_rate > paramsBT3.repMin 
                        crep_rate = crep_rate - paramsBT3.galtiRep ; % faster when fed -> allows recovery! 
                    end 
                    clear temp_galtirate
                elseif metabolites(1).conc(cpos(1),cpos(2)) > 0 
                    temp_galacrate = MM(metabolites(1).conc(cpos(1),cpos(2)),paramsBT3.galacVmax,paramsBT3.galacKm) ; 
                    metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) - temp_galacrate ; 
                    cbiomass = cbiomass + temp_galacrate ;
                    clear temp_galacrate
                else 
                    cbiomass = cbiomass - iterationrate*paramsBT3.biomassLost ;
                    crep_rate = crep_rate - paramsBT3.galtiRep ; % replication slows at no food 
                end 

                if cbiomass <= paramsBT3.starve % dies 
                    [dish, dish_labels] = Adjust(cn, dish, dish_labels, 0, cpos(1), cpos(2)) ; 
                    bt3 = bt3 - 1 ; 
                    BT3_die.starve = BT3_die.starve + 1 ; 
                    agent{cn} = BT3(crep_rate,clast_r, cbiomass, [], 17) ; % not anchored and not alive
                    continue 
                end 

                % Nitrate boost
                if metabolites(12).conc(cpos(1),cpos(2)) > 0
                    metabolites(12).conc(cpos(1),cpos(2)) = 0 ; 
                    if crep_rate > paramsBT3.repMin
                        crep_rate = crep_rate - paramsBT3.galtiRep ; % faster when fed -> allows recovery!
                    end
                end 

                % release antibiotics if there are BT1 nearby
                [check, posECs] = Cell_neighbors(cpos(1),cpos(2),1,dish,tot,8) ; % Look if there is an BT1 next to the pathogen (can be a max of 8)
                if check == 0   
                    numECs = size(posECs,1) ; % row number represents number of surrounding pathogens
                    metabolites(6).conc(cpos(1),cpos(2)) = metabolites(6).conc(cpos(1),cpos(2)) + numECs*paramsBT3.antibio ; % scale added antipathogens to pathogen number
                end 

                % breakdown mucus
                if randsample(paramsBT3.densemucDegr,1) == 1 
                    [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9.5,dish,tot,1) ;
                    if check == 0 % mucus is near
                        metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox ; % galactose released from mucus
                        dish(M_spot(1,1),M_spot(1,2)) = 0 ; % no more mucus at spot
                    end
                end 

                if randsample(paramsBT3.loosemucDegr,1) == 1 
                    [check, M_spot] = Cell_neighbors(cpos(1),cpos(2),9,dish,tot,1) ;
                    if check == 0 % mucus is near
                        metabolites(1).conc(cpos(1),cpos(2)) = metabolites(1).conc(cpos(1),cpos(2)) + paramsGC.MucusBox*2 ; % galactose released from mucus
                        dish(M_spot(1,1),M_spot(1,2)) = 0 ; % no more mucus at spot
                    end
                end 

                % Movement
                currspeed = paramsBT3.speed ; 
                while currspeed > 0 
                    cpos = curr.position ; %reset cpos to fit current location
                    [check,opencoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,8) ; 
                    [loosecheck,~] = Cell_neighbors(cpos(1),cpos(2),9.5,dish,tot,1) ; 
                    [densecheck,~] = Cell_neighbors(cpos(1),cpos(2),9,dish,tot,1) ; 
                    [check3,opencoors3] = Cell_neighbors(cpos(1),cpos(2),3,dish,tot,8) ; 
                    
                    if check == 0
                        if densecheck == 0 
                            currspeed = currspeed - paramsBT3.densePen ; 
                        elseif loosecheck == 0
                            currspeed = currspeed - paramsBT3.loosePen ; 
                        else 
                            currspeed = currspeed - 1 ; 
                        end
                        
                        pickone = randi(size(opencoors,1)) ; % randomly pick a row
                        choiceCoors(1,:) = opencoors(pickone,:) ;  
                        [dish, dish_labels] = Adjust(cn,dish,dish_labels, 0, cpos(1),cpos(2)) ; %get out of old space
                        agent{cn}.position = [ choiceCoors(1,1) choiceCoors(1,2) ] ; %move into random nearby open space to the left
                        cpos = agent{cn}.position ; %rewrite cpos
                        curr = agent{cn} ; % reset curr 
                        [dish, dish_labels] = Adjust(cn,dish,dish_labels, 3, cpos(1),cpos(2)) ; 
                    
                    elseif size(opencoors3,1) == 8 
                        break 
                    
                    else 
                        if densecheck == 0 
                            currspeed = currspeed - paramsBT3.densePen ; 
                        elseif loosecheck == 0
                            currspeed = currspeed - paramsBT3.loosePen ; 
                        else 
                            currspeed = currspeed - 1 ; 
                        end
                        
                        % dish_labelscheck = dish_labels ; 
                        if cpos(2) >=2 
                            pos1 = cpos ; 
                            pos2 = [cpos(1) cpos(2)-1] ; 

                            agent1 = dish_labels(pos1(1,1),pos1(1,2)) ;
                            agent2 = dish_labels(pos2(1,1),pos2(1,2))  ;
                            
                            if agent2 ~= 0 % only swap with cells

                                dish_labels(pos1(1,1),pos1(1,2)) = agent2 ;
                                dish_labels(pos2(1,1),pos2(1,2)) = agent1 ;

                                agent{agent1}.position = [pos2(1,1) pos2(1,2)] ;
                                agent{agent2}.position = [pos1(1,1) pos1(1,2)] ;

                                dummysave = dish(pos1(1,1),pos1(1,2)) ;
                                dish(pos1(1,1),pos1(1,2)) = dish(pos2(1,1),pos2(1,2)) ;
                                dish(pos2(1,1),pos2(1,2)) = dummysave ;
                            end 

                            cpos = agent{cn}.position ; %rewrite cpos
                            curr = agent{cn} ; % reset curr

                        end 
                    end 
                end

                % Replicate if its been enough time 
                if tick-clast_r >= crep_rate
                    [isnew, newcoors] = Cell_neighbors(cpos(1),cpos(2),0,dish,tot,1) ; % Looking for empty space to go into
                    if isnew == 0 % space to replicate! 
                        actrep = tick-clast_r ; 
                        BT3Reps = vertcat(BT3Reps, actrep) ; 
                        crep_rate = crep_rate - paramsBT1.parentpoleRep ; % parental poles for OG cell ; 
                        clast_r = tick ; % rewrite last rep
                        agnts_new = agnts_new + 1 ; 
                        split = cbiomass/2 ; 
                        cbiomass = split ; 
                        posnew = [newcoors(1) newcoors(2)] ; 
                        bt3= bt3 + 1 ;
                        agent{prev_agnts+agnts_new} = BT3(paramsBT3.rep_rate,tick,cbiomass, posnew, 1) ; 
                        [dish, dish_labels] = Adjust(prev_agnts+agnts_new,dish,dish_labels, 3, newcoors(1),newcoors(2)) ; 
                    end 
                end

                % End Matters 
                agent{cn} = BT3(crep_rate,clast_r, cbiomass, cpos, 1) ;
            else 
                continue
            end
        end 
    end
%% Final Iteration Counts

    Ncounts(1) = nnz(dish(:) == 9) ; 
    Ncounts(2) = nnz(dish(:) == 9.5) ; 
    DenseMucusThickness(tick) = Ncounts(1)/tot ; 
    LooseMucusThickness(tick) = Ncounts(2)/tot  ; 

    % update datasets with current status for tick 
    BT1_TOTPOP(tick) = bt1 ; 
    BT2_TOTPOP(tick) = bt2; 
    BT3_TOTPOP(tick) = bt3 ; 
    DISH(tick).dishvals = dish ; 

    for metcount = 1:1:mets
        METABOLITES(metcount,tick).conc = metabolites(mets).conc ; 
    end 

    if flag_PrePath == 1 && bt3 == 0
        AllPathDead = tick ;
        flag_PrePath = 0 ;
        antipathDegrade = 1 ; 
    end

if rem(tick,inter) == 0 || tick == 1 
   
    %% Plot  
    W = 1 ;
    H = tot ; 

    figure(dish_graph)
    set(gcf,'color', 'white')
    sgtitle(DishgraphTitles)
    
    subplot(figrows,2,plotcount)
    hold on
    myplot(dish(W:H,W:H) == 1, 'b') 
    myplot(dish(W:H,W:H) == 10, 'r') 
    myplot(dish(W:H,W:H) == 2, 'c')  
    myplot(dish(W:H,W:H) == 9, [0.17 0.37 0.14]) 
    myplot(dish(W:H,W:H) == 9.5, [0.41 0.88 0.33])  
    myplot(dish(W:H,W:H) == 3, [0.4940 0.1840 0.5560 ]) 

    title(tick)
    hold off
   
    plotcount = plotcount + 1 ; 
    
end 

end 

%% Plotting

DEADPATH(countinput) = AllPathDead ; 

t = iterationrate*(1:1:nsteps) ;
halfmet = mets/2 ; 
METABOLITES1 = METABOLITES(1:halfmet) ; 
METABOLITES2 = METABOLITES(halfmet+1:mets) ; 

if saveit == 1  
    save(workspacetitle,"DISH", "BT1_TOTPOP","BT2_TOTPOP","BT3_TOTPOP","DenseMucusThickness","LooseMucusThickness","AllPathDead", "DEADPATH", "METABOLITES1","METABOLITES2")
    dishfilename = fullfile(parentfoldername, folder,figname) ;
    saveas(figure(dish_graph), dishfilename) ;
end 

figure(final_graph)
set(gcf,'color', 'white')
subplot(3,2,1)
hold on 
plot(t,BT1_TOTPOP, 'DisplayName', linelabels, 'Color',[myPlotColors(plotcolornum)],'LineWidth',2) 
title('a. BT1')
set(gca,'TitleHorizontalAlignment', 'left','FontSize',15) ; 
xlabel('time (min)')
ylabel('Concentration (units)')
lgd1 = legend ; 
lgd1.Title.String = legendtitles;

figure(final_graph)
subplot(3,2,2)
hold on 
plot(t, BT2_TOTPOP, 'DisplayName', linelabels, 'Color', [myPlotColors(plotcolornum)],"LineWidth",2)
title('b. BT2')
set(gca,'TitleHorizontalAlignment', 'left','FontSize',15) ; 
xlabel('time (min)')
ylabel('Concentration (units)')
lgd2 = legend ; 
lgd2.Title.String = legendtitles; 


figure(final_graph)
subplot(3,2,3)
hold on 
plot(t, BT3_TOTPOP, 'DisplayName', linelabels ,'Color', [myPlotColors(plotcolornum)],"LineWidth",2)
title('c. BT3')
set(gca,'TitleHorizontalAlignment', 'left','FontSize',15) ; 
xlabel('time (min)')
ylabel('Concentration (units)')
lgd2 = legend ; 
lgd2.Title.String = legendtitles; 

figure(final_graph)
subplot(3,2,5)
hold on 
plot(t, DenseMucusThickness, 'DisplayName', linelabels ,'Color', [myPlotColors(plotcolornum)],"LineWidth",2)
title('d. Dense Mucus Thickness')
set(gca,'TitleHorizontalAlignment', 'left','FontSize',15) ; 
xlabel('time (min)')
ylabel('Avg Thickness (\mum)')
lgd2 = legend ; 
lgd2.Title.String = legendtitles; 
ylim([0 100])
sgtitle(Finalgraphtitle)

figure(final_graph)
subplot(3,2,6)
hold on 
plot(t, LooseMucusThickness,':', 'DisplayName', linelabels ,'Color', [myPlotColors(plotcolornum)],"LineWidth",2)
title('e. Loose Mucus Thickness')
set(gca,'TitleHorizontalAlignment', 'left','FontSize',15) ; 
xlabel('time (min)')
ylabel('Avg Thickness (\mum)')
lgd2 = legend ; 
lgd2.Title.String = legendtitles; 
ylim([0 100])
sgtitle(Finalgraphtitle)
%%
if saveit == 1 
    finalgraphsname  = fullfile(parentfoldername, folder ,'PopulationPlots.fig') ;
    saveas(figure(final_graph), finalgraphsname) ;
end 
%% 

% Output final death counts
BT1_die
BT2_die
BT3_die        

if AllPathDead == 0 
    disp("Pathogens Surviving: " + num2str(bt3))
else 
    disp("Pathogens Died at tick = " + num2str(AllPathDead)) 
end 

disp("BT1 avg Reps: " + num2str(mean(BT1Reps)))
disp("BT2 avg Reps: " + num2str(mean(BT2Reps)))
disp("BT3 avg Reps: " + num2str(mean(BT3Reps)))
disp("Goblet avg Release: " + num2str(mean(GC_MUCProd)))

end 