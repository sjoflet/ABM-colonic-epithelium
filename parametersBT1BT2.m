function [paramsBT1, paramsBT2, paramsGC, paramsBT3] = parametersBT1BT2(itertime)

%% parameters

% iteration time based on 1 = 1 min 

%% Faculative BT1
paramsBT1.rep_rate = 35/itertime ; % want 40 min, tick rate based on iteration time
paramsBT1.O2 = 0 ; % min amount to start aerobic resp
paramsBT1.O2consume = 0.1489*itertime ; % ug per cell per min (BB 29) 
paramsBT1.starve = -60*48*itertime; % how long it can stay dormant
paramsBT1.galtiKm = 6.0117e-4 ; % [g/L] galactitol Km 
paramsBT1.galtiVmax = (13.4)*itertime ; % Vmax per min * min per tick
paramsBT1.galacKm = 6.0117 ; % galactose Km
paramsBT1.galacVmax = (5)*itertime ; % vmax per min * min per tick
paramsBT1.biomassLost = 1 ; 
paramsBT1.noAnchorFood = 0.05*itertime ; 
paramsBT1.degradationOdds = 1*itertime ; 
paramsBT1.repMin = paramsBT1.rep_rate ; 
paramsBT1.antibio =  50 ; % concentration of antibio that bacteria can withstand
paramsBT1.antipath =  4; % amount of antipathogen released
paramsBT1.proteo = .15 ;       
paramsBT1.MAMPs = 1 ; 
paramsBT1.PAMPs = 1 ; 

% Replication modifier factors 
paramsBT1.galtiRep = 1*itertime ; %food = quicker
paramsBT1.galacRep = 2*itertime ; % galactose = slower
paramsBT1.galacRepMax = 10*itertime+paramsBT1.rep_rate ; 
paramsBT1.starveRep = 10*itertime ; % slower with no food
paramsBT1.O2Rep = 1 ; % faster replication when aerobic resp 
paramsBT1.AnchorRep = -10*itertime ; % no anchor = slows
paramsBT1.AnchorRepMax = paramsBT1.rep_rate*2 ; % hard stop on slowing without anchoring
paramsBT1.parentpoleRep = -1 ; % older DNA = slows


%% Anaerobic BT2 
paramsBT2.rep_rate = 15/itertime ; % want 20 min, tick rate based on iteration time 
paramsBT2.O2 = 1e-15 ; % max amount tolerable
paramsBT2.but = 2*itertime ; % butyrate production rate 
paramsBT2.galti = 10*itertime ; % galactitol production rate 
paramsBT2.galacKm = 7.0117e-6 ; % galactitol Km 
paramsBT2.galacVmax = (20)*itertime ; % Vmax per min * min per tick
paramsBT2.dfKm = 4e-4 ; % dietary fiber Km
paramsBT2.dfVmax = .010 ; % dietary fiber Vmax
paramsBT2.biomassLost = 10 ; 
paramsBT2.starve = -60*24*itertime ; 
paramsBT2.degradationOdds = 8  ; %%%%
paramsBT2.repMin = 20/itertime ; 
paramsBT2.proteo  =  .13 ; 
paramsBT2.MAMPs = 1 ; 

% Replication modifier factors 
paramsBT2.galacRep = 1 ; 
paramsBT2.parentpoleRep = -1 ; 
paramsBT2.antibio = 4.5 ; 
paramsBT2.antipath = 2 ; 
paramsBT2.fiberdegrade = 5000; 
paramsBT2.fiberMin = 5 ; 
%% Goblet Cells 

paramsGC.butKm = 4 ; 
paramsGC.butVmax = 20 ;  
paramsGC.butMin = 60 ;
paramsGC.butScale = 10 ; 
paramsGC.noBut = 60 ; 
paramsGC.MucusBox = 5000*itertime ; 
paramsGC.O2 = 32*(32.5)*10^(-15)*itertime ; 
paramsGC.O2rec = 10e-15 ; 
paramsGC.mucMin = 50*itertime ; % minimum mucus production rate
paramsGC.mucProd = 100*itertime ; % standard mucus production rate

paramsGC.proteo = 4 ; 
paramsGC.antibioScale = 1 ; 
paramsGC.antibioReg = 50 ; 
paramsGC.epiMin = 5 ; 
paramsGC.antipath = 10 ; 
paramsGC.IgAMax = 100 ; 
paramsGC.IgAMin = 30 ; 
paramsGC.IgAscale = 5 ; 
paramsGC.igaMAMPs = 2 ; 
paramsGC.NitReg = 2 ; 

%% Pathogenic BT3

paramsBT3.repMin = 30/itertime ; 
paramsBT3.O2consume = 10*itertime;
paramsBT3.O2Rep = 1 ; 
paramsBT3.galtiRep = 2*itertime ; 
paramsBT3.galtiKm = 6.0117e-4 ; % galactitol Km 
paramsBT3.galtiVmax = (13.4)*itertime ;
paramsBT3.galacKm = 6.0117 ; % galactose Km
paramsBT3.galacVmax = (5.4)*itertime ; % vmax per min * min per tick
paramsBT3.biomassLost = 1 ; 
paramsBT3.speed = 6*itertime ; 
paramsBT3.rep_rate = 6/itertime ; 
paramsBT3.antipath  = 5 ; % max [anti-pathogens]
paramsBT3.antibio = 0.3 ;
paramsBT3.proteo = 1 ; 
paramsBT3.starve = -50/itertime ; 
paramsBT3.galtiRep  = 20/itertime ; 
paramsBT3.nitrateRep = 30/itertime ; 
paramsBT3.PAMPs = 5 ; 
paramsBT3.loosePen = 2 ;
paramsBT3.densePen = 5 ; 


paramsBT3.densemucDegr = 100 ; 
paramsBT3.loosemucDegr = 50 ; 
end 