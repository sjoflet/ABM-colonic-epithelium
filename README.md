# ABM-colonic-epithelium
MATLAB written agent-based model of a cross-sectional niche of the colonic epithelium 

Samantha Fletcher, Carly Ching, Mark Suprenant, Darash Desai, and Muhammad Zaman. 
Boston University 2025

MAIN.m is the primary file to run ABM simulations. 
  Currently there are two test types: 
    1 - varying dietary fiber and 2 - varying pathogen seeding amounts
  Iteration actionables occur within this script

ADJUSTABLE SCRIPTS - parameters can be tuned to different model bacteria types or adjusted based on specific questions, and         additional components can be easily added
MetabolitesStruct.m contains small molecule diffusion and mucus-trapping parameters
parametersBT1BT2.m contains agent parameters for all agent types: BT1, BT2, BT3, and Goblet cells
myPlotColors.m is a list of MATLAB triplet code colors that are applied to simulation runs in input order

FUNCTIONAL SCRIPTS
  Objects: 
    BT1.m 
    BT2.m
    BT3.m 
    Goblet.m

  General Functions: 
    Adjust.m - swaps dish values and dish label values
    Cell_neighbors.m - examines surrounding environment
    Diffuse_TDMAvec.m - Vectorized Alternating-Direction Implicit method diffusion solved by Thomas Diagonal Matrix Algorithm
        with boundary condition options: Neumann or Dirichlet 
    MM.m - Michaelis-Menten solve function
    addpathogen.m - adds pathogens at the lumen end and displaces other agents where necessary
    mucus_environment.m - controls dense mucus expansion into loose mucus based on proteolytic factor levels
    myplot.m - creates scatter plot with filled dots for specific dish matrix values

    
    
