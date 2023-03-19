## Introduction

This program simulates a population with SNP data, for the purpose of use in a genetic association study, using the Python programming language. The user can set parameters including the number of individuals to simulate, number of generations to simulate, the probability of breeding. The program outputs a .ped and a .map file that can be used in subsequent genetic association analysis.

## Required packages
This program requires the following packages to be installed:

pandas

random

matplotlib

tkinter

Also, the program requires PLINK installations (PLINK v1.90b6.21) as it passes bash commands with bash version 5.0.17(1).

## Functions

    1- add_sex(pop): takes a population and changes the sex of the individuals to males and females 
    with a uniform distribution of P(m)=P(f)=0.5.

    2- offspring(pop,genTag,ID, indM, indF) takes a population as a dataframe, new generation tag,   
    and IDs of the parents, and produces one offspring. The offspring's sex is randomly 
    selected to be male or female.

    3- select_breeders(pop, breedPer) takes a population and selects males and females to breed based on 
    a float value between 0 and 1, breedPer, which represents the percentage of the population 
    that will breed. The number of breeding events (couples) is set by the percentage of the minimum
    of sex counts.

    4- produce_generation(pop,breedPer,genTag,maxOffs,numGenerations) takes a population, breeding percentageentage,
    maximum number of offspring a couple can have, and the number of generations to be produced. 
    The function produces new generations based on the breeding percentage and the maximum number 
    of offspring a couple can have, and the total number of new generations is given by 
    numGenerations. 

    5- out_files(pedDFrame,pedName, mapName) takes a dataframe of the final population and outputs it in 
    two files. The PED file contains information on each individual's ID, sex, and genotype for each 
    SNP. The MAP file contains information on the chromosomal position of each SNP. Also, it uses 
    Plink to perform frequency analysis and produce a statistics file.

## GUI Functions:

    1- Function runSimulation(): used as the function of the button 'Run' to perform the simulations
    and plot the results. The function execute plink bash commands to run the simulation. Then, it 
    calls the produce_generation function to generate the new offspring population. The number
    of possible offsprings is given at 3. It can be changed by editing line 260 and replacing 3 with
    the desired maximum.

    2- Function: select_file(): A function for the 'Select . sim file' to capture the path of the 
    simulation parameters file.

## Inputs
This tool requires the following inputs from the user:

1- Number of control individuals.

2- Number of cases individuals.

3- Number of generations to simulate.

4- Percentage of individuals that will breed, expressed as a float between 0 and 1.

5- Simulation file name prefix (for initial population).

6- Simulation parameter file (.sim) given in the plink format.
 https://zzz.bwh.harvard.edu/plink/simulate.shtml
 
## Outputs
This tool outputs the following files:

Two .ped files containing the simulated and last resulting generation individual SNPs.

Two .map files containing the SNPs of initial and last population.

Two .frq files containing MAF statistiscs for both generations.

Two .log files from running plink on both generations.

A plot of SNPs and MAF for both generations on the canvas.

## Useage:
Run the script 'sim_popgen' with an environment with python installation (v 3.9.13) with bash command:
         'python sim_popgen.py'

1- The tool uses GUI to prompt the user for the input parameters and .sim (plink parameters file).

2- The user clicks on the Run button and the simulation starts.

3- A plot is generated in the canvas with MAF for every SNP for both the initial population and the 
    last possible generation.
4- The output file are written in the same directory. PED, LOG, MAP and FRQ files can be found for both
    generations.



![My Image]Screenshot from 2023-03-19 13-11-16.png

