#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Description:
    
This tool is a set of functions that simulate the genetic inheritance of traits 
(represented as single nucleotide polymorphisms, or SNPs) in a population of individuals.
The simulation is based on the principles of Mendelian genetics, with parents passing on randomly
selected alleles to their offspring. 
The program produces an initial population with plink based on a simulation parameter file, 
randomly selects breeding pairs, produces offsprings, and iterates over multiple generations 
to produce a new population. The final population is output in a standard format for genetic analysis,
the PED format, which contains information on each individual's ID, sex, and genotype for each SNP.
Finally, it runs frequency analysis with plink on both the initial population and the last population
and plots the MAF for each SNP for both populations. This shows how MAF have changed over generations.

Functions:

1- The required packages are imported.
2- There are five functions defined in this code.
    - add_sex(pop): takes a population and changes the sex of the individuals to males and females 
    with a uniform distribution of P(m)=P(f)=0.5.
    - offspring(pop,genTag,ID, indM, indF) takes a population as a dataframe, new generation tag,   
    and IDs of the parents, and produces one offspring. The offspring's sex is randomly 
    selected to be male or female.
    - select_breeders(pop, breedPer) takes a population and selects males and females to breed based on 
    a float value between 0 and 1, breedPer, which represents the percentage of the population 
    that will breed. The number of breeding events (couples) is set by the percentage of the minimum
    of sex counts.
    - produce_generation(pop,breedPer,genTag,maxOffs,numGenerations) takes a population, breeding percentageentage,
    maximum number of offspring a couple can have, and the number of generations to be produced. 
    The function produces new generations based on the breeding percentage and the maximum number 
    of offspring a couple can have, and the total number of new generations is given by 
    numGenerations. 
    - out_files(pedDFrame,pedName, mapName) takes a dataframe of the final population and outputs it in 
    two files. The PED file contains information on each individual's ID, sex, and genotype for each 
    SNP. The MAP file contains information on the chromosomal position of each SNP. Also, it uses 
    Plink to perform frequency analysis and produce a statistics file.

3- GUI window is defined and initialized with lables and entries.
4- GUI Functions:
    - Function runSimulation(): used as the function of the button 'Run' to perform the simulations
    and plot the results. The function execute plink bash commands to run the simulation. Then, it 
    calls the produce_generation function to generate the new offspring population. The number
    of possible offsprings is given at 3. It can be changed by editing line 260 and replacing 3 with
    the desired maximum.
    - Function: select_file(): A function for the 'Select . sim file' to capture the path of the 
    simulation parameters file.

Useage:
Run the script 'sim_popgen' with an environment with python installation (v 3.9.13) with bash command:
         'python sim_popgen.py'
1- The tool uses GUI to prompt the user for the input parameters and .sim (plink parameters file).
2- The user clicks on the Run button and the simulation starts.
3- A plot is generated in the canvas with MAF for every SNP for both the initial population and the 
    last possible generation.
4- The output file are written in the same directory. PED, LOG, MAP and FRQ files can be found for both
    generations.


@author: Tamim AlMurad with ChatGPT support.
"""

#############################################################################################
######################################## Packages ###########################################
#############################################################################################


from tkinter import *
from tkinter import filedialog
from tkinter import font
import os as bash
import pandas as pd
import random
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox

#############################################################################################
######################################## Functions###########################################
#############################################################################################



def add_sex(pop):
    '''
    Funcrion add_sex(pop) : Change the sex of the simulated population (usually all females)
    to males and females with uniform distribution of P(m)=P(f)=0.5. 
    '''
    #Create list of random 1s and 2s (males and females) with size of the population (pop).
    sex = random.choices([1,2],k=pop.shape[0])
    #Change the "sex" column in the a population to have both sexes.
    pop[4]=sex
    return 
 
def offspring(pop,genTag,ID, indM, indF):
    '''
    Function offspring(pop,genTag, indM, indF): Produces one offspring given
    a population as a datafraem, new generation tag and IDs of the parents.
    The offspring will be randomly selected to be male or female.
    '''
    #Create a list for the offspring individual filling the first 6 columns of 
    #its entry in a ped format.
    oS = ['offSpring'+genTag,'os'+str(ID),'0','0',random.randint(1,2),'N']
    #Get the number of SNPs.
    numSNPs = int((pop.shape[1]-6)/2)
    
    #Get random allel position from parents.
    randomAllelM = random.randint(0,1)
    randomAllelF = random.randint(0,1)
    
    #Get the indices of the parents from their IDs.
    maleIndex = pop.loc[pop[1]==indM].index[0]
    femaleIndex = pop.loc[pop[1]==indF].index[0]

    #Get the alleles for every SNP for the offspring individual from parents.
    for i in range(6,5+numSNPs*2,2):
        oS.append(pop[i+randomAllelM][maleIndex])
        oS.append(pop[i+randomAllelF][femaleIndex])
        
    return oS

def select_breeders(pop, breedPer):
    '''
    Function select_breeders(pop, breedPer): out of a population, it selects male and females 
    to breed. The breedPer is a float between 0 and 1 that is the value of the percentage of the 
    population that will breed.
    The number of breed events (couples) is set by the percentage of the minimum of sex counts.
    '''
    #Get the number of couples to breed from the percentage and the population.
    maleCount = (pop[4]==1).sum()
    femaleCount = (pop[4]==2).sum()
    breedsCount = int(min(maleCount,femaleCount)*breedPer)
    
    #In two list, get the IDs of the breeding males and females.
    breedMales = random.sample(list(pop[pop[4]==1][1]),breedsCount)
    breedFemales = random.sample(list(pop[pop[4]==2][1]),breedsCount) 
    
    return breedMales, breedFemales

def produce_generation (pop,breedPer,genTag,maxOffs,numGenerations):
    '''
    Function produce_generation (pop,breedPer): Out of a population, produce new generation and 
    subsequenct generations based on the breeding percentage provided in breedPer and the maximum 
    number of offsprings a couple can have maxOffs. The total number of new generations is given 
    by numGenerations.The number of springs will be random between 1 and maxOffs.

    '''
    
    #Iterate to produce multiple generations.
    for k in range(0,numGenerations):
        #Get the IDs of the breeders.
        bMales, bFemales = select_breeders(pop, breedPer)
        
        #Initialize a list to hold a generation data.
        generation=[]
        
        #To show which generation is being produced.
        print('################### Starting generation number'+str(k)+"#######################")
        
        #Check if there is still males and females to breed.
        if (len(bMales)>0 and len(bFemales)>0):
        #Nested loop to iterate over the breeds, produce offsprings and add it to the generation list.
            for i in range(0,len(bMales)):
                for j in range(0,random.randint(1,maxOffs)):
                    generation.append(offspring(pop,genTag+str(k+1),str(k+1)+str(i)+str(j),bMales[i],bFemales[i]))
                
            
            #Make a dataframe out of the list.
            generationDF = pd.DataFrame(generation)
            #To produce the next generation, the new generation becomes the population.
            pop=generationDF
        
        #Not enough males and females to breed. Stop and output the last generation and generation number.
        else:
            print('The last Generation is '+str(k+1)+'th generation.')
            return generationDF, k+1  
    #Only the last generation is returned with the generation number starting from 0 for the original population.
    return generationDF, k+1    

def out_files(pedDFrame,pedName, mapName):
    '''
    Function out_files(pedDFrame,pedName, mapName): takes a dataframe representing generation SNP data
    , the file name prefix and .map file name (same map file generated in the simulation). As an output
    it produces .ped file, copies the .map file with a name as per the prefix and run allele frequency
    stats with plink producing .frq file.
    '''
    #Write a generation dataframe into ped file.
    pedDFrame.to_csv(pedName+".ped",index=False,sep=' ',header=None)
    #Copy the map file with a new name.
    bash.system("cp "+ mapName+".map" + " " + pedName+".map")
    #Produce statistics of the new generation with plink.
    bash.system("plink --file " + pedName +" --freq --out " + pedName)


#############################################################################################
######################################## GUI Part ###########################################
#############################################################################################

# Tkinter object.
root = Tk()
root.title("BreedSNP")
root.geometry("900x600")

#Lables and Entries
nCasesLab = Label(root, text='Enter the number of cases',font=("Arial",14))
nCasesLab.pack(anchor=W)
nCasesEnt = Entry(root,width=10)
nCasesEnt.pack(anchor=W,padx=50)

nControlLab = Label(root, text='Enter the number of controls',font=("Arial",14)).pack(anchor=W)
nControlEnt = Entry(root,width=10)
nControlEnt.pack(anchor=W,padx=50)

outputFileLab = Label(root, text='Enter the Prefix of the output file',font=("Arial",14)).pack(anchor=W)
outputFileEnt = Entry(root,width=10)
outputFileEnt.pack(anchor=W,padx=50)

genLab = Label(root, text='Enter the Number of Gernerations',font=("Arial",14)).pack(anchor=W)
genEnt = Entry(root,width=10)
genEnt.pack(anchor=W,padx=50)

breedPerLab = Label(root, text='Enter breeding percentage (value must be between 0 and 1)',font=("Arial",14)).pack(anchor=W)
breedPerEnt = Entry(root,width=10)
breedPerEnt.pack(anchor=W,padx=50)

#To store the path of the file chosen.
filePath = ''

# Run button.
def runSimulation():
    
    '''
    Function runSimulation(): used as the function of the button 'Run' to perform the simulations
    and plot the results.
    '''
    # get inputs from GUI and validate them.
    try:
       
        nCases = int(nCasesEnt.get())
        nControls = int(nControlEnt.get())
        outFileName = outputFileEnt.get()
        numGenerations = int(genEnt.get())
        breedPer = float(breedPerEnt.get())
        
    except ValueError:
        messagebox.showerror("Error", "Invalid input. Please enter valid numbers.")
        
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

    #Simulate with plink.   
    bash.system("plink --simulate " + filePath + " acgt --recode --out " + outFileName + " --simulate-ncases "+ 
              str(nCases) +" --simulate-ncontrols " +  str(nControls))
    #Get statistics with plink for the simulated population.
    bash.system("plink --file " + outFileName + " --freq --out " + outFileName)
    
    #Read the result ped file into a dataframe.
    data = pd.read_csv(outFileName+'.ped',sep=' ',header=None)
    #Change the sex of the simulated population to contain both makes and females.
    add_sex(data)

    #Produce new generations and return the last generation in a datafrome with the actual resulting number of generations.
    generation,numGenActuall = produce_generation(data, float(breedPer), 'G',3,int(numGenerations))
    
    #Write ped file named as offSFileName and produce a file with the statistics.
    offSFileName = 'Generation'+str(numGenActuall)
    out_files(generation, offSFileName,outFileName )

    #Read the frequency statistics files of the original population and the last generation.
    gen1Freq=pd.read_fwf(outFileName+'.frq')
    gen2Freq=pd.read_fwf(offSFileName+'.frq')
    
    #Label to indicate the simulation is done.
    doneLabe = Label(root, text='Simulation is done',font=("Arial",14)).pack()
    
    #Plot the results in the canvas of the GUI.
    fig = plt.Figure(figsize=(9, 4), dpi=100)
    gen1Freq.plot(x='SNP',y='MAF', label='Generation 0', ax=fig.gca(),marker='o')
    gen2Freq.plot(x='SNP',y='MAF', label='Generation '+str(numGenerations), ax=fig.gca(),marker='o')
    fig.gca().set_ylabel('MAF')
    plt.legend()
   
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack()
       
    
def select_file():
    
    '''
    Function: select_file(): A function for the 'Select . sim file' to capture the path of the simulation
    parameters file.
    '''
    global filePath 
    filePath = filedialog.askopenfilename()
    fileLab = Label(root, text='The selected file is: '+filePath,font=("Times",10)).pack()
    #Allow the 'Run' file to be clicked.
    runButton.configure(state="normal")

# Create a button to open the file dialog.
fileButton = Button(root, text="Select a .sim File",font=("Arial",11),command=select_file,height = 5, width=18)
fileButton.place(x=700, y=10)

# Button to start the simulation.
runButton = Button(root, text="Run",command=runSimulation,font=("Arial",11),height= 5, width=18)
runButton.place(x=700, y=100)
#No running without a file.
runButton["state"]= DISABLED

root.mainloop()







