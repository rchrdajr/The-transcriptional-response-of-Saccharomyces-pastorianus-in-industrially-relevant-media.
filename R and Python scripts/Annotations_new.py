# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:52:51 2024

@author: z38709rr
"""

## import libraries
import os
import csv






## LOAD ETH INTERSECTION DATA
#%%
## change directory
os.getcwd()
os.chdir("\\\\nask.man.ac.uk\\home$\\Documents\\R_Proj\\1_Feature_counts\\unnorm_data_1")
#%%
##---------------------------------

## CREATE DICTIOMARY FROM ANNOTATIONS
## LOAD ANOTATIONS DATA
#%%

new_dico = {} # to define dictionary
with open("uniprot_annotation.csv", "r") as uniprot_annotation:
    header = next(uniprot_annotation)
    for line in uniprot_annotation:
        print(line)
        l = line.split("\n")[0] ## removes "\n" from the end of line, that creates a list, with what is before the "\n", and what's after, so we chose the first element
        l = l.replace('"','')
        l = l.split("\t") # creates a list with everything that was separated by "\t", in our case it's the columns of the file 
        print(l[0])
        print(l[1])
        parent = l[0]
        geneID = l[1]
        protID = l[2]
        entry = l[3]
        entryName = l[4]
        GeneName = l[5]
        primaryGname = l[6]
        synonymGname = l[7]
        orderedLocusGname = l[8]
        orfGname = l[9]
        protName = l[10]
        status = l[11]
        organism = l[12]
        length = l[13]
        ECnumber = l[14]
        GObioProcess = l[15]
        GOmolecFunction = l[16]
        GOcellComponent = l[17]
        function = l[18]
        catalyticActivity = l[19]
        absorb = l[20]
        cofact = l[21]
        enzRegulation = l[22]
        kinetics = l[23]
        pathw = l[24]
        redoxPot = l[25]
        tempDepend = l[26]
        pHDepend = l[27]
        subcellLoc = l[28]
        intramemb = l[29]
        transmemb = l[30]
        subUStruc = l[31]
        interact = l[32]
        tissueSpecif = l[33]
        
        new_dico[geneID] = [parent, geneID, protID, entry, entryName, GeneName, primaryGname, synonymGname, orderedLocusGname,
                orfGname, protName, status, organism, length, ECnumber, GObioProcess, GOmolecFunction, GOcellComponent,
                function, catalyticActivity, absorb, cofact, enzRegulation, kinetics, pathw, redoxPot, tempDepend,
                pHDepend, subcellLoc, intramemb, transmemb, subUStruc, interact, tissueSpecif] ## exerice check how to add values ! # ALMOST THERE, MISSING SOMETHING
     
#%%
## example of what you want: expression_dico13eth["SPGP00000000"] = ["-1.5\t0.001"]

list_conditions = ["eth","leu","wort"]
temp_conditions = ["13", "22", "30"]
for tempcond in temp_conditions:
    for cond in list_conditions:
        ## create empty dictionaries for each condiion created in these loops
        globals()["expression_dico_"+str(tempcond)+str(cond)] = {}
        with open("analysis\\fc2_aa_vs_"+str(cond)+"\\subset_fc2_"+str(tempcond)+"aa_vs_"+str(tempcond)+str(cond)+".csv", "r") as data: #open file
            next(data)
            for line in data:
                line = line.split("\n")[0]
                line = line.replace('"','')
                line = line.replace("SPGP_","SPGP")
                line = line.split(",")
                gene = line[0]
                infos = str(line[2])+"\t"+str(line[6]) ### "LFC\tpadjValue"
                globals()["expression_dico_"+str(tempcond)+str(cond)][gene] = infos 
        
        
##### SOOKIE
        
        
#%%

annotation_dico = {} # to define dictionary
with open("uniprot_annotation.csv", "r") as uniprot_annotation:
    header = next(uniprot_annotation)
    for line in uniprot_annotation:
        l = line.replace('"','')
        l = l.split("\t") # creates a list with everything that was separated by "\t", in our case it's the columns of the file 
        geneID = l[1]
        annotation_dico[geneID] = line
        

list_conditions = ["eth","leu","wort"]

###

for cond in list_conditions:
    with open("analysis\\fc2_aa_vs_"+str(cond)+"\\intersect_13aa_13"+str(cond)+"_vs_22aa_22"+str(cond)+"_vs_30aa_30"+str(cond)+".csv", "r") as data, open("analysis\\fc2_aa_vs_"+str(cond)+"\\annotated_intersect_13aa_13"+str(cond)+"_vs_22aa_22"+str(cond)+"_vs_30aa_30"+str(cond)+".csv","w") as output: #open file
        next(data)
        newheader = "LFC 13C"+"\t"+"padj-value 13C"+"\t"+"LFC 22C"+"\t"+"padj-value 22C"+"\t"+"LFC 30C"+"\t"+"padj-value 30C"+"\t"+str(header)
        output.write(newheader)
        for line in data: # loop through lines
            line = line.split("\n")[0]
            line = line.replace('"','')
            line = line.replace("SPGP_","SPGP")
            line = line.split(",") # split lines into 2 columns by the comma
            gene = line[1]
            expr = globals()["expression_dico_13"+str(cond)][gene]+"\t"+globals()["expression_dico_22"+str(cond)][gene]+"\t"+globals()["expression_dico_30"+str(cond)][gene]
            if gene in annotation_dico.keys():
                newLine = annotation_dico[gene]
                output.write(str(expr)+"\t"+str(newLine)) ### to do: add the expression and padj value in first columns ( output.write("LFC 13C"+"\t"+"padj-value 13C"+"\t"+"LFC 22C"+"\t"+"padj-value 22C"+"\t"+"LFC 30C"+"\t"+"padj-value 30C"+"\t"+sr(newLine)))
            else:
                output.write(str(expr)+"\t"+"S. pastorianus\t"+str(gene)+"\n")
   
## LOAD OTHOLOGIES DATA
#%%
with open("orthologies_FM1318orthoCBS1513_Scerevisiae_clean.csv", "r") as data: # save as dico key = S.cerevisiae protid, value = ortholog protid
    next(data)
    liness = data.readlines()
orth_dico = {line.split(",")[0]: line.split(",")[1] for line in liness} #split into 2 values + save into dico (:???????)

for cond in list_conditions: # replacing values with orthologs
    with open("analysis\\fc2_aa_vs_"+str(cond)+"\\annotated_intersect_13aa_13"+str(cond)+"_vs_22aa_22"+str(cond)+"_vs_30aa_30"+str(cond)+".csv", "r") as intersect_data, open("analysis\\fc2_aa_vs_"+str(cond)+"\\orth_annotated_intersect_13aa_13"+str(cond)+"_vs_22aa_22"+str(cond)+"_vs_30aa_30"+str(cond)+".csv", "w") as orth_intersect_data:
        #next(intersect_data)
        lines = intersect_data.readlines()
        #new_orth_line = [] 
        
        for line in lines:
            #split line based on tab characters
           columns = line.strip().split("\t")
           if len(columns) < 9:
               # line is unchanged
               pass
           elif columns[8] in orth_dico:
               columns[8] = orth_dico[columns[8]].strip() # update column strips \n
          
           newlineval = "\t".join(columns)+ "\n"
           orth_intersect_data.writelines(newlineval)