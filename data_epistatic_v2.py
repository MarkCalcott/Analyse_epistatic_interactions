#############################################################################
#
# Author: Mark Calcott
# Date: 20190402
# Overview: Analyse epistatic interactions
#
#############################################################################

import numpy as np
from scipy import stats

def create_IC50_dictionary_raw_data(filename):
    '''
    Takes in a file containing a column of sequences and a column of IC50's
    Returns a dictionary of sequence:IC50 combinations
    '''
    dictIC50 = {}
    fh = open(filename)
    assert(fh.readline().startswith("Sequence"))
    for line in fh:
        if line == "": break
        data = line.strip().split(',')
        sequence = data[0].split('_')[1] #removes the numbering for each sequence
        dictIC50[sequence] = [float(x) for x in data[1:] if not x ==""]
    fh.close()
    return dictIC50.copy()

def wrapper(filename, WT):
    '''
    For each mutation, calculates whether there is a statistically significant change in activity.
    '''

    originalData = create_IC50_dictionary_raw_data(filename)

    #Create a data structure for analysis
    sortedData = {"WT":{}, "Mt":{}}
    for i in range(len(originalData.keys()[0])):
        #Add all variants at position 1 to the corresponding WT or mutant dictionary
        ## Do the same for variants at position 2, 3, etc...
        sortedData["WT"][i] = {}
        sortedData["Mt"][i] = {}
        for key in originalData.keys():
            #Creates a key with the current position marked with a '*'
            mutant_key = key[:i] + "*" + key[i+1:]
            if key[i] == WT[i]:
                #If the '*' residue is identical to the wildtype sequence, add to a WT dictionary
                sortedData["WT"][i][mutant_key] = originalData[key]
            else:
                #Otherwise add to the mutant dictionary
                sortedData["Mt"][i][mutant_key] = originalData[key]

    #Save data
    saveFile = open("EpistaticData_raw.csv", 'w')
    saveFile.write("Mutation,Position,WT Ave,WT StDev,Mt Ave,Mt StDev,Ratio,Effect,pvalue\n")
    print "Position\tPositive\tNegative\tNull"
    for i in range(len(originalData.keys()[0])):
        count_pos = 0
        count_neg = 0
        count_neu = 0
        
        for key in sortedData["WT"][i].keys():
            #For each position in the sequence, compare the WT sequence to the corresponding mutant
            WT_raw = sortedData["WT"][i][key]
            WT_IC50 = np.average(WT_raw)
            WT_StDev = np.std(WT_raw, ddof=0) #standard deviation of a population
            Mt_raw = sortedData["Mt"][i][key]
            Mt_IC50 = np.average(Mt_raw)
            Mt_StDev = np.std(Mt_raw, ddof=0) #standard deviation of a population
            ratio = Mt_IC50/WT_IC50
            #Perform t-test and label positive or negative
            pvalue = stats.ttest_ind(WT_raw, Mt_raw, equal_var=False)[1] #unpaired, two-tailed t-test, unequal variances
            if pvalue > 0.05:
                effect = "Null"
                count_neu += 1  
            elif ratio >= 1:
                effect = "Positive"
                count_pos += 1
            elif ratio <= 1:
                effect = "Negative"
                count_neg += 1

            saveFile.write(",".join([str(x) for x in [key, i+1, WT_IC50, WT_StDev, Mt_IC50, Mt_StDev, ratio, effect, pvalue]]) + "\n")
        print "\t".join([str(x) for x in [i+1, count_pos, count_neg, count_neu]])
    saveFile.close()


wrapper("IC50dataerror.csv", 'SHTKSRF')



