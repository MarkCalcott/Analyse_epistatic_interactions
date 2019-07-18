#############################################################################
#
# Author: Mark Calcott
# Date: 20181207
# Overview: Analyse evolutionary pathways
#
#############################################################################

def create_IC50_dictionary(filename):
    '''
    Takes in a file containing a column of sequences and a column of IC50's
    Returns a dictionary of sequence:IC50 combinations
    '''

    dictIC50 = {}
    fh = open(filename)
    assert(fh.readline().startswith("Sequence"))
    for line in fh:
        data = line.strip().split(',')
        sequence = data[0].split('_')[1]
        dictIC50[sequence] = float(data[1])

    fh.close()
    return dictIC50.copy()

def stopcondition(mutantActivity, parentActivity, increaseTolerance):
    return mutantActivity <= increaseTolerance*parentActivity #stop if the activity of the mutant is less than or equal to the parent

def creategraph(start, end, dictIC50, increaseTolerance, checkIC50 = True, stoppingRule = stopcondition):
    '''
    Returns a list containing all the paths from
    the start to the end on the condition that each
    node is separated by 1 mutation and increases IC50
    '''
    graphs = []

    #Look at each residue, and if it differs from the end residue then add to a list
    newNodes = []
    for i in range(len(start)):
        if start[i] != end[i]:
            newstart = start[:]
            newstart[i] = end[i]
            
            ##If the checkIC50 and stop condition are met, then add this node as the end point
            stop = stopcondition(dictIC50[''.join(newstart)], dictIC50[''.join(start)], increaseTolerance)
            
            if checkIC50 and stop:
                pass
            else:
                newNodes.append(newstart[:])
    #Return current node if there are no new nodes
    if len(newNodes) == 0:
        return [[''.join(start)]]

    #Create new graphs with nodes that increase IC50
    for x in newNodes:
        newgraphs = creategraph(x, end, dictIC50, increaseTolerance, checkIC50, stoppingRule)
        for y in newgraphs:
            newList = []
            newList = [''.join(start)]
            newList.extend(y)
            graphs.append(newList[:])
    return graphs

def outputPathsAsCSV(paths, filename, dictIC50, length = 9):
    '''
    Creates an output file of paths and IC50 values
    '''
    outputfile = open(filename, 'w')
    for path in paths:
        IC50values = []
        currentPath = ['' for x in range(length)]
        for i in range(len(path)):
            currentPath[i] = path[i]
            IC50values.append(dictIC50[path[i]])
        line = currentPath[:]
        line.extend([str(x) for x in IC50values])
        outputfile.write(','.join(line) + '\n')
    outputfile.close()

def wrapper(filename, start, end, increaseTolerance = 1):

    dictIC50 = create_IC50_dictionary(filename)

    #Output all paths in csv format
    filename = '_'.join([start, end, str(increaseTolerance), 'allpaths.csv'])
    paths = creategraph(list(start), list(end), dictIC50, increaseTolerance, False)
    outputPathsAsCSV(paths, filename, dictIC50)

    #Output truncated paths in csv format
    filename = '_'.join([start, end, str(increaseTolerance), 'truncatedpaths.csv'])
    paths = creategraph(list(start), list(end), dictIC50, increaseTolerance)
    outputPathsAsCSV(paths, filename, dictIC50)
    
wrapper("IC50data.csv", 'SHTKSRF', 'YNYRYDH', 1.2)

