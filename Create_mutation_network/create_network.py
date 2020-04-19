#############################################################################
#
# Author: Mark Calcott
# Date: 20200419
# Overview: Analyse evolutionary pathways
#
#############################################################################

def create_EC50_dictionary(filename):
    '''
    Takes in a file containing a column of sequences and a column of EC50's
    Returns a dictionary of sequence:EC50 combinations
    '''
    dictEC50 = {}
    fh = open(filename)
    assert(fh.readline().startswith("Sequence"))
    for line in fh:
        data = line.strip().split(',')
        sequence = data[0].split('_')[1]
        dictEC50[sequence] = float(data[1])

    fh.close()
    return dictEC50.copy()


def create_network_dictionary(keys):
    '''
    Creates a dictionary of the number of mutations from wt as the key,
    and all the sequences with that number of differences as values
    '''
    dictNetwork = {}
    for i in range(len(keys[0])+1):
        dictNetwork[i] = []

    for j in keys:
        distance = j.count('1')
        newData = dictNetwork[distance]
        newData.append(j)
        dictNetwork[distance] = newData
            
    return dictNetwork.copy()


def hamming_distance(x, y):
    '''
    Returns the number of differences between two strings
    '''
    distance = len(x)
    for i in range(distance):
        if x[i] == y[i]:
            distance -=1
    return distance

def saveEdges(name, dictEC50, dictNetwork):
    '''
    Outputs file for the edges in cytoscape
    '''
    filename = name + "_edges.csv"
    outputfile = open(filename, 'w')
    
    keys = sorted(dictNetwork.keys()[:-1])
    
    for key in keys:
        for node1 in dictNetwork[key]:
            for node2 in dictNetwork[key+1]:
                if hamming_distance(node1,node2) == 1:
                    data = [node1, node2, dictEC50[node2]/float(dictEC50[node1])]
                    outputfile.write(','.join([str(x) for x in data]) + '\n')
    outputfile.close()

def saveNodes(name, dictEC50, dictNetwork):
    '''
    Outputs file for the nodes in cytoscape
    '''
    filename = name + "_nodes.csv"
    outputfile = open(filename, 'w')
    
    for key in dictEC50.keys():
        data = [key, round(dictEC50[key]/float(dictEC50['0000000']),1)]
        outputfile.write(','.join([str(x) for x in data]) + '\n')
    outputfile.close()

def stopcondition(mutantActivity, parentActivity, increaseTolerance = 1.16):
    return mutantActivity <= increaseTolerance*parentActivity #stop if the activity of the mutant is less than or equal to the parent


def findpaths(start, end, dictEC50, checkEC50 = True, stoppingRule = stopcondition):
    '''
    Returns a list containing all the paths from
    the start to the end on the condition that each
    node is separated by 1 mutation and increases EC50
    '''
    graphs = []

    #Look at each residue, and if it differs from the end residue then add to a list
    newNodes = []
    for i in range(len(start)):
        if start[i] != end[i]:
            newstart = start[:]
            newstart[i] = end[i]
            
            ##If the checkEC50 and stop condition are met, then add this node as the end point
            stop = stopcondition(dictEC50[''.join(newstart)], dictEC50[''.join(start)])
            
            if checkEC50 and stop:
                pass
            else:
                newNodes.append(newstart[:])
    #Return current node if there are no new nodes
    if len(newNodes) == 0:
        return [[''.join(start)]]

    #Create new graphs with nodes that increase EC50
    for x in newNodes:
        newgraphs = findpaths(x, end, dictEC50, checkEC50, stoppingRule)
        for y in newgraphs:
            newList = []
            newList = [''.join(start)]
            newList.extend(y)
            graphs.append(newList[:])
    return graphs

def outputPathsAsCSV(paths, name, dictEC50, length = 9):
    '''
    Creates an output file of paths and EC50 values
    '''
    filename = name + "_paths.csv"
    outputfile = open(filename, 'w')
    for path in paths:
        EC50values = []
        currentPath = ['' for x in range(length)]
        for i in range(len(path)):
            currentPath[i] = path[i]
            EC50values.append(dictEC50[path[i]])
        line = currentPath[:]
        line.extend([str(x) for x in EC50values])
        outputfile.write(','.join(line) + '\n')
    outputfile.close()

def wrapper(filename, start = '0000000', end = '1111111'):
    # Create dictionaries for network graph
    dictEC50 = create_EC50_dictionary(filename)
    dictNetwork = create_network_dictionary(dictEC50.keys())
    
    # Filename specific to data set
    name = filename.split()[0]
    
    #Save the edges of a graph for cytoscape
    saveEdges(name, dictEC50, dictNetwork)
    
    #Save the nodes of a graph for cytoscape
    saveNodes(name, dictEC50, dictNetwork)
    
    # Output all paths in whcih each step increases by 16%
    paths = findpaths(list(start), list(end), dictEC50)
    outputPathsAsCSV(paths, name, dictEC50)
    

wrapper("20_39 EC50s.csv")
wrapper("36_37 EC50s.csv")