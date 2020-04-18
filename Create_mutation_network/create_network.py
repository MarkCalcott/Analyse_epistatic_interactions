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

def wrapper(filename, start = '0000000', end = '1111111', increaseTolerance = 1):
    # Create dictionaries for network graph
    dictEC50 = create_EC50_dictionary(filename)
    dictNetwork = create_network_dictionary(dictEC50.keys())
    
    # Filename specific to data set
    name = filename.split()[0]
    
    #Save the edges of a graph for cytoscape
    saveEdges(name, dictEC50, dictNetwork)
    
    #Save the nodes of a graph for cytoscape
    saveNodes(name, dictEC50, dictNetwork)
    

wrapper("20_39 EC50s.csv")
wrapper("36_37 EC50s.csv")