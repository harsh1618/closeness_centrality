import networkx as nx
from math import sin, cos, sqrt, atan2, radians

def geographicalDistance(lat1,lon1,lat2,lon2):
    '''
    Compute the distance along the earth's surface between
    the two given coordinates.
    '''
    R = 6373.0
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c
    return distance

def getGraph():
    '''
    Read the California road network graph from the data files.
    '''
    G = nx.Graph()
    f=open("data/cal.cedge.txt")
    elist=[]
    for line in f:
        l = line.split(" ")[1:]
        elist.append(tuple((int(l[0]),int(l[1]),float(l[2]))))
    G.add_weighted_edges_from(elist)

    nlist=[]
    latDict={}
    lonDict={}
    f=open("data/cal.cnode.txt")
    for line in f:
        l = line.split(" ")
        latDict[int(l[0])]=float(l[1])
        lonDict[int(l[0])]=float(l[2])
    nx.set_node_attributes(G, 'lat', latDict)
    nx.set_node_attributes(G, 'lon', lonDict)
    return G
