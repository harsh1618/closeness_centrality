#!/usr/bin/env python
'''
Implementation of efficient top-k closeness centrality search
algorithm given by Olsen, Labouseur and Hwang, 2014.

Authors:
Harshvardhan Sharma
Malay Kumar Singh
Dhruv Anand
'''

import networkx as nx
import heapq
from copy import copy
from random import randint
import sys
import matplotlib.pyplot as plt
from time import time

sys.setrecursionlimit(100000)

k = 100

if len(sys.argv) < 2:
    print "Usage:", sys.argv[0], "<type>"
    print "type should be 0 for brute force and 1 for delta-PFS method"
    exit(1)

if sys.argv[1] == '0':
    print "Running brute force algo"
    brute_force = True
else:
    brute_force = False

SHOW_GRAPHS = False
DIRECTED = False
RANDOM_GRAPH = True

if RANDOM_GRAPH:
    random_graph = nx.fast_gnp_random_graph(30000, 0.00008)
    # select largest connected component
    connected_components = list(nx.connected_component_subgraphs(random_graph))
    connected_components.sort(key = len, reverse = True)
    G = connected_components[0]
    # add weights
    for source, dest in G.edges():
        G[source][dest]['weight'] = randint(1, 100)
else:
    G = nx.Graph()
    G.add_nodes_from([0,1,2,3,4])
    G.add_weighted_edges_from([(0,1,1), (0,2,1), (0, 3, 1), (0,4,1), (0,5,1)])
V = G.number_of_nodes()
print "Number of vertices", V

if SHOW_GRAPHS:
    nx.draw(G)
    plt.show()

L = {} #vertex-level map

A = [] #priority queue

phi = {} # min of upper bound on centrality - threshold for each vertex
pruned = {}
upper_bounds = {}
lower_bounds = {}
schedule = nx.DiGraph()

# wrapper functions to implement max heap of (centrality, vertex)
def heappush(i):
    heapq.heappush(A, (-i[0], i[1]))

def heappop():
    top = heapq.heappop(A)
    return (-top[0], top[1])

def nlargest(k):
    return [(-x[0], x[1]) for x in heapq.nsmallest(k, A)]


def process(vertex, s, delta):
    '''
    Recursively process the graph and update centrality values.

    'vertex' is the vertex to be processed.
    s is the sum of distances to all vertices reachable from 'vertex'.
    'delta' is the distance upper bound for the previous vertex.
    '''
    centrality = float((len(L) - 1)**2) / (s * (V - 1))
    #centrality = float((len(filter(lambda x: L[x] is not None, L)) - 1)**2) / (s * (V - 1))
    heappush((centrality, vertex))

    # prune if upper bound of centrality is less than k'th largest
    #if len(A) >= k:
    #    prune_threshold = nlargest(k)[-1][0]
    #    prune(vertex, s, prune_threshold, delta)
    #if vertex in pruned: return

    for v in schedule[vertex]:
        if v in pruned: continue
        s_new, delta_new, old_L = deltaPFS(v, vertex, s, delta)
        process(v, s_new, delta_new)
        #rollback
        for key in old_L:
            L[key] = old_L[key]

def PFS(v):
    '''
    Perform Priority First Search from vertex v.
    PFS is a modified version of Dijkstra's algorithm that calculates
    additional values required for centrality computation and pruning.

    Returns
    s: sum of distances to all reachable vertices
    delta: distance upper bound
    '''
    L[v] = 0
    s = 0
    delta = 0
    Q = [(0, v)]
    popped = {}
    while len(Q) > 0:
        l, n = heapq.heappop(Q)
        # skip if node already encountered (invalid due to a previous decrement key)
        if n in popped: continue
        else: popped[n] = True
        s += l
        delta = max(delta, l)
        for v_prime in G[n]:
            l_prime = l + G[n][v_prime]['weight']
            if L.get(v_prime) is None:
                L[v_prime] = l_prime
                heapq.heappush(Q, (l_prime, v_prime))
            elif l_prime < L[v_prime]:
                L[v_prime] = l_prime
                # no decrement key operation in heapq, push the new value
                heapq.heappush(Q, (l_prime, v_prime))
    return s, delta

def deltaPFS(v, previous, s, delta):
    '''
    Perform incremental PFS from vertex v using distance computations
    done during the PFS or deltaPFS of 'previous'.

    'previous' is the vertex processed immediately before v.
    's' and 'delta' are the sum of distances to reachable nodes
    and the distance upper bound for the previous vertex respectively.

    Returns
    s: sum of distances to all reachable vertices
    delta: distance upper bound
    old_L: preserved vertex levels to be used for rollback.
    '''
    weight = G[v][previous]['weight']
    alpha_p = L[previous]
    alpha_v = alpha_p - weight
    old_L = {}
    old_L[v] = L[v]
    Q = [(alpha_v, v)]
    if not DIRECTED: s-= 2
    s += weight * len(L)
    #s += weight * len(filter(lambda x: L[x] is not None, L))
    L[v] = alpha_v
    delta += weight
    popped = {}
    while len(Q) > 0:
        l, n = heapq.heappop(Q)
        if n in popped: continue
        else: popped[n] = True

        if old_L.get(n) is not None:
            s += (l - alpha_v)
            delta = max(delta, l-alpha_v)
        else:
            s -= (old_L[n] - l)

        for v_prime in G[n]:
            l_prime = l + G[n][v_prime]['weight']
            l_prime2 = L.get(v_prime)
            if l_prime2 is None or l_prime < l_prime2:
                L[v_prime] = l_prime
                heapq.heappush(Q, (l_prime, v_prime))
                if v_prime not in old_L:
                    old_L[v_prime] = l_prime2
    return s, delta, old_L

def prune(vertex, s, threshold, delta):
    '''
    Check if deltaPFS computation for the vertex can be skipped.
    s is the sum of distances to all reachable nodes.
    threshold is the k'th largest value in the closeness centrality priority queue.
    delta is the distance upper bound.
    
    Adds vertices that can be pruned to the global 'pruned' dict
    '''
    Q = [(0, vertex)]
    popped = {}
    while len(Q) > 0:
        l, f = heapq.heappop(Q)
        if f in popped: continue
        else: popped[f] = True

        # simplifications for undirected graphs
        if not DIRECTED:
            upper_bound = V
            lower_bound = V
        else:
            upper_bound = upper_bounds.get(f, V)
            lower_bound = lower_bounds.get(f, 0)

        s_prime = s - (V - lower_bound)*delta - l * upper_bound
        c_prime = float(upper_bound - 1)**2 / (V-1)*(s_prime)
        if f not in phi:
            phi[f] = c_prime - threshold
        elif c_prime - threshold < phi[f]:
            phi[f] = c_prime - threshold
            for f_prime in G[f]:
                heapq.heappush(Q, (l + G[f][f_prime]['weight'], f_prime))
            if len(schedule[f]) == 0 and c_prime < threshold:
                pruned[f] = True

def getSchedule(G):
    '''
    Generate a schedule for processing the nodes so as to minimize
    total running time.

    Returns:
    start_vertices: a list of vertices with in degree 0 in the schedule
    schedule: a directed networkx graph in which presence of an edge
              (u, v) implies that v is processed after u.
    '''
    mst = nx.minimum_spanning_tree(G)
    schedule = nx.DiGraph()
    for a,b in mst.edges():
        schedule.add_edge(a, b)
    if SHOW_GRAPHS:
        nx.draw(schedule)
        plt.show()
    start_vertices = []
    for v, in_degree in schedule.in_degree_iter():
        if in_degree == 0:
            start_vertices.append(v)
    return start_vertices, schedule

if brute_force:
    start = time()
    i = 0
    for vertex in G.nodes():
        i += 1
        if i%1000 == 0: print i
        L = {}
        s, delta = PFS(vertex)
        centrality = float((len(L) - 1)**2) / (s * (V - 1))
        heappush((centrality, vertex))
    print "Time taken:", int(time() - start)

else:
    start_vertices, schedule = getSchedule(G) # schedule is a networkx graph, start_vertices is a list

    start = time()
    for start_vertex in start_vertices:
        L = {}
        phi = {}
        pruned = {}
        upper_bounds = {}
        lower_bounds = {}
        s, delta = PFS(start_vertex)
        process(start_vertex, s, delta)
    print "Time taken:", int(time() - start)

for i in nlargest(k):
    print i[1], i[0]

