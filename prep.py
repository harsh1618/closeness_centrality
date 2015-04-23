from networkx import *
from math import sin, cos, sqrt, atan2, radians
def dist(lat1,lon1,lat2,lon2):
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
G=Graph()
f=open("cal.cedge.txt")
elist=[]
for line in f:
	l = line.split(" ")[1:]
	elist.append(tuple((int(l[0]),int(l[1]),float(l[2]))))
G.add_weighted_edges_from(elist)

if not is_connected(G):
	exit()
print("Closeness centrality")
count=0
nlistk=[]
for i in G.nodes():
	c=closeness_centrality(G,i,distance=True)
	# print i,c
	nlistk.append((i,c))
	count+=1
	if count>100:
		break
k=10

nlistk.sort(key=lambda tup: tup[1],reverse=True)
for i in nlistk:
	print i[1]
exit()
# for i in G.nodes():
# 	c=closeness_centrality(G,i)
# 	count+=1
# 	if count%100==0:
# 		print count
# k=10
# x=[tuple(v,c[v]) for v in G.nodes()]
# topk=x.sort(key=lambda tup: tup[1],reverse=True)[:k]

nlist=[]
latDict={}
lonDict={}
f=open("cal.cnode.txt")
for line in f:
	l = line.split(" ")
	latDict[int(l[0])]=float(l[1])
	lonDict[int(l[0])]=float(l[2])
set_node_attributes(G, 'lat', latDict)
set_node_attributes(G, 'lon', lonDict)
count=0
for i in G.nodes():
	# print str(i)+" "+str(G.node[i]['lat'])+" "+str(G.node[i]['lon'])
	count+=1
	if count>100:
		break
