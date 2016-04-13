import sys
import matplotlib.pyplot as plt
from numpy import matrix,copy
from graphviz import Digraph
from PIL import Image
import intervalgraph as ig
import json
import heapq

# Checks to see if the epsilon nbhd of value of t2 intersects the epsilon nbhd of value of t2
def BoolIntersect(t1,t2,ts,epsilon):
	if (ts[t2] + epsilon) >= (ts[t1] - epsilon) and (ts[t2] - epsilon) <= (ts[t1] + epsilon):
		return True

# Grows epsilon components for a given time t. It does this by simultaneiously checking both neighboring points,
# if they exist, and only adding them if they intersect time t and all other times in the component, as well as
# intersecting eachother.
def GrowComponent(t,ts,epsilon,compList):
	index = 1
	if t == 0:
		while not(t+index > len(ts) - 1):
			for time in compList[t]:
				if not(BoolIntersect(t+index,time,ts,epsilon)):
					return
			compList[t].append(t+index)
			index += 1
	elif t == len(ts) - 1:
		while not(t-index < 0):
			for time in compList[t]:
				if not(BoolIntersect(t-index,time,ts,epsilon)):
					return
			compList[t].insert(0, t-index)
			index += 1
	elif not(BoolIntersect(t,t-index,ts,epsilon)) and BoolIntersect(t,t+index,ts,epsilon):
		while not(t+index > len(ts) - 1):
			for time in compList[t]:
				if not(BoolIntersect(t+index,time,ts,epsilon)):
					return
			compList[t].append(t+index)
			index += 1
	elif BoolIntersect(t,t-index,ts,epsilon) and not(BoolIntersect(t,t+index,ts,epsilon)):
		while not(t-index < 0):
			for time in compList[t]:
				if not(BoolIntersect(t-index,time,ts,epsilon)):
					return
			compList[t].insert(0, t-index)
			index += 1
	else:
		while BoolIntersect(t,t-index,ts,epsilon) and BoolIntersect(t,t+index,ts,epsilon):
			if BoolIntersect(t+index,t-index,ts,epsilon):
				for time in compList[t]:
					if not(BoolIntersect(t-index,time,ts,epsilon) and BoolIntersect(t+index,time,ts,epsilon)):
						return
				compList[t].append(t+index)
				compList[t].insert(0, t-index)
				if t+index == len(ts) - 1:
					index += 1
					while not(t-index < 0):
						for time in compList[t]:
							if not(BoolIntersect(t-index,time,ts,epsilon)):
								return
						compList[t].insert(0, t-index)
						index += 1
					return
				if t-index == 0:
					index += 1
					while not(t+index > len(ts) - 1):
						for time in compList[t]:
							if not(BoolIntersect(t+index,time,ts,epsilon)):
								return
						compList[t].append(t+index)
						index += 1
					return
				index += 1
				if not(BoolIntersect(t,t-index,ts,epsilon)) and BoolIntersect(t,t+index,ts,epsilon):
					while not(t+index > len(ts) - 1):
						for time in compList[t]:
							if not(BoolIntersect(t+index,time,ts,epsilon)):
								return
						compList[t].append(t+index)
						index += 1
					return
				if BoolIntersect(t,t-index,ts,epsilon) and not(BoolIntersect(t,t+index,ts,epsilon)):
					while not(t-index < 0):
						for time in compList[t]:
							if not(BoolIntersect(t-index,time,ts,epsilon)):
								return
						compList[t].insert(0, t-index)
						index += 1
					return
			else:
				return

# For a given epsilon, build component list for each t
# Inputs: ts = time series, epsilon
# Outputs: compList = list of components indexed by time
def BuildCompList(ts,epsilon):
	compList = []
	for t in range(0,len(ts)):
		compList.append([t])
		GrowComponent(t,ts,epsilon,compList)
	return compList

# Normalize time series and find global min/max
# Inputs: ts
# Outputs: newts = normalized time series, minVal = global min, maxVal = global max
def Normalize(ts):
	maxVal = ts[0]
	minVal = ts[0]
	for value in ts:
		if value > maxVal:
			maxVal = value
		if value < minVal:
			minVal = value
	newts = []
	for value in ts:
		newValue = float((value - minVal)) / (maxVal - minVal)
		newts.append(newValue)
	return newts,minVal,maxVal

# Build epsilon indexed list whose entries are lists of components generated at that epsilon
# Inputs: ts = time series, step = user defined parameter that specifies the increase in epsilon 
# Outputs: eiList = epsilon indexed list
def BuildEIList(ts,step):
	newts,minVal,maxVal = Normalize(ts)
	eiList = []
	epsilon = 0
	while epsilon <= 0.55:
		eiList.append(BuildCompList(newts,epsilon))
		epsilon += step
	return eiList

# Concatenate a list of lists
def Concatenate(list):
	newList = []
	for item in list:
		newList += item
	return newList

# Uniqify a list, i.e. turn it into a set
def Uniqify(inputList):
	for item1 in inputList:
		tempList = list(inputList)
		tempList.remove(item1)
		for item2 in tempList:
			if item1 == item2:
				inputList.remove(item1)
	return inputList

# Sorts the unique components grown into two lists, minList and maxList
# Inputs: eiList, ts
# Outputs: minList = list containing all components that are minima, maxList = similar
def MinMaxLabel(eiList,ts):
	compList = Uniqify(Concatenate(eiList))
	minList = []
	maxList = []
	for comp in compList:
		if 0 in comp:
			if not(len(ts)-1 in comp):
				if ts[comp[-1]+1] > ts[comp[-1]]:
					minList.append(comp)
				else:
					maxList.append(comp)
		elif len(ts)-1 in comp:
			if not(0 in comp):
				if ts[comp[0]-1] > ts[comp[0]]:
					minList.append(comp)
				else:
					maxList.append(comp)
		else:
			if (ts[comp[0]-1] > ts[comp[0]]) and (ts[comp[-1]+1] > ts[comp[-1]]):
				minList.append(comp)
			if (ts[comp[0]-1] < ts[comp[0]]) and (ts[comp[-1]+1] < ts[comp[-1]]):
				maxList.append(comp)
	return minList,maxList

# Return component-chain list indexed by time and a labeled(min/max/n/a) chain list
# Creates list indexed by time where each entry is an epsilon indexed list of components grown at that
# time. This is called chainList. Then I assign each component a label, 'min'/'max'/'n/a' and record
# this as labeledChains.
def BuildChains(eiList,ts,step):
	minList,maxList = MinMaxLabel(BuildEIList(ts,step),ts)
	chainList = []
	labeledChains = []
	for i in range(0,len(eiList[0])):
		chainList.append([])
		labeledChains.append([])
	for list in eiList:
		ndx = 0
		for item in list:
			chainList[ndx].append(item)
			if item in minList:
				labeledChains[ndx].append('min')
			elif item in maxList:
				labeledChains[ndx].append('max')
			else:
				labeledChains[ndx].append('n/a')
			ndx+=1
	return chainList,labeledChains

# Count number of epsilon steps that mins/maxes persisted and return list of min/max lifetimes.
# Note, only considers those that are local mins/maxes at epsilon = 0. 
# Outputs: minLife = time indexed list initialized to all 0's, time's entry is # of consecutive 
# epsilon steps that time was a min, maxLife = similar
def EpsLife(labeledChains):
	minLife = [0] * len(labeledChains)
	maxLife = [0] * len(labeledChains)
	for t in range(0,len(labeledChains)):
		if labeledChains[t][0] == 'min':
			s = 0
			while labeledChains[t][s] == labeledChains[t][0]:
				minLife[t] += 1
				s+=1
		elif labeledChains[t][0] == 'max':
			s = 0
			while labeledChains[t][s] == labeledChains[t][0]:
				maxLife[t] += 1
				s+=1
	return minLife,maxLife

# Extract the n deepest lifetime mins and maxes. If there are ties, the sequentially first one is chosen
# return list of 2n events, ordered by highest min, highest max, second highest min, ...
def DeepLife(minLife,maxLife,ts,n):
	minLifeCopy = list(minLife)
	maxLifeCopy = list(maxLife)
	deepEventList = []
	for ndx in range(0,n):
		minimum = max(minLifeCopy)
		maximum = max(maxLifeCopy)
		minIndex = minLifeCopy.index(minimum)
		maxIndex = maxLifeCopy.index(maximum)
		deepEventList.append(minIndex)
		deepEventList.append(maxIndex)
		minLifeCopy[minIndex] = 0
		maxLifeCopy[maxIndex] = 0
	return deepEventList

# Given a min and max, find the maximum epsilon step where their components are disjoint
def FindEps(eiList, deepEventList):
	value = 0
	epsilon = 0
	maxEps = -1
	while value == 0:
		for ndx1 in range(0,len(deepEventList)):
			for ndx2 in range(0,len(deepEventList)):
				if len(set(eiList[epsilon][deepEventList[ndx1]]).intersection(eiList[epsilon][deepEventList[ndx2]])) != 0 and ndx1 != ndx2:
					value = 1
		if value == 0:
			maxEps = epsilon
		epsilon += 1
	return maxEps

# Process a list of time series and output a list of time series info. Each item in the list will correspond to time series
# and will be a list of the form [eiList, minTime, maxTime, eps]
# eiList = epsilon-indexed list of components, deepMin/MaxList = list of times of first n mins/maxes
# eps = highest step of epsilon at which the min component and max component are disjoint
def ProcessTS(tsList, n, step):
	sumList = []
	for ts in tsList:
		eiList = BuildEIList(ts,step)
		chainList, labeledChains = BuildChains(eiList,ts,step)
		minLife,maxLife = EpsLife(labeledChains)
		deepEventList = DeepLife(minLife,maxLife,ts,n)
		eps = FindEps(eiList,deepEventList)
		sumList.append([eiList,deepEventList,eps])
	return sumList

# Find the minimum value of eps such that min/max interval of every time series being processed
# is disjoint. So take min of all eps
def FindMaxEps(sumList):
	epsilon = sumList[0][2]
	for ts in sumList:
		if ts[2] < epsilon:
			epsilon = ts[2]
	return epsilon

# For each ts, pull min/max components for each step up to eps indexed by ts
# then highest min, highest max, second highest min, ...
def PullEventComps(sumList,maxEps,step,n):
	eventCompList = []
	ndx = 0
	for tsList in sumList:
		eventCompList.append([])
		for event in range(0,2*n):
			eventCompList[ndx].append(tsList[0][maxEps][tsList[1][event]])
		ndx += 1
	return eventCompList

# Build partial order list indexed by epsilon. Each PO is a list indexed by 1st ts highest min, 1st ts highest max, 
# 1st ts second highest min, 1st ts second highst max, ..., last ts nth highest min, last ts nth highest max
# Each entry is a list of all mins/maxes that occur before the given min/max. 
def BuildPO(eventCompList,maxEps,step,n):
	PO = []
	for ts in range(0,len(eventCompList)):
		for event in range(0,2*n):
			PO.append([])
			for ndx in range(0,len(eventCompList)):
				for ndx1 in range(0,2*n):
					fixed = eventCompList[ts][event]
					checker = eventCompList[ndx][ndx1]
					intSize = len(set(fixed).intersection(checker))
					if intSize == 0:
						if fixed[0] < checker[0]:
							PO[2*n*ts + event].append(2*n*ndx + ndx1)
	return PO

# Convert PO's to graph class
def POToGraph(PO,TSLabels,n):
	G = ig.Graph()
	for value in range(0,len(PO)):
		G.add_vertex(value,TSLabels[value/(2*n)])
	for i in range(0,len(PO)):
		for j in PO[i]:
			G.add_edge(i,j)
	G = G.transitive_reduction()
	return G

# Convert graphs to digraphs
def GraphToDigraph(G):
	DG = Digraph(comment = 'PO')
	for v in G.vertices():
		DG.node(str(v),G.vertex_label(v))
	for e in G.edges():
		DG.edge(str(e[0]),str(e[1]))
	DG.render('graph.gv',view=True)

# Labeling function for the genes
def CreateLabel(sumList):
	d = len(sumList)
	label = [0]*(2*d)
	for ndx in range(0,d):
		lastEvent = max(sumList[ndx][1])
		indexOfLE = sumList[ndx][1].index(lastEvent)
		if indexOfLE%2 == 1:   # last event was min, since ordered min,max,min,max,...
			label[2*d - 1 - ndx] = 1
			label[d - 1 - ndx] = 0
		else:
			label[2*d - 1 - ndx] = 0
			label[d - 1 - ndx] = 1
	label = '0b' + ''.join([str(x) for x in label])	#Cast to binary
	label = int(label,2)
	return label

# Create JSON string for each graph
def ConvertToJSON(graph,sumList,TSLabels):
	G = graph
	output = {}
	output["poset"] = [ list(G.adjacencies(i)) for i in G.vertices() ]
	output["events"] = [ TSLabels.index(G.vertex_label(i)) for i in G.vertices() ]
	output["label"] = CreateLabel(sumList)
	output["dimension"] = len(TSLabels)
	with open('pattern.json', 'w') as fp:
	  json.dump(output, fp)

# Parse file where genes are in columns output list of TS and list of TS labels 
def ParseColFile(fileName):
	f = open(fileName)
	TSData = []
	TSLabels = []
	timeStepList = []

	value = 0 												# indicates if line of gene labels has been reached
	for line in f:
		if line.split()[0] != '#':
			if value == 0:
				value = 1
				for ndx in range(0,len(line.split()) - 1):
					TSData.append([])
					TSLabels.append(line.split()[ndx+1])
			else:
				for ndx in range(0,len(line.split()) - 1):
					TSData[ndx].append(float(line.split()[ndx+1]))
					timeStepList.append(float(line.split()[0]))
	f.close()
	return TSData,TSLabels,timeStepList

# Parse file with genes in row format
def ParseRowFile(fileName):
	f = open(fileName)
	TSList = []
	TSLabels = []

	value = 0 											# indicates if time_points line has been reached
	for line in f:
		if line.split()[0] != '#':
			if value == 0:
				timeStepList = [float(time) for time in line.split()[1:]]
				value = 1
			else:
				TSList.append([float(item) for item in line.split()[1:]])
				TSLabels.append(line.split()[0])
	f.close()
	return TSList,TSLabels, timeStepList

# Parse network file to know which time series to pick off
def ParseNetworkFile(fileName):
	f = open(fileName)
	chosenTSList = []
	for line in f:
		chosenTSList.append(line.split()[0])
	f.close()
	return chosenTSList

# Pull the data corresponding to those TS returned by ParseNetworkFile
def PickNetworkTS(TSList,TSLabels,chosenTS):
	newTSList = []
	for ts in chosenTS:
		tsIndex = TSLabels.index(ts)
		newTSList.append(TSList[tsIndex])
	return newTSList

# Truncate timeseries. Keep all data before timeCutOff
def TruncateTS(newTSList,timeStepList,timeCutOff):
	indexOfCutOff = timeStepList.index(timeCutOff)
	truncatedTSList = []
	for ts in newTSList:
		truncatedTSList.append(ts[:(indexOfCutOff+1)])
	print(truncatedTSList)
	return truncatedTSList

# The arguments are dataFileName, fileType = 'row' or 'col', networkFileName,
# n = number of mins/maxes to pull, timeCutOff = ignore data after this time ( = -1 if no cutOff), step (default to 0.01)
def main():
	dataFileName = sys.argv[1]
	fileType = sys.argv[2]
	networkFileName = sys.argv[3]
	n = int(sys.argv[4])
	timeCutOff = float(sys.argv[5])
	if len(sys.argv) == 6:
		step = 0.01
	else:
		step = float(sys.argv[6])

	if fileType == 'col':
		TSList = ParseColFile(dataFileName)[0]
		TSLabels = ParseColFile(dataFileName)[1]
		timeStepList = ParseColFile(dataFileName)[2]
	elif fileType == 'row':
		TSList = ParseRowFile(dataFileName)[0]
		TSLabels = ParseRowFile(dataFileName)[1]
		timeStepList = ParseRowFile(dataFileName)[2]

	newTSLabels = ParseNetworkFile(networkFileName)
	newTSList = PickNetworkTS(TSList,TSLabels,newTSLabels)

	if timeCutOff != float(-1):
		newTSList = TruncateTS(newTSList,timeStepList,timeCutOff)

	sumList = ProcessTS(newTSList,n,step)
	maxEps = FindMaxEps(sumList)
	eventCompList = PullEventComps(sumList,maxEps,step,n)
	PO = BuildPO(eventCompList,maxEps,step,n)
	graph = POToGraph(PO,newTSLabels,n)
	ConvertToJSON(graph,sumList,newTSLabels)

	# # Prints the PO's from the conversion to S.H.'s graph class
	# GraphToDigraph(graph)
	
main()
