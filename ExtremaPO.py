
import sys
import matplotlib.pyplot as plt
from numpy import matrix,copy
from graphviz import Digraph
from PIL import Image

# Checks to see if the epsilon nbhd of value of t2 intersects the epsilon nbhd of value of t2
def BoolIntersect(t1,t2,ts,epsilon):
	if (ts[t2] + epsilon) >= (ts[t1] - epsilon) and (ts[t2] - epsilon) <= (ts[t1] + epsilon):
		return True

# Grows epsilon components for a given time t
def GrowComponent(t,ts,epsilon,compList):
	index = 1
	if t == 0:
		while not(t+index > len(ts) - 1) and BoolIntersect(t,t+index,ts,epsilon):
			for time in compList[t]:
				if not(BoolIntersect(t+index,time,ts,epsilon)):
					return
			compList[t].append(t+index)
			index += 1
	elif t == len(ts) - 1:
		while not(t-index < 0) and BoolIntersect(t,t-index,ts,epsilon):
			for time in compList[t]:
				if not(BoolIntersect(t-index,time,ts,epsilon)):
					return
			compList[t].insert(0, t-index)
			index += 1
	elif not(BoolIntersect(t,t-index,ts,epsilon)) and BoolIntersect(t,t+index,ts,epsilon):
		while not(t+index > len(ts) - 1) and BoolIntersect(t,t+index,ts,epsilon):
			for time in compList[t]:
				if not(BoolIntersect(t+index,time,ts,epsilon)):
					return
			compList[t].append(t+index)
			index += 1
	elif BoolIntersect(t,t-index,ts,epsilon) and not(BoolIntersect(t,t+index,ts,epsilon)):
		while not(t-index < 0) and BoolIntersect(t,t-index,ts,epsilon):
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
					while not(t-index < 0) and BoolIntersect(t,t-index,ts,epsilon):
						for time in compList[t]:
							if not(BoolIntersect(t-index,time,ts,epsilon)):
								return
						compList[t].insert(0, t-index)
						index += 1
					return
				if t-index == 0:
					while not(t+index > len(ts) - 1) and BoolIntersect(t,t+index,ts,epsilon):
						for time in compList[t]:
							if not(BoolIntersect(t+index,time,ts,epsilon)):
								return
						compList[t].append(t+index)
						index += 1
					return
				index += 1
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

# Extract deepest lifetime min and max
def DeepLife(minLife,maxLife):
	value = 0
	for t in range(0,len(minLife)):
		if minLife[t] > value:
			deepMin = t
			value = minLife[t]
	value = 0
	for t in range(0,len(maxLife)):
		if maxLife[t] > value:
			deepMax = t
			value = maxLife[t]
	return deepMin,deepMax

# Given a min and max, find the maximum epsilon step where their components are disjoint
def FindEps(eiList, minTime, maxTime):
	for ndx in range(0,len(eiList)):
		if len(set(eiList[ndx][minTime]).intersection(eiList[ndx][maxTime])) == 0:
			eps = ndx
	return eps

# Process a list of time series and output a list of time series info. Each item in the list will correspond to time series
# and will be a list of the form [eiList, minTime, maxTime, eps]
# eiList = epsilon-indexed list of components, min/maxTime = time of chosen min/max, eps = highest step of epsilon at which
# the min component and max component are disjoint
# param specifies how the min/max are chosen: 0 = deepest, 1 = first (to be implemented)
def ProcessTS(tsList, param, step):
	sumList = []
	for ts in tsList:
		eiList = BuildEIList(ts,step)
		chainList, labeledChains = BuildChains(eiList,ts,step)
		minLife,maxLife = EpsLife(labeledChains)
		if param == 0:
			minTime,maxTime = DeepLife(minLife,maxLife)
#		elif param == 1:
#			minTime,maxTime = FirstLife(minLife,maxLife)
		eps = FindEps(eiList,minTime,maxTime)
		sumList.append([eiList,minTime,maxTime,eps])
	return sumList

# Find the minimum value of eps such that min/max interval of every time series being processed
# is disjoint. So take min of all eps
def FindMaxEps(sumList):
	epsilon = sumList[0][3]
	for ts in sumList:
		if ts[3] < epsilon:
			eps = ts[3]
	return eps

# For each ts, pull min/max components for each step up to eps indexed by ts then step
# min component comes before max component
def PullMinMaxComps(sumList,maxEps,step):
	minMaxCompList = []
	ndx = 0
	for tsList in sumList:
		minMaxCompList.append([])
		for eps in range(0,maxEps+1):
			minMaxCompList[ndx].append([tsList[0][eps][tsList[1]],tsList[0][eps][tsList[2]]])
		ndx += 1
	return minMaxCompList

# Build partial order list indexed by epsilon. Each PO is a list indexed by 1st ts min, 1st ts max, ..., last ts min, last ts max
# Each entry is a list of all mins/maxes that occur before the given min/max. 
def BuildPO(minMaxCompList,maxEps,step):
	maxStep = maxEps
	POsumList = []
	for eps in range(0,maxStep+1):
		PO = []
		for ts in range(0,len(minMaxCompList)):
			for exType in range(0,2):
				PO.append([])
				for ndx in range(0,len(minMaxCompList)):
					for ndx1 in range(0,2):
						fixed = minMaxCompList[ts][eps][exType]
						checker = minMaxCompList[ndx][eps][ndx1]
						intSize = len(set(fixed).intersection(checker))
						if intSize == 0: #intSize < 2
#							if (intSize == 1 and len(fixed) > 1 and len(checker) > 1) or intSize == 0:
								if fixed[0] < checker[0]:
									PO[2*ts + exType].append(2*ndx + ndx1)
		POsumList.append(PO)
	return POsumList

# For each epsilon, for each PO in POsumList, sum number of inequalities in PO and divide by number in total order
# So find how close it is to a linear order
def ConvertPOsumList(POsumList):
	n = len(POsumList[0])
	countTO = (n*(n-1))/2.0
	percentPOList = []
	for PO in POsumList:
		count = 0
		for ndx in range(0,n):
			count += len(PO[ndx])
		percentPOList.append(round(count/countTO,4))
	return percentPOList

# Plot the percent of total order vs epsilon
def PlotPercent(percentPOList,step):
	epsList = []
	for ndx in range(0,len(percentPOList)):
		epsList.append(ndx*step)
	plt.plot(epsList, percentPOList, 'ro')
	plt.axis([0, 0.6, 0, 1])
	plt.ylabel('percent of linear order')
	plt.xlabel('epsilon')
	plt.savefig('percentPlot.png')
	f = Image.open("percentPlot.png").show()

# Pick off a few time series if do not want to look at all of them
def PickTS(TSList,chosenTS,TSLabels):
	newTSList = []
	newTSLabels = []
	ndx = 0
	for ts in TSList:
		if ndx in chosenTS:
			newTSList.append(TSList[ndx])
			newTSLabels.append(TSLabels[ndx])
		ndx += 1
	return newTSList,newTSLabels

# Convert each PO in POsumList into adjacency matrix
def ConvertPO(POsumList):
	matrixPOsumList = []
	for PO in POsumList:
		POmatrixList = []
		for ndx in range(0,len(PO)):
			POmatrixList.append([])
			for ndx1 in range(0,len(PO)):
				if ndx1 in PO[ndx]:
					POmatrixList[ndx].append(1)
				else:
					POmatrixList[ndx].append(0)
		matrixPOsumList.append(matrix(POmatrixList))
	return matrixPOsumList

# Reduce each PO matrix as much as posible. i.e. transitive reduction
# For graphing purposes
def ReducePOmatrix(matrixPOsumList):
	reducedSumList = []
	for PO in matrixPOsumList:
		mPO = matrix(copy(PO))
		colSum = sum(mPO)
		for ndx1 in range(0,mPO.shape[0]):
			if colSum[0,ndx1] > 1:
				for ndx2 in range(0,mPO.shape[0]):
					if mPO[ndx2,ndx1] == 1:
						for ndx3 in range(0,mPO.shape[0]):
							if mPO[ndx1,ndx3] == 1:
								mPO[ndx2,ndx3] = 0
		reducedSumList.append(mPO)
	return reducedSumList

# Graph Partial Orders
def GraphPO(reducedSumList,TSLabels):
	graphNum = 0
	for PO in reducedSumList:
		graph = Digraph(comment = graphNum)
		for value in range(0,PO.shape[0]):
			if value%2 == 0:
				label = ' min'
			else:
				label = ' max'
			graph.node(str(value),TSLabels[value/2] + label)
		for row in range(0,PO.shape[0]):
			for col in range(0,PO.shape[0]):
				if PO[row,col] == 1:
					graph.edge(str(col),str(row))
		graph.render('graph'+str(graphNum)+'.gv',view=True)
		graphNum += 1

# Print partial orders
def PrintPO(reducedSumList,newTSLabels):
	start = raw_input('To plot range of partial orders, enter starting epsilon step, else enter "x":\n')
	if not(start=='x'):
		end = input('Enter ending epsilon step:\n')
		newSumList = []
		for ndx in range(int(start),int(end+1)):
			newSumList.append(reducedSumList[ndx])
		GraphPO(newSumList,newTSLabels)

# Parse RawData.tsv file output list of TS and list of TS labels
def ParseFile(fileName):
	f = open(fileName)
	TSData = []
	for ndx in range(0,43):
		TSData.append([])
	lineNDX = 0
	TSLabels = []
	for line in f:
		lineList = line.split()
		if lineNDX == 5:
			for ndx in range(0,len(lineList)-1):
				TSLabels.append(lineList[ndx+1])
		elif lineNDX >= 6:
			for ndx in range(0,len(lineList)-1):
				TSData[ndx].append(float(lineList[ndx+1]))
		lineNDX += 1
	f.close()
	return TSData,TSLabels

def main():
	step = float(sys.argv[1])
	if len(sys.argv) == 2:
		param = 0
	else:
		param == sys.argv[2]

	# To do: Put this into sys.argv inputs. need to parse string to get it into list
	chosenTS = [2,3,11,14,16,17]
	
	TSList = ParseFile('RawData.tsv')[0]
	TSLabels = ParseFile('RawData.tsv')[1]
	newTSList,newTSLabels = PickTS(TSList,chosenTS,TSLabels)
	sumList = ProcessTS(newTSList,param,step)
	maxEps = FindMaxEps(sumList)
	minMaxCompList = PullMinMaxComps(sumList,maxEps,step)
	POsumList = BuildPO(minMaxCompList,maxEps,step)
	percentPOList = ConvertPOsumList(POsumList)
	matrixPOsumList = ConvertPO(POsumList)
	reducedSumList = ReducePOmatrix(matrixPOsumList)

	PlotPercent(percentPOList,step)
	PrintPO(reducedSumList,newTSLabels)
	

main()
