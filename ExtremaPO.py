
import matplotlib.pyplot as plt

# For a given time, check to see if epsilon nbhd of next time steps image intersects
# epsilon nbhd of time's image. If so, add to time's component list and check next time, etc.
def ForwardCheck(t,ts,epsilon,compList):
	if t < len(ts)-1:
		next = t+1
		while (ts[next] + epsilon) >= (ts[t] - epsilon) and (ts[next] - epsilon) <= (ts[t] + epsilon):
			var = 1
			for time in compList[t]:
				if not ((ts[next] + epsilon) >= (ts[time] - epsilon) and (ts[next] - epsilon) <= (ts[time] + epsilon)):
					var = 0
					return
			if var:
				compList[t].append(next)
			if next == len(ts) - 1:
				return
			next += 1

# For a given time, check to see if epsilon nbhd of prev time steps image intersects
# epsilon nbhd of time's image. If so, add to time's component list and check prev time, etc.
def ReverseCheck(t,ts,epsilon,compList):
	if t > 0:
		prev = t-1
		while (ts[prev] + epsilon) >= (ts[t] - epsilon) and (ts[prev] - epsilon) <= (ts[t] + epsilon):
			var = 1
			for time in compList[t]:
				if not ((ts[prev] + epsilon) >= (ts[time] - epsilon) and (ts[prev] - epsilon) <= (ts[time] + epsilon)):
					var = 0
					return
			if var:
				compList[t].insert(0, prev)
			if prev == 0:
				return
			prev -= 1

# For a given epsilon, build component list for each t
def BuildCompList(ts,epsilon):
	compList = []
	for t in range(0,len(ts)):
		compList.append([t])
		ForwardCheck(t,ts,epsilon,compList)
		ReverseCheck(t,ts,epsilon,compList)
	return compList

	# Normalize time series and find min/max
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

# Build epsilon indexed list of components based on percent of max-min
#def BuildEIList(ts):
#	newts,minVal,maxVal = Normalize(ts)
#	eiList = []
#	step = 0.05*(abs(maxVal-minVal))
#	epsilon = 0
#	while epsilon <= 0.55*(maxVal-minVal):
#		eiList.append(BuildCompList(ts,epsilon))
#		epsilon += step
#	return eiList

# Build epsilon indexed list of components based on normalizing
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

# Uniqify a list
def Uniqify(inputList):
	for item1 in inputList:
		tempList = list(inputList)
		tempList.remove(item1)
		for item2 in tempList:
			if item1 == item2:
				inputList.remove(item1)
	return inputList

# Create a list of min and max components
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

# Count epsilon-life of mins/maxes and return list of min/max lifetimes
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

# Extract deepest lifetime mins/maxes
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

# Given a min and max, find the maximum epsilon where their components are disjoint
def FindEps(eiList, minTime, maxTime):
	for ndx in range(0,len(eiList)):
		if len(set(eiList[ndx][minTime]).intersection(eiList[ndx][maxTime])) == 0:
			step = ndx
	return step

# Process a list of time series and output a list of time series info. Each item in the list will correspond to time series
# and will be a list of the form [eiList, minTime, maxTime, step]
# eiList = epsilon-indexed list of components, min/maxTime = time of chosen min/max, step = highest step of epsilon at which
# the min component and max component are disjoint
# param specifies how the min/max are chosen: 0 = deepest, 1 = first
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

# Find the minimum of all steps from each ts
def FindMaxEps(sumList):
	eps = sumList[0][3]
	for ts in sumList:
		if ts[3] < eps:
			eps = ts[3]
	return eps

# For each ts, pull min/max components for each step up to eps indexed by ts then step
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

def ParseFile(fileName):
	f = open(fileName)
	TSData = []
	for ndx in range(0,6): #range(0,len(lineList)-1):
		TSData.append([])
	for line in f:
		lineList = line.split()
		for ndx in range(0,6): #range(0,len(lineList)-1):
			TSData[ndx].append(float(lineList[ndx+1]))
	f.close()
	return TSData

def PlotPercent(percentPOList,step):
	epsList = []
	for ndx in range(0,len(percentPOList)):
		epsList.append(ndx*step)
	plt.plot(epsList, percentPOList, 'ro')
	plt.axis([0, 0.6, 0, 1])
	plt.ylabel('percent of linear order')
	plt.xlabel('epsilon')
	plt.show()

def main():
	step = 0.01
	tsList = ParseFile('SampleData.tsv')
	sumList = ProcessTS(tsList,0,step)
	maxEps = FindMaxEps(sumList)
	minMaxCompList = PullMinMaxComps(sumList,maxEps,step)
	POsumList = BuildPO(minMaxCompList,maxEps,step)
	percentPOList = ConvertPOsumList(POsumList)
	PlotPercent(percentPOList,step)

main()
