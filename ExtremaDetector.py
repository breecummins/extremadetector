# Extrema Detector First run

# Parser
def ParseFile(fname):
	f=open(fname,'r')
	for _ in range(5):
		f.readline()
	genelist=f.readline().split()[1:]
	times=[[float(w) for w in l.split()[1:]] for l in f.readlines()]
	timeseries=numpy.array(times).transpose()
	f.close()
	return timeseries

# Compute Candidate list
def CandidateList(ts):
	candList = [0]
	for t in range(1,len(ts) - 1):
		if (ts[t+1] >= ts[t] and ts[t-1] >= ts[t]) or (ts[t+1] <= ts[t] and ts[t-1] <= ts[t]):
			candList.append(t)
	candList.append(len(ts) - 1)
	return candList

# For a given time, check to see if epsilon nbhd of next time steps image intersects
# epsilon nbhd of time's image. If so, add to time's component list and check next time, etc.
def ForwardCheck(t,ts,epsilon,compDict):
	if t <= len(ts)-2:
		next = t+1
		while (ts[next] + epsilon) >= (ts[t] - epsilon) and (ts[next] - epsilon) <= (ts[t] + epsilon):
			compDict[t].append(next)
			next += 1
			print(next,compDict)
			if next == len(ts):
				return

# For a given time, check to see if epsilon nbhd of prev time steps image intersects
# epsilon nbhd of time's image. If so, add to time's component list and check prev time, etc.
def ReverseCheck(t,ts,epsilon,compDict):
	if t > 0:
		prev = t-1
		while (ts[prev] + epsilon) >= (ts[t] - epsilon) and (ts[prev] - epsilon) <= (ts[t] + epsilon):
			compDict[t].insert(0,prev)
			prev -= 1
			print(compDict)
			if prev == 0:
				return

# For a given epsilon, use candidate list to build component dictionary
def BuildCompDict(ts,epsilon,candList):
	compDict = {}
	for t in candList:
		compDict[t] = [t]
		ForwardCheck(t,ts,epsilon,compDict)
		ReverseCheck(t,ts,epsilon,compDict)
	return compDict

# Main
def main():
	#data = ParseFile('TestData.tsv')
	#print(data)
	ts = [2,1,2,4,5,4,5]
	candList = CandidateList(ts)
	print(BuildCompDict(ts,1,candList))

main()
