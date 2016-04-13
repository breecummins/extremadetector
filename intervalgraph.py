#!/usr/bin/python

import subprocess, copy, json
from collections import defaultdict

################################
# Directed Acyclic Graph Class #
################################

class Graph:
  def __init__(self):
    """ Initialize an empty graph object """
    self.vertices_ = set()
    self.adjacency_lists_ = {}
    self.vertex_labels_ = {}
    self.edge_labels_ = {}
  def add_vertex(self, v, label = ''):
    """ Add the vertex v to the graph and associate a label if one is given """
    if v in self.vertices_: return
    self.vertices_.add(v)
    self.adjacency_lists_[v] = set ()
    self.vertex_labels_[v] = label
    self.edge_labels_[v] = {}
  def add_edge(self, u, v, label = ''):
    """ Add the edge u -> v to the graph and associate a label if one is given """
    self.add_vertex(u)
    self.add_vertex(v)
    self.adjacency_lists_[u].add(v)
    self.edge_labels_[u][v] = label
  def remove_edge(self, u, v):
    """ Remove the edge u -> v from the graph """
    self.adjacency_lists_[u].discard(v)
    self.edge_labels_[u].pop(v, None)
  def vertex_label(self, v):
    """ Return the label on the vertex v """
    return self.vertex_labels_[v]
  def edge_label(self, u, v):
    """ Return the label on the edge u -> v """
    return self.edge_labels_[u][v]
  def vertices(self):
    """ Return the set of vertices in the graph """
    return self.vertices_
  def edges(self):
    """ Return a complete list of directed edges (u,v) in the graph """
    return [(u,v) for u in self.vertices() for v in self.adjacencies(u)]
  def adjacencies(self, v):
    """ Return the set of adjacencies of v, i.e. { u : v -> u } """
    return self.adjacency_lists_[v]
  def clone(self):
    """ Return a copy of this graph """
    return copy.deepcopy(self)
  def transpose(self):
    """ Return a new graph with edge direction reversed. """
    G = Graph ()
    for v in self.vertices(): G.add_vertex(v,self.vertex_label(v))
    for (u,v) in self.edges(): G.add_edge(v,u,self.edge_label(u,v))
    return G
  def transitive_closure(self):
    """ Return a new graph which is the transitive closure """
    G = self.clone ()
    for w in self.vertices():
      for u in self.vertices():
        for v in self.vertices():
          if w in G.adjacencies(u) and v in G.adjacencies(w):
            G . add_edge(u,v)
    return G
  def transitive_reduction(self):
    """ Return a new graph which is the transitive reduction """
    TC = self.transitive_closure ()
    G = self.clone ()
    for (u,v) in TC.edges():
      for w in TC.adjacencies(v):
        G.remove_edge(u,w)
    return G
  def graphviz(self):
    """ Return a graphviz string describing the graph and its labels """
    gv = 'digraph {\n'
    indices = { v : str(k) for k,v in enumerate(self.vertices())}
    for v in self.vertices(): gv += indices[v] + '[label="' + self.vertex_label(v) + '"];\n'
    for (u,v) in self.edges(): gv += indices[u] + ' -> ' + indices[v] + ' [label="' + self.edge_label(u,v) + '"];\n'
    return gv + '}\n'

###############################
# Partially Ordered Set Class #
###############################

class Poset:
  def __init__(self, graph):
    """  Create a Poset P from a DAG G such that x <= y in P iff there is a path from x to y in G """
    self.vertices_ = set(graph.vertices())
    self.descendants_ = graph.transitive_closure()
    self.ancestors_ = self.descendants_.transpose()
    self.children_ = graph.transitive_reduction()
    self.parents_ = self.children_.transpose()
  def vertices(self):
    """ Return the set of elements in the poset """
    return self.vertices_
  def parents(self, v):
    """ Return the immediate predecessors of v in the poset """
    return self.parents_.adjacencies(v)
  def children(self, v):
    """ Return the immediate successors of v in the poset """
    return self.children_.adjacencies(v)
  def ancestors(self, v):
    """ Return the set { u : u < v } """
    return self.ancestors_.adjacencies(v)
  def descendants(self, v):
    """ Return the set { u : v < u } """
    return self.descendants_.adjacencies(v)
  def less(self, u, v):
    """ Return True if u < v, False otherwise """
    return u in self.ancestors(v)
  def maximal(self, subset):
    """ Return the set of elements in "subset" which are maximal """
    return frozenset({ u for u in subset if not any ( self.less(u,v) for v in subset ) })

###################################
# Poset-to-PatternGraph Algorithm #
###################################

def PosetToPatternGraph(poset):
  """ Generate from poset the Hasse diagram of the poset of down-sets of "poset" ordered by inclusion """
  pattern_graph = Graph()
  recursion_stack = [poset.maximal(poset.vertices())]
  while recursion_stack:
    clique = recursion_stack.pop()
    for v in clique:
      parent_clique = poset.maximal(clique.difference([v]).union(poset.parents(v)))
      if parent_clique not in pattern_graph.vertices():
        recursion_stack.append (parent_clique)
      pattern_graph.add_edge (parent_clique,clique, str(v))
  return pattern_graph

###################################
# Create Interval Graph from Data #
###################################

def IntervalGraph(events):
  """ Generate a labelled interval graph from a list of "events". 
      An "event" is a two-element list [label, interval],
      where interval is two-element list of floats. e.g. ["X", [1.0, 3.4] ]
      Each event becomes a vertex in the graph and there is an edge
      from u to v whenever the interval associated with u is disjoint
      from and bounded above by the interval associated with v
  """
  N = len(events)
  G = Graph()
  for i in range(0, N):
    G.add_vertex(i,events[i][0])
  for i in range(0, N):
    for j in range(0, N):
      if events[i][1][1] <= events[j][1][0]:
        G.add_edge(i,j)
  return G

# ##############
# # 7D Example #
# ##############

# # Names of genes in order given in network spec file
# genes = ["FKH1", "SPT2", "PLM2", "WTM2", "SWI4", "NDD1", "HCM1" ]

# # Location of extreme
# data = {}
# data["FKH1"] = [[10.0,20.0], [40.0,50.0], [80.0,90.0], [115.0,125.0]]    #MIN
# data["SPT2"] = [[10.0,90.0], [95.0, 150.0]]                              #MIN
# data["PLM2"] = [[10.0,15.0], [25.0,35.0], [70.0,80.0], [95.0,105.0]]     #MIN
# data["WTM2"] = [[25.0,35.0], [65.0,75.0], [100.0,110.0], [145.0,150.0]]  #MIN
# data["SWI4"] = [[15.0,25.0], [50.0,60.0], [120.0,130.0], [145.0,150.0]]  #MAX
# data["NDD1"] = [[10.0,15.0], [35.0,45.0], [75.0,85.0], [105.0,115.0]]    #MIN
# data["HCM1"] = [[20.0,30.0], [60.0,70.0], [95.0,105.0], [140.0,150.0] ]  #MAX

# # List of min/max events
# events = [ [gene, datum] for gene in genes for datum in data[gene] ]

# # Final state label:
# #  2*D bits, the ith bit is 1 if the final state is decreasing in variable i
# #        the (i+D)th bit is 1 if the final state is increasing in variable i
# label = 0b10100000101111   # variables 0,1,2,3,5 decreasing, variables 4,6 increasing

# # Create the interval graph
# G = IntervalGraph(events)

# # Produce a JSON file suitable for the C++ pattern matching program
# HasseDiagram = G.transitive_reduction();
# output = {}
# output["poset"] = [ list(HasseDiagram.adjacencies(i)) for i in HasseDiagram.vertices() ]
# output["events"] = [ genes.index(HasseDiagram.vertex_label(i)) for i in HasseDiagram.vertices() ]
# output["label"] = label
# output["dimension"] = len(genes)
# with open('pattern.json', 'w') as fp:
#   json.dump(output, fp)

# # Output in graphviz format
# #print HasseDiagram.graphviz()

# P = Poset(HasseDiagram)
# PG = PosetToPatternGraph(P)
# print PG.graphviz()

