//  graph.h <Starter Code>
//  < Richard Luong >
//
//  Basic graph class using adjacency list representation.
//
//  University of Illinois at Chicago
//  CS 251: Fall 2021
//  Project #7 - Openstreet Maps
//

#pragma once

#include <unordered_map>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
	private:
		//
		//  Adjacency list implementation where each key in a map represents
		//  a vertice and each verticess maps to a set of vertices that are
		//  its neighbors.
		//
		unordered_map<VertexT, unordered_map<VertexT, WeightT>> adjList;
		int numEdges;

	public:
	//
	//  constructor:
	//
	//  Constructs an empty graph.
	//
	graph() {
		numEdges = 0;
	}

	//
	//  copy constructor
	//
	graph(const graph& other) {
		this->adjList = other.adjList;
		this->numEdges = other.numEdges;
	}

	//
	//  operator =
	//
	graph& operator=(const graph& other) {
		if (this == &other)
			return *this;

		this->adjList.clear();
		this->adjList = other.adjList;
		this->numEdges = other.numEdges;

		return *this;
	}

	//
	//  NumVertices
	//
	//  Returns the # of vertices currently in the graph.
	//
	int NumVertices() const {
		return adjList.size();
	}

	//
	//  NumEdges
	//
	//  Returns the # of edges currently in the graph.
	//
	int NumEdges() const {
		return numEdges;
	}

	//
	//  addVertex
	//
	//  Adds the vertex v to the graph if there's room, and if so
	//  returns true.  If the vertex already
	//  exists in the graph, then false is returned.
	//
	bool addVertex(VertexT v) {;
		//
		//  is the vertex already in the graph?  If so, we do not
		//  insert again otherwise Vertices may fill with duplicates:
		//
		if (adjList.count(v))
			return false;

		//
		//  if we get here, vertex does not exist so insert.
		//
		adjList[v];

		return true;
	}

	//
	//  addEdge
	//
	//  Adds the edge (from, to, weight) to the graph, and returns
	//  true.  If the vertices do not exist false is returned.
	//
	//  NOTE: if the edge already exists, the existing edge weight
	//  is overwritten with the new edge weight.
	//
	bool addEdge(VertexT from, VertexT to, WeightT weight) {
		//
		//  1. from vertice does not exist
		//  2. to vertice does not exist
		//
		if (!adjList.count(from) || !adjList.count(to))
			return false;

		//
		// from and to vertices exist so we're adding
		// a new edge
		//
		numEdges += !adjList[from].count(to);
		adjList[from][to] = weight;

		return true;
	}

	//
	// getWeight
	//
	// Returns the weight associated with a given edge.  If
	// the edge exists, the weight is returned via the reference
	// parameter and true is returned.  If the edge does not
	// exist, the weight parameter is unchanged and false is
	// returned.
	//
	bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
		//
		//  1. from vertice does not exist
		//  2. to vertice does not exist
		//  3. edge does not exist
		//
		if (!adjList.count(from) || !adjList.count(to) || !adjList.at(from).count(to))
			return false;

		//
		//  Get weight
		//
		weight = adjList.at(from).at(to);

		return true;
	}

	//
	// neighbors
	//
	// Returns a set containing the neighbors of v, i.e. all
	// vertices that can be reached from v along one edge.
	// Since a set is returned, the neighbors are returned in
	// sorted order; use foreach to iterate through the set.
	//
	set<VertexT> neighbors(VertexT v) const {		
		set<VertexT> n;

		if (adjList.count(v)) {
			for (auto adj : adjList.at(v))
				n.insert(adj.first);
		}
		return n;
	}

	//
	// getVertices
	//
	// Returns a vector containing all the vertices currently in
	// the graph.
	//
	vector<VertexT> getVertices() const {
		vector<VertexT> keys;

		for (auto k : adjList)
			keys.push_back(k.first);

		return keys;
	}

	//
	// dump
	//
	// Dumps the internal state of the graph for debugging purposes.
	//
	// Example:
	//    graph<string,int>  G(26);
	//    ...
	//    G.dump(cout);  // dump to console
	//
	void dump(ostream& output) const {
		int i = 0;

		output << "***************************************************" << endl;
		output << "********************* GRAPH ***********************" << endl;

		output << "**Num vertices: " << this->NumVertices() << endl;
		output << "**Num edges: " << this->NumEdges() << endl;

		output << endl;
		output << "**Vertices:" << endl;
		for (auto keys : adjList) {
		  output << " " << i << ". " << keys.first << endl;
		  ++i;
		}

		output << endl;
		output << "**Edges:" << endl;
		for (auto keys : adjList) {
			output << " " << keys.first << ": ";

			for (auto neighbors : keys.second) {
				output << "(" << neighbors.first << "," << neighbors.second << ") ";
			}
			output << endl;
		}

		output << "**************************************************" << endl;
	}
};