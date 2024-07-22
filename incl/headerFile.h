#ifndef HEADERFILE_H
#define HEADERFILE_H


#include <boost/graph/random.hpp>
#include <boost/random.hpp>

//included by developer

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_map> 
#include <cmath>
#include <ctime>
#include <bits/stdc++.h>


#include <boost/graph/adjacency_list.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/unordered_map.hpp>

/* Edge Property  (all properties: https://www.boost.org/doc/libs/1_85_0_beta1/libs/graph/doc/property.html) */
typedef boost::property<boost::edge_weight_t, int> EdgeProperty;

/* Graph (https://www.boost.org/doc/libs/1_61_0/libs/graph/doc/graph_concepts.html) */
typedef boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeProperty> Graph;

/* Vertex descriptor */
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

/* Edge descriptor */
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

/* Iterators */
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;
typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIterator;

/* Create Edge Property Map */
typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

/* Fibonacci Heap */
struct CompareVertices 
{
    bool operator()(const std::pair<Vertex, int>& p1, const std::pair<Vertex, int>& p2) const
    {
        return p1.second > p2.second; 
    }
};

typedef boost::heap::fibonacci_heap<std::pair<Vertex, int>, boost::heap::compare<CompareVertices> > PriorityQueue;

/* Handle for the fibonacci heap. You should use it to decrease the priority of a pair in the PriorityQueue */
typedef PriorityQueue::handle_type handle_t;

#endif