// All required headers are included in headerfile
#include "headerFile.h"


// Create a shortcut for new line.
#define newline std::cout << "\n"


using namespace boost;

// Function to print the graph along with edge weights
void printGraph(const Graph& G, const WeightMap& weightMap)
{
    // Iterate over all vertices in the graph
    VertexIterator vi,viend ;
    for (tie(vi,viend)=vertices(G); vi!=viend;++vi)
    {
        std::cout << "Vertex " << *vi << ": ";            
        // Iterate over all outgoing edges from the current vertex
        OutEdgeIterator oei, oeiend;
        for (tie(oei,oeiend)=out_edges(*vi,G); oei!=oeiend; ++oei)
            std::cout << "Edge:" << *oei << " Weight: " << weightMap[*oei] << "  ";
        std::cout << std::endl;

    }

}

// Function to get the opposite vertex of a given edge and vertex
Vertex opposite(const Graph& G, Edge e , Vertex u){
    if(source(e, G) == u){
        return target(e, G);
    }
    else{
        return source(e,G);
    }
}

// Implementation of Dijkstra's algorithm
void dijkstra(const Graph& G, Vertex s, const WeightMap& cost, std::vector<int>& dist, std::vector<Edge>& pred, bool reverse){
    PriorityQueue minHeap;
    handle_t handleType;
    unordered_map<Vertex, handle_t> handleMap;
    
    EdgeIterator ei, ei_End;
    tie(ei, ei_End) = edges(G);
    Edge nillEdge = *ei_End;

    // Initialize the predecessors for all edges
    for(tie(ei, ei_End) = edges(G); ei != ei_End; ++ei){
        pred.push_back(nillEdge);
    }
    Vertex v;
    // Initialize the distance for the source vertex
    dist[int(s)] = 0;
    std::pair<Vertex, int> pair;
    pair.first = s;
    pair.second  = 0;

    handleType = minHeap.push(pair);
    handleMap[pair.first] = handleType;
    // handleMap.insert(pair.first, handleType);  


    
    while(!minHeap.empty()){
        Vertex u = minHeap.top().first;
        minHeap.pop();

        int du = dist[int(u)];

        if(!reverse){
            OutEdgeIterator out_edge_begin, out_edge_end;
            // Iterate over all outgoing edges from the current vertex
            for(tie(out_edge_begin, out_edge_end) = out_edges(u, G); out_edge_begin != out_edge_end; ++out_edge_begin){
                
                v = opposite(G, *out_edge_begin, u);   
                int c = du + cost[*out_edge_begin];

                if(pred[int(v)] == nillEdge && v != s){
                    pair.first = v;
                    pair.second = c;
                    handleType = minHeap.push(pair);
                    handleMap[pair.first] = handleType;
                }
                else if(c < dist[int(v)]){
                    pair.first = v;
                    pair.second = c;
                    handle_t retrievedHandle;
                    if (handleMap.find(v) != handleMap.end()) {
                        retrievedHandle = handleMap[v];
                    }

                    minHeap.decrease(retrievedHandle, pair);
                }
                else{
                    continue;
                }

                dist[int(v)] = c;
                pred[int(v)] = *out_edge_begin;
            }

        }
        else{
            InEdgeIterator in_edge_begin, in_edge_end;
            // Iterate over all incoming edges to the current vertex
            for(tie(in_edge_begin, in_edge_end) = in_edges(u, G); in_edge_begin != in_edge_end; ++in_edge_begin){
                
                v = opposite(G, *in_edge_begin, u);   
                int c = du + cost[*in_edge_begin];

                if(pred[int(v)] == nillEdge && v != s){
                    pair.first = v;
                    pair.second = c;
                    handleType = minHeap.push(pair);
                    handleMap[pair.first] = handleType;
                }
                else if(c < dist[int(v)]){
                    pair.first = v;
                    pair.second = c;
                    handle_t retrievedHandle;
                    if (handleMap.find(v) != handleMap.end()) {
                        // std::cout << "test1" << std::endl;

                        retrievedHandle = handleMap[v];
                    }

                    minHeap.decrease(retrievedHandle, pair);
                    // std::cout << "test2" << std::endl;
                }
                else{
                    continue;
                }

                dist[int(v)] = c;
                pred[int(v)] = *in_edge_begin;
            }

        }

        
    }
    handleMap.clear();
    
}
// Function to find the shortest path using Dijkstra's algorithm from source to target
void dijkstraSP(const Graph& G, Vertex s, Vertex t, WeightMap& cost, std::vector<int>& dist, std::vector<Edge>& pred, std::vector<Vertex>& result){
    PriorityQueue minHeap;
    handle_t handleType;
    unordered_map<Vertex, handle_t> handleMap;
    
    EdgeIterator ei, ei_End;
    tie(ei, ei_End) = edges(G);
    Edge nillEdge = *ei_End;

    for(tie(ei, ei_End) = edges(G); ei != ei_End; ++ei){
        pred.push_back(nillEdge);
    }


    Vertex v;
    dist[int(s)] = 0;
    std::pair<Vertex, int> pair;
    pair.first = s;
    pair.second  = 0;

    handleType = minHeap.push(pair);
    handleMap[pair.first] = handleType;
    
    while(!minHeap.empty()){
        Vertex u = minHeap.top().first;
        minHeap.pop();
        result.push_back(u);
        if(u == t){
            // std::cout << "Found the node: " << u << std::endl;
            break;
        }
        int du = dist[int(u)];
        OutEdgeIterator out_edge_begin, out_edge_end;

        for(tie(out_edge_begin, out_edge_end) = out_edges(u, G); out_edge_begin != out_edge_end; ++out_edge_begin){
            
            v = opposite(G, *out_edge_begin, u);   
            int c = du + cost[*out_edge_begin];

            if(pred[int(v)] == nillEdge && v != s){
                pair.first = v;
                pair.second = c;
                handleType = minHeap.push(pair);
                handleMap[pair.first] = handleType;
            }
            else if(c < dist[int(v)]){
                continue;
                pair.first = v;
                pair.second = c;
                handle_t retrievedHandle;
                if (handleMap.find(v) != handleMap.end()) {
                    retrievedHandle = handleMap[v];
                }
                minHeap.decrease(retrievedHandle, pair);
            }
            else{
                continue;
            }

            dist[int(v)] = c;
            pred[int(v)] = *out_edge_begin;
            
        }
    }
    handleMap.clear();    
}
// Function to create a grid graph with specified rows and columns
Graph createGrid(int rows, int cols){
    Graph G;
    
    std::vector<Vertex> nodes;


    for(int i = 0; i<cols*rows; ++i){
        nodes.push_back(add_vertex(G));
    }
    

    for(int node=0; node<nodes.size(); ++node){
        
        if(node<nodes.size()-cols){
            add_edge(nodes[node], nodes[node+cols],G);
            add_edge(nodes[node+cols], nodes[node],G);
        }

        if(node<nodes.size()-1){
            if(!(node%cols == cols-1)){
                add_edge(nodes[node], nodes[node+1],G);
                add_edge(nodes[node+1], nodes[node],G);
            }
            

            
        }
    }
    
    nodes.clear();
    return G;
}

// Function to calculate the Euclidean distance between two points  
int distance(int x1, int y1, int x2, int y2) {
    return std::sqrt(std::pow((x2 - x1), 2) + std::pow((y2 - y1), 2));
}

// Function to recalculate weights based on the target vertex and grid layout
void recalculateWeight(const Graph& G, Vertex t, const WeightMap& cost, int cols){

    std::vector<int> tarDist(num_vertices(G),-1);

    VertexIterator v_current,v_end ;
    EdgeIterator ei, ei_end;
    Vertex next_node;

    int tarNode = int(t);
    // Calculate new weights based on the Euclidean distance to the target vertex
    for (tie(v_current,v_end)=vertices(G); v_current!=v_end;++v_current)
    {
        int current_node = int(*v_current);
        tarDist[current_node] = distance(current_node%cols, current_node/cols, tarNode%cols, tarNode/cols);
        OutEdgeIterator out_edge_begin, out_edge_end;
        for(tie(out_edge_begin, out_edge_end) = out_edges(*v_current, G); out_edge_begin != out_edge_end; ++out_edge_begin){
            next_node = opposite(G, *out_edge_begin, *v_current);
            if(tarDist[next_node] == -1){
                tarDist[int(next_node)] = distance(int(next_node)%cols, int(next_node)/cols, tarNode%cols, tarNode/cols);                
                // std::cout<< "edge: "<< *out_edge_begin << "\tcost " << cost[*out_edge_begin] << "\tj dist: " << tarDist[next_node] << "\t i dist: " << tarDist[current_node] << std::endl; 
                int new_cost = cost[*out_edge_begin] + tarDist[next_node] - tarDist[current_node];
                cost[*out_edge_begin] = new_cost;
            }
        }

    }

    tarDist.clear();

} 
// Function to calculate l1 and l2 dists
void calculate_L1_L2(const Graph& G, const WeightMap& cost, std::vector<int>& distL1, std::vector<int>& distL1rev,  std::vector<int>& distL2,  std::vector<int>& distL2rev){
    std::vector<Edge> pred1;
    
    Vertex source_L1, source_L2;
    VertexIterator vI,vIend;
    srand(time(0));
    int ran = rand() % num_vertices(G);
    tie(vI,vIend)=vertices(G);
    source_L1 = *vI+ran;  
    
    dijkstra(G, source_L1, cost, distL1, pred1, 0);
    pred1.clear();
    dijkstra(G, source_L1, cost, distL1rev, pred1, 1);
    pred1.clear();
    int maxElementIndex = std::max_element(distL1rev.begin(),distL1rev.end()) - distL1rev.begin();
    

    tie(vI,vIend)=vertices(G);
    source_L2 = *vI+maxElementIndex;
  
    dijkstra(G, source_L2, cost, distL2, pred1, 0);
    pred1.clear();
    dijkstra(G, source_L2, cost, distL2rev, pred1, 1);
    pred1.clear();

}

// Return max to calculate ht(i)
int returnMax(std::vector<int>& distL1, std::vector<int>& distL1rev,  std::vector<int>& distL2,  std::vector<int>& distL2rev, Vertex sI, Vertex tI){

    int a = distL1[int(tI)] - distL1[int(sI)];
    int b = distL2rev[int(sI)] - distL2rev[int(tI)];    
    return std::max(a,b);

}


//reculculate the weights using the formula of w'
void recalculateWeight_L1_L2(const Graph& G, Vertex t, const WeightMap& cost){
    
    std::vector<int> distL1(num_vertices(G));
    std::vector<int> distL1rev(num_vertices(G));
    std::vector<int> distL2(num_vertices(G));
    std::vector<int> distL2rev(num_vertices(G));

    calculate_L1_L2(G,cost,distL1,distL1rev,distL2,distL2rev);

    VertexIterator v_current,v_end ;
    EdgeIterator ei, ei_end;
    Vertex next_node;
    std::vector<int> tarDist(num_vertices(G),-1);

    int tarNode = int(t);
    for (tie(v_current,v_end)=vertices(G); v_current!=v_end;++v_current)
    {
        int current_node = int(*v_current);
        OutEdgeIterator out_edge_begin, out_edge_end;
        for(tie(out_edge_begin, out_edge_end) = out_edges(*v_current, G); out_edge_begin != out_edge_end; ++out_edge_begin){
            next_node = opposite(G, *out_edge_begin, *v_current);
            
            int ht_curr = returnMax(distL1,distL1rev,distL2,distL2rev,*v_current,t);
            tarDist[int(current_node)] = ht_curr;
            
            if(tarDist[int(next_node)] == -1){
                int ht_next = returnMax(distL1,distL1rev,distL2,distL2rev,next_node,t);
                tarDist[int(next_node)] = ht_next;
                int new_cost = cost[*out_edge_begin] + ht_next - ht_curr;
                cost[*out_edge_begin] = new_cost;
            
            }
            

        }

    }
    
}

int main()
{   

    // Clear terminal window.
    system("clear");

    // Graph G1;
    // WeightMap weightMap1;

    // int rows = 80;
    // int cols = 1000;
    
    // G1 = createGrid(rows,cols);
    // std::cout << "Grid Graph" << std::endl;
    // std::cout << "Rows: " << rows << "\tCols: " << cols << std::endl;
    // std::cout << "Verticles " << num_vertices(G1) << "\tEdges " << num_edges(G1) << std::endl;

    // newline;

    EdgeIterator ei, eiEnd;        
    // srand(time(0)); // Use current time as seed for random generator 
    // for (tie(ei, eiEnd) = edges(G1); ei != eiEnd; ++ei) weightMap1[*ei] = (rand() % 10000)+1;

    // printGraph(G1,weightMap1);
    
    // newline;
    // newline;

    Graph G2;
    WeightMap weightMap2;

    random::mt19937 rng;
    rng.seed(static_cast<unsigned int>(std::time(0))); // Seed the generator with the current time

    //generation  of random graph in boost 
    generate_random_graph(G2, 60000, 120000 ,rng, false, true);

    std::cout << "Random Generated Graph" << std::endl;
    std::cout << "Num of verticles " << num_vertices(G2) << std::endl;
    std::cout << "Num of edges " << num_edges(G2) << std::endl;

    newline;
    srand(time(0)); // Use current time as seed for random generator 
    for (tie(ei, eiEnd) = edges(G2); ei != eiEnd; ++ei) weightMap2[*ei] = (rand() % 10000)+1;

    // printGraph(G2,weightMap);


    //store the vertex of first and last col to each vector in order to use to to find a random vertex
    // remove comments to work

    // std::vector<Vertex> firstCol;
    // std::vector<Vertex> lastCol;

    // VertexIterator vBegin, vEnd;
    // tie(vBegin, vEnd) = vertices(G1);

    // for(tie(vBegin, vEnd); vBegin != vEnd; ++vBegin){

    //     int k = *vBegin;
    //     if((k%cols) == 0){
    //         firstCol.push_back(*vBegin);
    //     }
    //     if((k%cols) == cols-1){
    //         lastCol.push_back(*vBegin);
    //     }
    // }

    srand(time(0));

    // random Vertexs for grid //remove comments to work

    // Vertex randSource = firstCol[rand()%firstCol.size()];
    // Vertex randTarget = lastCol[rand()%lastCol.size()];

    VertexIterator iterV, iterVend;
    tie(iterV, iterVend) = vertices(G2);

    Vertex randSource = *iterV + rand()%num_vertices(G2); 
    Vertex randTarget = *iterV + rand()%num_vertices(G2);

    


    //initialize the vector with the number of nodes
    std::vector<int> dist1(num_vertices(G2));
    std::vector<Edge> pred1;
    // store the vertex which poped from heap
    std::vector<Vertex> result1;

    long counter;
    float secondsRes1;
    float secondsRes2;
    clock_t start, finish;
    int mulNum = 10;    
    
    // Call my dijkstraSP checker and the timer calculate its runtime.
    start = clock();
    counter = 0;
    while(clock() - start < mulNum*CLOCKS_PER_SEC){
        counter+=1;
        dist1.clear();
        dist1.resize(num_vertices(G2));
        pred1.clear();
        result1.clear();
        dijkstraSP(G2, randSource, randTarget, weightMap2, dist1, pred1, result1);
    }
    finish = clock();
    secondsRes1 = float(finish - start) / CLOCKS_PER_SEC;
    secondsRes1 = secondsRes1 / counter;


    //print the results
    std::cout << "DijkstraSP runtime: " << std::fixed << secondsRes1 << std::setprecision(9);
    std::cout << " seconds" << std::endl;

    std::cout << "Number of elements extracted from last heap: " << result1.size() << std::endl;
    newline;

    //recalculateWeight(G1, randTarget, weightMap1, cols);
    recalculateWeight_L1_L2(G2, randTarget, weightMap2);

    std::vector<int> dist2(num_vertices(G2));
    std::vector<Edge> pred2;
    std::vector<Vertex> result2;

   
    // Call my dijkstraSP checker and the timer calculate its runtime.
    start = clock();
    counter = 0;
    while(clock() - start < mulNum*CLOCKS_PER_SEC){
        counter+=1;
        dist2.clear();
        dist2.resize(num_vertices(G2));
        pred2.clear();
        result2.clear();
        dijkstraSP(G2, randSource, randTarget, weightMap2, dist2, pred2, result2);
    }
    finish = clock();
    secondsRes2 = float(finish - start) / CLOCKS_PER_SEC;
    secondsRes2 = secondsRes2 / counter;

    std::cout << "DijkstraSP runtime: " << std::fixed << secondsRes2 << std::setprecision(9);
    std::cout << " seconds" << std::endl;

    std::cout << "Number of elements extracted from last heap: " << result2.size() << std::endl;
    newline;
    
    // recalculateWeight(G,*iter2-1, weightMap, cols);

    
    // std::vector<int> dist2(num_vertices(G));
    // std::vector<Edge> pred2;
    // std::vector<Vertex> result2;
    

    // dijkstraSP(G, *iter1, *iter2-1, weightMap, dist2, pred2, result2);
    // newline;

    // printGraph(G,weightMap);


    // std::cout << "#####RESULTS#####" << std::endl;
    // std::cout << "res1: " << result1.size() << "\tres2: " << result2.size() << std::endl; 

    return 0;
}