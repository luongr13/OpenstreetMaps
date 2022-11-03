// application.cpp <Starter Code>
// < Richard Luong >
//
// Interface for Openstreet Maps program that provides directions for
// the shortest path to a destination that is in the middle
// of two locations.
//
// University of Illinois at Chicago
// CS 251: Fall 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> /*setprecision*/
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <vector>

#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

//
//  functor for priority_queue
//
//  Smaller distances have priority
//
class prioritize {
 public:
  bool operator()(const pair<long long, double> p1,
                  const pair<long long, double> p2) const {
    return p1.second > p2.second;
  }
};

//
// Implement your creative component application here
// TO DO: add arguments
//
void creative() {}

//  _getEdges
//
//  Calculates edge weights between vertices and inserts edge.
//
void getEdges(vector<FootwayInfo>& Footways, map<long long, Coordinates>& Nodes,
              graph<long long, double>& G) {
  for (auto footway : Footways) {
    for (long unsigned int i = 0; i < footway.Nodes.size() - 1; ++i) {
      //
      //  get IDs
      //
      long long id1 = footway.Nodes.at(i);
      long long id2 = footway.Nodes.at(i + 1);

      //
      //  get coordinates of each id
      //
      Coordinates coord1 = Nodes.at(id1);
      Coordinates coord2 = Nodes.at(id2);

      double distance =
          distBetween2Points(coord1.Lat, coord1.Lon, coord2.Lat, coord2.Lon);

      //
      // add 2 way edge
      //
      if (!G.addEdge(id1, id2, distance)) cout << "Unable to add edge" << endl;
      if (!G.addEdge(id2, id1, distance)) cout << "Unable to add edge" << endl;
    }
  }
}

//  updateGraph
//
//  Adds Nodes as vertices and Footways as edges.
//
void updateGraph(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, graph<long long, double>& G) {
  //
  //  Add vertices.
  //
  for (auto key : Nodes) G.addVertex(key.first);

  //
  //  Add edges.
  //
  getEdges(Footways, Nodes, G);
}

//  printStats
//
//  Prints various stats across multiple structures.
//
void printStats(const map<long long, Coordinates>& Nodes,
                const vector<FootwayInfo>& Footways,
                const vector<BuildingInfo>& Buildings,
                const graph<long long, double>& G) {
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
}

//  findBuildings
//
//  Takes two string queries and searches the vector of buildings for
//  buildings that match the queries.
//
//  1. If an abbreviation is used, it must be exact.
//  2. If the full name is used, a substring will work.
//
void findBuildings(const vector<BuildingInfo>& Buildings, string query1,
                   string query2, BuildingInfo& person1,
                   BuildingInfo& person2) {
  size_t found;
  bool p1found = 0, p2found = 0;
  // bool p2found = 0;
  for (auto b : Buildings) {
    found = b.Fullname.find(query1);
    if (!p1found && (b.Abbrev == query1 || found != string::npos)) {
      person1 = b;
      p1found = 1;
    }

    found = b.Fullname.find(query2);
    if (!p2found && (b.Abbrev == query2 || found != string::npos)) {
      person2 = b;
      p2found = 1;
    }
  }
}

//  nearestBuilding
//
//  Finds the nearest building based on coordinates given.
//  Used to find the nearest building to the midpoint coordinates.
//
BuildingInfo nearestBuilding(const vector<BuildingInfo>& Buildings,
                             const set<string> unreachableBuildings,
                             const Coordinates center) {
  BuildingInfo ret;
  double min = numeric_limits<double>::max();

  for (auto building : Buildings) {
    Coordinates xy = building.Coords;
    double distance =
        distBetween2Points(center.Lat, center.Lon, xy.Lat, xy.Lon);
    if (distance < min && !unreachableBuildings.count(building.Fullname)) {
      min = distance;
      ret = building;
    }
  }
  return ret;
}

//  nearestNode
//
//  Finds the nearest node to a building based on the building's coordinates.
//  Return the building's ID.
//
long long nearestNode(map<long long, Coordinates>& Nodes,
                      vector<FootwayInfo>& Footways, const BuildingInfo b) {
  long long id;
  double distance;
  double min = numeric_limits<double>::max();

  for (auto footway : Footways) {
    for (auto ID : footway.Nodes) {
      distance = distBetween2Points(b.Coords.Lat, b.Coords.Lon, Nodes[ID].Lat,
                                    Nodes[ID].Lon);
      if (distance < min) {
        min = distance;
        id = ID;
      }
    }
  }

  cout << " " << id << endl;
  cout << " (" << Nodes[id].Lat << ", " << Nodes[id].Lon << ")";

  return id;
}

//	Dijkstra
//
//	Generates a map with shortest distance from starting vertex to
//	every other vertex in the graph.
//
void Dijkstra(graph<long long, double>& G, long long startV,
              map<long long, double>& distances,
              map<long long, long long>& predecessors) {
  double INF = numeric_limits<double>::max();
  priority_queue<pair<long long, double>,          //  (key,value) pair
                 vector<pair<long long, double>>,  //  stores pair into vector
                 prioritize>
      unvisited_q;
  set<long long> visited;

  for (long long vertex : G.getVertices()) {
    unvisited_q.push({vertex, INF});
    distances[vertex] = INF;
  }

  distances[startV] = 0;
  unvisited_q.push({startV, 0});

  while (unvisited_q.size() != 0) {
    pair<long long, double> vertex = unvisited_q.top();
    unvisited_q.pop();

    long long id = vertex.first;

    if (distances.at(id) == INF)
      break;
    else if (visited.count(id))
      continue;
    else
      visited.insert(id);

    for (auto neighbor : G.neighbors(id)) {
      double weight;
      G.getWeight(id, neighbor, weight);

      double altDistance = distances[id] + weight;

      if (altDistance < distances[neighbor]) {
        distances[neighbor] = altDistance;
        predecessors[neighbor] = id;
        unvisited_q.push({neighbor, altDistance});
      }
    }
  }
}

//  getPath
//
//  Given the predecessor map and am end vertex,
//  print the path to the end vertex in the correct order
//  using its predecessors.
//
vector<long long> getPath(map<long long, long long>& predecessors,
                          long long endVertex) {
  stack<long long> path;
  vector<long long> inOrderPath;
  long long currV = endVertex;

  while (currV != 0) {
    path.push(currV);
    currV = predecessors[currV];
  }
  while (path.size() != 0) {
    currV = path.top();
    path.pop();
    inOrderPath.push_back(currV);
  }
  return inOrderPath;
}

//  printPaths
//
//  Prints the path for both persons
//
void printPaths(map<long long, double>& p1Distance,
                map<long long, double>& p2Distance, vector<long long>& p1Path,
                vector<long long>& p2Path, long long node_nearby_center) {
  cout << endl;
  cout << "Person 1's distance to dest: " << p1Distance[node_nearby_center]
       << " miles" << endl;
  cout << "Path: ";
  for (long unsigned int i = 0; i < p1Path.size(); ++i) {
    cout << p1Path[i];
    if (i != p1Path.size() - 1) cout << "->";
  }
  cout << endl;

  cout << endl;

  cout << "Person 2's distance to dest: " << p2Distance[node_nearby_center]
       << " miles" << endl;
  cout << "Path: ";
  for (long unsigned int i = 0; i < p2Path.size(); ++i) {
    cout << p2Path[i];
    if (i != p2Path.size() - 1) cout << "->";
  }
}

//  printBuildingInfo
//
//  Prints a building's full name and it's coordinates
//
void printBuildingInfo(BuildingInfo building, Coordinates coord) {
  cout << " " << building.Fullname << endl;
  cout << " (" << coord.Lat << ", " << coord.Lon << ")";
}

void application(map<long long, Coordinates>& Nodes,
                 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
                 graph<long long, double>& G) {
  double INF = numeric_limits<double>::max();  //	Max double value
  string person1Building, person2Building;     //	Queries

  cout << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);

  while (person1Building != "#") {
    Coordinates coord1, coord2, center;                //	Coordinates
    BuildingInfo building1, building2, nearby_center;  //	Buildings

    map<long long, double> p1Distance, p2Distance;  //	Distance Maps
    map<long long, long long> p1Pred, p2Pred;       //	Predecessor Maps

    set<string> unreachableBuildings;

    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    cout << endl;
    //
    //  Getting each person's Building Info.
    //
    findBuildings(Buildings, person1Building, person2Building, building1,
                  building2);

    //
    //  Check if the buildings were found based on input query.
    //  Return false to reloop if not found.
    //
    if (building1.Fullname == "") {
      cout << "Person 1's building not found" << endl;
    } else if (building2.Fullname == "") {
      cout << "Person 2's building not found" << endl;
    } else {
      //
      //	Find the middle point between two building coordinates.
      //
      Coordinates coord1 = building1.Coords;
      Coordinates coord2 = building2.Coords;

      bool foundDestination = 0, repeatLoop = 0;
      while (!foundDestination) {
        Coordinates center = centerBetween2Points(coord1.Lat, coord1.Lon,
                                                  coord2.Lat, coord2.Lon);
        //
        //	Find the closest building to the midpoint / center coordinates.
        //
        nearby_center =
            nearestBuilding(Buildings, unreachableBuildings, center);

        //
        //	Print Building Information
        //
        if (!repeatLoop) {
          cout << "Person 1's point:" << endl;
          printBuildingInfo(building1, coord1);
          cout << endl;
          cout << "Person 2's point:" << endl;
          printBuildingInfo(building2, coord2);
          cout << endl;
          cout << "Destination Building:" << endl;
          printBuildingInfo(nearby_center, nearby_center.Coords);
          cout << endl;
          cout << endl;
        } else {
          cout << "New destination building:" << endl;
          printBuildingInfo(nearby_center, nearby_center.Coords);
          cout << endl;
        }

        //
        //	Since building might not be on Footways, find the nearest node
        //	for each location since nodes are on footways and print
        //	node ID with it's coordinates
        //
        long long node_building1, node_building2;
        if (!repeatLoop) {
          cout << "Nearest P1 node:" << endl;
          node_building1 = nearestNode(Nodes, Footways, building1);
          cout << endl;
          cout << "Nearest P2 node:" << endl;
          node_building2 = nearestNode(Nodes, Footways, building2);
          cout << endl;
        }
        cout << "Nearest destination node:" << endl;
        long long node_nearby_center =
            nearestNode(Nodes, Footways, nearby_center);

        //  MILESTONE 10
        Dijkstra(G, node_building1, p1Distance, p1Pred);
        Dijkstra(G, node_building2, p2Distance, p2Pred);

        vector<long long> p1Path = getPath(p1Pred, node_nearby_center);
        vector<long long> p2Path = getPath(p2Pred, node_nearby_center);

        //  MILESTONE 11
        if (p1Distance[node_building2] == INF) {
          cout << endl;
          cout << endl;
          cout << "Sorry, destination unreachable." << endl;
          break;
        } else if (p1Distance[node_nearby_center] == INF ||
                   p2Distance[node_nearby_center] == INF) {
          cout << endl;
          cout << endl;
          cout << "At least one person was unable to reach the destination "
                  "building. ";
          cout << "Finding next closest building..." << endl;
          cout << endl;
          unreachableBuildings.insert(nearby_center.Fullname);
          repeatLoop = 1;
        } else {
          //  PRINT PATHS
          cout << endl;
          printPaths(p1Distance, p2Distance, p1Path, p2Path,
                     node_nearby_center);
          cout << endl;
          foundDestination = 1;
        }
      }
    }
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates> Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo> Footways;
  // info about each building, in no particular order
  vector<BuildingInfo> Buildings;
  XMLDocument xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  //
  //  graph containing nodes as vertices and footways as edges
  //
  graph<long long, double> G;
  updateGraph(Nodes, Footways, G);

  cout << endl;

  //
  // Stats
  //
  printStats(Nodes, Footways, Buildings, G);

  cout << endl;

  //
  // Menu
  //
  string userInput;
  cout << "Enter \"a\" for the standard application or "
       << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    application(Nodes, Footways, Buildings, G);
  }
  
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
