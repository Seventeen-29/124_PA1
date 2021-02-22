/**
*Programming Assignemnt 1, Source Code
*Chris Zhu, Rakesh Nori
**/
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <time.h>
#include <chrono>

using namespace std;
const double scalar = 3.0;

double dist(vector<double> pt1, vector<double> pt2){
  double sum = 0;
  for (int i = 0; i < pt1.size(); ++i){
  	double delta = (pt1[i] - pt2[i]);
  	sum += delta * delta;
  }
  return sqrt(sum);
}

vector<double> randPoint(int numDimensions){
	vector<double> pt;
	for (int i = 0; i < numDimensions; i++){
		pt.push_back((double) rand()/RAND_MAX);
	}
	return pt;
}

int getBoxCoord(vector<double> pt, double weight_threshold, int numDimensions){
	int grid_size = ceil (1 / weight_threshold);
	int index_accumulator = 0;
	for(int j = 0; j < numDimensions; j++){ //iterate across currPoint to get box index
		index_accumulator += (int) ((floor (pt[j] / weight_threshold)) * pow(grid_size, j));
	}
	return index_accumulator;
}

vector<int> getNearbyBoxes(int currBox, int gridSize, int numDimensions){
	int max = ((int) pow(gridSize, numDimensions));
	vector<int> temp;
	temp.push_back(0);
	for(int i = 0; i < numDimensions; i++){
		vector<int> curr;
		for(int elt : temp){
			curr.push_back(elt);
			curr.push_back(elt + (int) pow(gridSize, i));
			curr.push_back(elt - (int) pow(gridSize, i));
		}
		temp = curr;
	}
	vector<int> output;
	for(int elt : temp){
		if(elt + currBox < max && elt + currBox >= 0){
			output.push_back(elt + currBox);
		}
	}
	return output;
}

double prims(int numVertices, int numDimensions){
	if (numDimensions <= 0){
		throw invalid_argument("bad call prims_nonzero on dim <= 0");
	}
	vector<vector<double>> graphVertices;
	vector<bool> visited(numVertices, false);
	vector<double> distances(numVertices, numDimensions * 1.0);
	distances[0] = 0.0; //pick starting vertex to be index 0
	double totalEdgeWeight = 0;
	double weight_threshold = sqrt(numDimensions); //worst case distance bound, unused for n = 0, 2, 3, 4

	if(numDimensions == 2){
		weight_threshold = scalar * 0.7 * pow(numVertices, -0.5);
	}
	else if(numDimensions == 3){
		weight_threshold = scalar * pow(numVertices, -0.37);
	}
	else if(numDimensions == 4){
		weight_threshold = scalar * pow(numVertices, -0.29) + pow(numVertices, -1.7);
	}

	int grid_size = ceil (1 / weight_threshold);
	vector<vector<int>> boxes(pow(grid_size, numDimensions));

	for (int i = 0; i < numVertices; ++i){
		vector<double> currPoint = randPoint(numDimensions);
		graphVertices.push_back(currPoint); //generate the graph
		int currIndex = getBoxCoord(currPoint, weight_threshold, numDimensions);
		boxes[currIndex].push_back(i);
	}

	for (int i = 0; i < numVertices; ++i){
		int minIndex = -1;
		double minDist = INT_MAX;
		//find minimal distance vertex that is not yet visited, get index
		for(int k = 0; k < numVertices; ++k){
			if (!visited[k] && (distances[k] < minDist)){
				minIndex = k;
				minDist = distances[k];
			}
		}

		totalEdgeWeight += distances[minIndex];
		visited[minIndex] = true;

		for (int nearbyBox : getNearbyBoxes(getBoxCoord(graphVertices[minIndex],
		 weight_threshold, numDimensions), grid_size, numDimensions)){
			for (int nearbyVertexIndex : boxes[nearbyBox]){
				if (!visited[nearbyVertexIndex]){
					double currDist = dist(graphVertices[nearbyVertexIndex], graphVertices[minIndex]);
					if (currDist < distances[nearbyVertexIndex]){
						distances[nearbyVertexIndex] = currDist;
					}
				}
			}
		}
	if (i % 150 == 0){
			printf ("%2.2f%%\n", 100 * ((double) i) / numVertices);
		}
	}

	return totalEdgeWeight;
}

int main(int argc, char *argv[]){	
	srand((unsigned)time(NULL));
	int flag = stoi(argv[1]);
	int n = stoi(argv[2]);
	int numTrials = stoi(argv[3]);
	int numDimensions = stoi(argv[4]);

	auto start = chrono::high_resolution_clock::now(); 

	double results = 0.0;
	for (int i = 0; i < numTrials; i++)
	{
		results +=  prims(n, numDimensions);
	}
	results = results / numTrials;

	auto stop = chrono::high_resolution_clock::now(); 
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 

	cout << "N = " << n << ". Dim = " << numDimensions << ". AVG MST weight: " << results << " [" <<
  duration.count() << " ms]" << endl;
}

// n ^(1 - 1/d) --> RUNTIME