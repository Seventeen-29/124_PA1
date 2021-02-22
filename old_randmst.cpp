/**
*Programming Assignemnt 1, Source Code
*Chris Zhu, Rakesh Nori
# Note: dimension 0 is the same thing as dimension 1
**/
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <time.h>
#include <chrono>

int numDimensions;
int numTrials;
double maxi;
using namespace std;
vector<bool> vis;
vector<double> currPoints;
vector<vector<double> > graphVertices;

double dist(vector<double> points1, vector<double> points2){
  double sum = 0;
  for (int i = 0; i < points1.size(); ++i)
  {
  	double curr = (points1[i] - points2[i]);
  	sum += curr * curr;
  }
  //cout << "SUM: " << sqrt(sum) << endl;
  return sqrt(sum);
}

vector<double> generatePoints(int numDimensions){
	vector<double> points;
	for (int i = 0; i < numDimensions; i++)
	{
		points.push_back((double) rand()/RAND_MAX);
		//cout << points[i] << endl;
	}
	return points;
}

int findMin(vector<double> vec, vector<bool> visited){
	int mindex = -1;
	double min = INT_MAX;
	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] < min && visited[i] == 0)
		{
			min = vec[i];
			mindex = i;
		}
	}
	return mindex;
}

double prims(int n, int numDimensions){
	double totalSum = 0;
	vector<bool> vis(n, 0);
	vector<double> distances(n, 10000000);
	double totalEdgeWeight = 0;
	for (int i = 0; i < n; i++)
	{
		graphVertices.push_back(generatePoints(numDimensions));  // generate graph
	}
	distances[0] = 0;
	for (int i = 0; i < n; i++)
	{
		int currIndex = findMin(distances, vis);
		totalEdgeWeight += distances[currIndex];
		distances[currIndex] = INT_MAX;
		vis[currIndex] = 1;
		currPoints = graphVertices[currIndex];
		for (int j = 0; j < n; j++)
		{
			if (!vis[j])
			{
				double distancia;
				if (numDimensions == 0)
					distancia = ((double) rand()/RAND_MAX);
				else
					distancia = dist(graphVertices[j], currPoints);
				if (distancia < distances[j])
				{
					distances[j] = distancia;
				}
			}
		}
	}
	return totalEdgeWeight;
}



int main(int argc, char *argv[]){

	//auto start = chrono::high_resolution_clock::now(); 

	srand( (unsigned)time( NULL ) );
	int flag = stoi(argv[1]);
	int n = stoi(argv[2]);
	numTrials = stoi(argv[3]);
	numDimensions = stoi(argv[4]);

	double results = 0.0;
	for (int i = 0; i < numTrials; i++)
	{
		results +=  prims(n, numDimensions);
	}
	results = results / numTrials;

	//auto stop = chrono::high_resolution_clock::now(); 
	//auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start); 

	//cout << duration.count() << " ms" << endl; 

	cout << "N = " << n << ". Average MST total weight for our trials: " << results << endl;

}

// n ^(1 - 1/d) --> RUNTIME