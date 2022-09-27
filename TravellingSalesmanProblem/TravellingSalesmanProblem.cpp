#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>


using namespace std;

class TSP_bruteforce1 {

};

vector<vector<double>> ComputeEuclideanDistanceMatrix(const vector<vector<double>>& locations, int n);
int read_tsp_file(const char* fname);

vector<double> distances;

void print_distances(vector<double>& dist, unsigned int n)
{
	for (unsigned int i = 0; i < dist.size(); i++)
	{
		cout << dist[i] << ", ";	
	}
	cout << endl;
}

int read_tsp_file(const char* fname)
{
	ifstream file(fname);
	vector<double> xs, ys;
	vector<vector<double>> temp;
	vector<vector<double>> distance_matrix;

	if (file.is_open())
	{
		string line;

		getline(file, line);
		getline(file, line);
		getline(file, line);
		getline(file, line);
		getline(file, line);
		getline(file, line);
		getline(file, line);

		while (std::getline(file, line)) {
			if (line[0] == 'E')
				break;

			stringstream sin(line);
			int id;
			double x, y;
			sin >> id >> x >> y;

			temp.push_back({ x,y });
		}

		unsigned int n = temp.size();

		distances.resize(n * n);

		//calculate distance matrix
		distance_matrix = ComputeEuclideanDistanceMatrix(temp, n);
		
		for (int i = 0; i < n; i++) {
			print_distances(distance_matrix[i], n);
			cout << '\n';
		}

		file.close();
	}
	else
	{
		cout << fname << " file not open" << endl;
		return 0;
	}

	//TODO: calculate shortest possible route

	//choosing 1st city
	int source = 0;
	vector<int> nodes;
	vector<vector<int>> path;

	for (int i = 0; i < nodes.size(); i++)
		cout << nodes[i] << 't';
	/// pushing the rest num_nodes-1 cities into a bundle
	for (int i = 0; i < 22; i++)
	{
		if (i != source)
		{
			nodes.push_back(i);
		}
	}
	int n = nodes.size();
	int shortest_path = INT_MAX;
	
	/// generating permutations and tracking the minimum cost
	while (next_permutation(nodes.begin(), nodes.end())) {
		//path.clear();
		int path_weight = 0;
		int j = source;

		for (int i = 0; i < n; i++)
		{
			path_weight += distance_matrix[j][nodes[i]];
			j = nodes[i];
			//path.push_back({ j,nodes[i] });
		}

		path_weight += distance_matrix[j][source];
		
		shortest_path = min(shortest_path, path_weight);

	}
	
	return shortest_path;
	
}

//euklidovska vzdalenost mezi dvema bodama
float distance(double x1, double y1, double x2, double y2)
{
    // Calculating distance
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

// @brief Generate distance matrix.
vector<vector<double>> ComputeEuclideanDistanceMatrix(const vector<vector<double>>& locations, int n)
{
	vector<vector<double>> distances = vector<vector<double>>(n, vector<double>(n, double{0}));
	for (int startingCity = 0; startingCity < n; startingCity++) 
	{
		for (int endingCity = 0; endingCity < n; endingCity++) 
		{
			if (startingCity != endingCity) {
				distances[startingCity][endingCity] = static_cast<double>(distance(locations[startingCity][0], locations[startingCity][1], locations[endingCity][0], locations[endingCity][1]));
			}
		}
	}
	return distances;
}

int main()
{
	cout << "shortest path: " << read_tsp_file("ulysses22.tsp.txt") << '\n';
	cin.get();

	return 0;
}
