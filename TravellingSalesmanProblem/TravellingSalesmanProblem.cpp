#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <omp.h> 


using namespace std;

class TSP_bruteforce_int {

	//print single row from a (distance) matrix
	void print_distances(vector<int>& dist, unsigned int n)
	{
		for (unsigned int i = 0; i < dist.size(); i++)
		{
			cout << dist[i] << ", ";
		}
		cout << endl;
	}

	//euclidean distance between 2 points
	float calculate_distance(int x1, int y1, int x2, int y2)
	{
		// Calculating distance
		int res = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
		return res;
	}

	// @brief Generate distance matrix.
	vector<vector<int>> ComputeEuclideanDistanceMatrix(const vector<vector<int>>& locations, int n)
	{
		vector<vector<int>> distances = vector<vector<int>>(n, vector<int>(n, int{ 0 }));
		for (int startingCity = 0; startingCity < n; startingCity++)
		{
			for (int endingCity = 0; endingCity < n; endingCity++)
			{
				if (startingCity != endingCity) {
					distances[startingCity][endingCity] = calculate_distance(locations[startingCity][0], locations[startingCity][1], locations[endingCity][0], locations[endingCity][1]);
				}
			}
		}
		return distances;
	}

public:

	int read_tsp_file(const char* fname)
	{
		ifstream file(fname);
		vector<int> xs, ys;
		vector<vector<int>> temp;
		vector<vector<int>> distance_matrix;

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
				sin >> id;
				sin	>> x >> y;

				temp.push_back({ static_cast<int>(x), static_cast<int>(y) });
			}

			unsigned int n = temp.size();

			//calculate distance matrix
			//distance_matrix = ComputeEuclideanDistanceMatrix(temp, n);

			// matrix representation of graph
			distance_matrix = { { 0, 10, 15, 20 },
							   { 10, 0, 35, 25 },
							   { 15, 35, 0, 30 },
							   { 20, 25, 30, 0 } };

			for (int i = 0; i < distance_matrix.size(); i++) {
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

		/// pushing the rest num_nodes-1 cities into a bundle
		for (int i = 0; i < distance_matrix.size(); i++)
		{
			if (i != source)
			{
				nodes.push_back(i);
			}
		}
		int n = nodes.size();
		int shortest_path = INT_MAX;

//#pragma omp parallel num_threads(4)

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

};

class TSP_bruteforce_double {

	void print_distances(vector<double>& dist, unsigned int n)
	{
		for (unsigned int i = 0; i < dist.size(); i++)
		{
			cout << dist[i] << ", ";
		}
		cout << endl;
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
		vector<vector<double>> distances = vector<vector<double>>(n, vector<double>(n, double{ 0 }));
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

public:

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

			//calculate distance matrix
			distance_matrix = ComputeEuclideanDistanceMatrix(temp, n);

			for (int i = 0; i < distance_matrix.size(); i++) {
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

		/// pushing the rest num_nodes-1 cities into a bundle
		for (int i = 0; i < distance_matrix.size(); i++)
		{
			if (i != source)
			{
				nodes.push_back(i);
			}
		}
		int n = nodes.size();
		int shortest_path = INT_MAX;

		/// generating permutations and tracking the minimum cost
		/**/while (next_permutation(nodes.begin(), nodes.end())) {
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
};

int main()
{
	TSP_bruteforce_int sol1;
	TSP_bruteforce_double sol2;

	cout << "_double---------------------------------------------------------------\n";
	cout << "shortest path: " << sol2.read_tsp_file("ulysses222.tsp.txt") << '\n';
	cout << "_double---------------------------------------------------------------\n";
	cout << endl;
	cout << "_int---------------------------------------------------------------\n";
	cout << "shortest path: " << sol1.read_tsp_file("ulysses22.tsp.txt") << '\n';
	cout << "_int---------------------------------------------------------------\n";
	cin.get();

	return 0;
}
