#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <omp.h> 
#include <chrono>

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
				sin >> x >> y;

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

	double path_cost(vector<int> cities, vector<vector<double>> distance_matrix) {
		double pathWeight = 0;

		for (int i = 0; i < cities.size() - 1; i++)
		{
			//int r = cities[i] - 1;
			//int c = cities[i + 1] - 1;
			pathWeight += distance_matrix[cities[i] - 1][cities[i + 1] - 1];
		}
		//path from last to first
		//int r = cities[cities.size() - 1] - 1;
		//int c = cities[0] - 1;
		pathWeight += distance_matrix[cities[cities.size() - 1] - 1][cities[0] - 1];
		return pathWeight;
	}

public:

	int read_tsp_file(const char* fname)
	{
		ifstream file(fname);
		vector<double> xs, ys;
		vector<vector<double>> temp;
		vector<vector<double>> distance_matrix;

		// Using time point and system_clock
		std::chrono::time_point<std::chrono::system_clock> start, end;

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
		start = std::chrono::system_clock::now();
		//TODO: calculate shortest possible route
		vector<int> cities;

		/// pushing the rest num_nodes-1 cities into a bundle
		for (int i = 1; i < distance_matrix.size() + 1; i++)
			cities.push_back(i);

		int min_dist = INT_MAX;
		vector<int> shortest_path{};

		/*parallel here */
		/// generating permutations and tracking the minimum cost
#pragma omp parallel
		{
			double dist = 0;
#pragma omp for
			for (int i = 0; i < cities.size(); ++i)
			{
				auto localPath = cities;

				rotate(localPath.begin(), localPath.begin() + i, localPath.begin() + i + 1);
				do
				{
					dist = path_cost(localPath, distance_matrix);

#pragma omp critical

					{
						if (dist < min_dist)
						{
							shortest_path = localPath;
							min_dist = dist;
						}
					}

				} while (next_permutation(localPath.begin() + 1, localPath.end()));

			}


		}

		end = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds = end - start;

		cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
		cout << "path: ";
		for (int i : shortest_path)
			cout << i << " --> ";
		cout << '\n';

		return min_dist;

	}
};

/*
@brief
Optimization of travelling salesman problem using Branch and bound algorithm
*/
class TSP_BranchAndBound {

	vector<vector<double>> readCoordinationDataFromFile(const char* fname) {
		ifstream file(fname);
		vector<double> xs, ys;
		vector<vector<double>> coordinationData;

		// Using time point and system_clock
		std::chrono::time_point<std::chrono::system_clock> start, end;

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

				coordinationData.push_back({ x,y });
			}

			file.close();
		}
		else
		{
			cout << fname << " file not open" << endl;
			return {};
		}

		return coordinationData;
	}

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

	// This function stores transpose
	// of A[][] in B[][]
	/*
	void transpose(int A[][N], int B[][N])
	{
		int i, j;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				B[i][j] = A[j][i];
	}
	*/

	vector<vector<double>> reduce_matrix(vector<vector<double>> m) {

		vector<double> min_elements = {};
		//rows
		for (int i = 0; i < m.size(); i++)
		{
			vector<double> row = m[i];
			auto min_element = *std::min_element(begin(row), end(row));
			min_elements.push_back(min_element);
			for (auto& val : row)
				val -= min_element;
			m[i] = row;
		}
		cout << "test\n";
		//cols 
		
		for (int cols = 0; cols < m.size(); cols++)
		{
			for (int rows = 0; rows < m.size(); rows++)
			{

			}
		}

		return {};
	}

public:
	void solve(const char* fileName) {
		
		//vector<vector<double>> cords = readCoordinationDataFromFile(fileName);
		//calculate distance matrix
		//vector<vector<double>> distance_matrix = ComputeEuclideanDistanceMatrix(cords, cords.size());
		vector<vector<double>> distance_matrix = {{INT_MAX, 5, 40, 11},
												  {5, INT_MAX, 9, 6 },
												  {40, 9, INT_MAX, 8 },
												  {11, 6, 8, INT_MAX }};
		for (int i = 0; i < distance_matrix.size(); i++) {
			print_distances(distance_matrix[i], distance_matrix.size());
			cout << '\n';
		}


		/*
		1. matrix reduction function
		2. 
		*/
		reduce_matrix(distance_matrix);
	}

};

int main()
{
	TSP_bruteforce_int sol1;
	TSP_bruteforce_double sol2;
	TSP_BranchAndBound sol3;

	//cout << "_double---------------------------------------------------------------\n";
	//cout << "shortest path: " << sol2.read_tsp_file("ulysses222.tsp.txt") << '\n';
	//cout << "_double---------------------------------------------------------------\n";
	cout << endl;
	//cout << "_int---------------------------------------------------------------\n";
	//cout << "shortest path: " << sol1.read_tsp_file("ulysses22.tsp.txt") << '\n';
	//cout << "_int---------------------------------------------------------------\n";
	cout << "Branch and Bound---------------------------------------------------------------\n";
	sol3.solve("ulysses222.tsp.txt");
	cout << "Branch and Bound---------------------------------------------------------------\n";
	cin.get();
	 
	return 0;

	//chovani next permutation:
	//dycky se na prvni misto dosadi neco jineho

	/*1,2,3,4
	1243
	1324
	1342
	1423
	1432
	*/
	//vsechny permutace kde na zacatku je 1 bude zpracovavat 1 vlakno
	//vsechny permutace kde na zacatku je 2 bude zpracovavat 2 vlakno
	/*
	1234
	2134
	3124
	4123
	- kazde vlakno muze permutovat 1 vector
	- pouzit rotate tak, abych dycky na zacatku mel to cislo ktere chci
	*/
	/*
	1234
	2134
	3124
	4123
	rotate(path.begin, path.begin+i, path.begin+i+1)
	-- 14
	*/
}
