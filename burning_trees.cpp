#include <omp.h>
#include <time.h>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <random>
#include <algorithm>


using namespace std;

// Function to generate random trees at a given seed with initial probability p
void generate_forest(std::vector<std::vector<int>>& forest, int seed, float p){
    int N = forest.size();
    // Initialise the random number generator outside the loop
    std::mt19937 generator(seed);
    // Initialise uniform distribution, U(0, 100) and convert p into integer P
    std::uniform_int_distribution<int> distribution(0, 100);
    int P = static_cast<int>(p * 100);
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                // generate random variable, X from the U(0, 100)
                int random_number = distribution(generator);
                // if X smaller than P, a tree is generated
                if (random_number < P){
                    forest[i][j] = 1;
                }
            }
        }
    }
}


// Function to print the current state of forest
void print_forest(std::vector<std::vector<int>> forest){
    for (const auto& row : forest) {
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;   
}


// Function to starts the burning process of forest fire 
void starts_to_burn(std::vector<std::vector<int>>& forest){
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < forest[0].size(); i++){
            // Burn only the top most part of the forest, count as step 0
            if (forest[0][i] == 1){
                forest[0][i] = 2;
            }
            else if (forest[0][i] == 0){
            }
        }
    }
}


// Function to check whether the input string is a false representation of wind directions
bool check_false_wind(string wind){
    bool check = true;
    // only accept left, right, up, down and none 
    std::vector<string> true_wind = {"left", "right", "up", "down", "none"};
    for (int i = 0; i < true_wind.size(); i++){
        // if the input string is the appropriate one return false
        if (wind == true_wind[i]){
            check = false;
        }
    }
    // else return true to show that the input string is false
    return check;
}


// Function to run the forest fire iteratively with consideration of wind direction
std::vector<std::vector<int>> parallel_burning_wind(std::vector<std::vector<int>> forest, string wind = "none"){
    int N = forest.size();
    // create a duplicate of the current forest state called new_forest
    std::vector<std::vector<int>> new_forest(N, std::vector<int>(N, 0));
    new_forest = forest;    
    // If the input string for wind direction is not correct, send error message
    if (check_false_wind(wind)){
        std::cerr << "Error: Function only accept 'left', 'right', 'up', 'down' or 'none' for wind argument. \n" << std::endl;
    }
    #pragma omp parallel
    {
        // Use static scheduling for multithreading
        #pragma omp for schedule(static) 
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                // If there is a burning tree with value 2 on current forest change new_forest to 3
                if (forest[i][j] == 2){
                    new_forest[i][j] = 3;
                    // If left wind let the fire to skip two tiles to the left
                    if (wind == "left"){
                        if (j > 1 && forest[i][j - 2] == 1){
                            new_forest[i][j - 2] = 2;
                        } 
                    }
                    // If right wind let the fire to skip two tiles to the right
                    if (wind == "right"){
                        if (j < N - 2 && forest[i][j + 2] == 1){
                            new_forest[i][j + 2] = 2;
                        }
                    }
                    // If up wind let the fire to skip two upper tiles
                    if (wind == "up"){
                        if (i > 1 && forest[i - 2][j] == 1){
                            new_forest[i - 2][j] = 2;
                        }  
                    }
                    // If down wind let the fire to skip two bottom tiles
                    if (wind == "down"){
                        if (i < N - 2 && forest[i + 2][j] == 1){
                            new_forest[i + 2][j] = 2;
                        }
                    }
                    // Burn the left neighbour trees
                    if (j != 0 && forest[i][j - 1] == 1){
                        new_forest[i][j - 1] = 2;
                    }
                    // Burn the right neighbour trees
                    if (j != N - 1 && forest[i][j + 1] == 1){
                        new_forest[i][j + 1] = 2;
                    }
                    // Burn the bottom neighbour trees
                    if (i != N - 1 && forest[i + 1][j] == 1){
                        new_forest[i + 1][j] = 2;
                    }
                    // Burn the upper neighbour trees
                    if (i != 0 && forest[i - 1][j] == 1){
                        new_forest[i - 1][j] = 2;
                    }
                }

            }
        }
    }
    // Return the updated version of burning forest
    return new_forest;
}


// Function to check whether two burning forest states are different
bool check_forest(std::vector<std::vector<int>> forest, std::vector<std::vector<int>> new_forest){
    // count enumerates the number of non-equal elements in the two random grids (forests)
    int count = 0;
    bool same = false;
    // If the size of both forests are different output error message
    if (forest.size() != new_forest.size() || forest[0].size() != new_forest[0].size()){
        std::cerr << "Error: The size of both forests are not equal." << std::endl;
    }
    else{
        #pragma omp parallel
        {
            // Using reduction to preserve the value of count over multiple threads
            #pragma omp for reduction(+:count)
            for (int i = 0; i < forest.size(); i++){
                for (int j = 0; j < forest.size(); j++){
                    // If the element of the forest is equal to the other forest enumerates count
                    if (forest[i][j] == new_forest[i][j]){
                        count++;
                    }
                }
            }
            // If count is not equal to the total elements in the forest grid return true
            // Means that there is a different in the values of elements between both forest grids
            if (count != forest.size() * forest.size()){
                same = true;
            }
        }
    }
    return same;
}


// Function to run forest fire in one random 'forest' grid, returns the running time and total steps of forest fire 
// 'show' argument is used to show all states of forest fire
std::vector<double> all_burn_wind(std::vector<std::vector<int>>& forest, string wind = "none", bool show = false){
    int N = forest.size();
    // Starts the timing 
    double start_time = omp_get_wtime();
    // Initiate the forest fire
    starts_to_burn(forest);
    // 'steps' to count the total steps
    double steps = 0;
    // Declare a 'new_forest' grid to represent the first step of forest fire
    std::vector<std::vector<int>> new_forest = parallel_burning_wind(forest, wind);
    // While the 'new_forest' do not equal to 'forest do
    while (check_forest(forest, new_forest)){
        #pragma omp parallel
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < N; i++){
                for (int j = 0; j < N; j++){
                    // Replace the elements in 'forest' with elements in 'new_forest'
                    forest[i][j] = new_forest[i][j];
                }
            } 
        }
        // Enumerate the forest fire steps
        steps++;
        // Continue to burn 'new_forest'
        new_forest = parallel_burning_wind(new_forest, wind);
        if (show){
            // If 'show' = true print out the current state of burning forest
            std::cout << "Forest state after burning for " << steps << " steps. \n" << std::endl; 
            print_forest(forest);
            std::cout << std::endl;
        }

    }
    // Calculate the total running time after the forest fire process
    double running_time = omp_get_wtime() - start_time;
    // Store the runtime and total steps into vector to be return 
    std::vector<double> step_time = {steps, running_time};
    if (show){
        // If 'show' = true print out the final state of forest, total steps and runtime
        std::cout << "Total steps to completely burnt : " << steps << std::endl;
        std::cout << "Total running time : " << running_time << std::endl;
        std::cout << "Forest state after burning completely. \n" << std::endl; 
        print_forest(forest);
    }
    return step_time;
}


// Function to create a random seeds for M random grids simulation 
int return_seed(){
    std::random_device rd;
    std:uniform_int_distribution<int> unif(INT_MIN, INT_MAX);

    return unif(rd);
}


// Function to estimate the probability of fire to reach the bottom of forest grid
double fire_bottom(std::vector<std::vector<int>>& forest){
    double reaches_bottom = 0;
    int n = forest.size();
    for (int j = 0; j < n; j++){
        if (forest[n - 1][j] == 3){
            reaches_bottom = 1;
            break;
        }
    }
    // If fire reaches bottom forest grid return 1
    return reaches_bottom;
}


// Function to run forest fire simulation on M randomly generated forest grids
// Takes input arguments number of random grids M, initial probability p, size of random grid N and direction of wind
// This function returns the average steps, runtime and probability of fire reaching bottom grids at initial probability, p
std::vector<double> M_model_wind(int M, float p, int N, string wind = "none", bool show = false){
    // Initialisation to store the total steps, total runtime and number of cases where fire reaches bottom grids
    double sum_step = 0;
    double sum_time = 0;
    double sum_bottom = 0;
    #pragma omp parallel
    {
        // Use reduction to preserve values sum_step, sum_time, sum_bottom
        #pragma omp for reduction (+ : sum_step, sum_time, sum_bottom) 
        for (int i = 0; i < M; i++){
            // Generate random seed
            int seed = return_seed();
            std::vector<std::vector<int>> forest(N, std::vector<int>(N, 0)); 
            // Randomly generate forest grid with the given random seed
            generate_forest(forest, seed, p);
            // Run the forest fire simulation until completes and collect the runtime and steps
            std::vector<double> info = {all_burn_wind(forest, wind, show)};
            // Sum up the runtime and steps from the simulation
            sum_step = sum_step + info[0];
            sum_time = sum_time + info[1];
            // Add one to sum_bottom if forest fire simulation reaches bottom
            sum_bottom = sum_bottom + fire_bottom(forest);
        }
    }
    // Return the initial probability p, average steps, average running time and the probability of fire reaches bottom
    std::vector<double> average_results = {p, sum_step / M, sum_time / M, sum_bottom / M};
    return average_results;
}


// Function to run forest fire simulation on M randomly generated forest grids of size N x N at all initial probability, p = [0,1)
// Takes arguments number of random grids M, size of random grid N and direction of wind
// Function returns vector containing the initial probability, average steps, runtime and probability of fire reaching bottom grids
// Addition 'print' argument to print out the collected data if true
std::vector<std::vector<double>> convergence_data_wind(int M, int N, string wind = "none", bool show = false, bool print = false){
    // Start the timing
    double start_time = omp_get_wtime();
    if (print){
        // if true, prints out collected data
        std::cout << "The average results for " << M << " runs and " << N << " grid size : " << std::endl;
        std::cout << std::endl;
        std::cout << "Probability|Steps|Run time|Probability reaches bottom" << std::endl;
    }
    // Initialise vector to store results
    std::vector<std::vector<double>> data; 
    for (int k = 0; k < 20; k++){
        // Run M random grids forest fire simulation for all initial probability p = {0, 0.05, 0.1, ...}
        double p = k * 0.05;
        std::vector<double> MN_results = M_model_wind(M, p, N, wind, show);
        // Store results into vector 'data'
        data.push_back({MN_results});
        if (print){
            std::cout << MN_results[0] << ", " << MN_results[1] << ", " << MN_results[2] << ", " << MN_results[3] <<  std::endl;
        }
    }
    // Calculate the runtime for the function to collect all the data
    double run_time = omp_get_wtime() - start_time;
    std::cout << std::endl;
    std::cout << "Data collection using " << omp_get_max_threads() << " threads." << std::endl;
    std::cout << "Running time to collect data for " << M << " runs and grid size " << N << 
    " with " << wind << " wind : " << run_time << " s " << std::endl; 
    return data;
}


// Function to write the collected data from simulation into a text file
void write_file(string filename, std::vector<std::vector<double>>& data){
    ofstream outfile;
    int nrow = data.size();
    int ncol = data[0].size();
    outfile.open(filename);
    for (int i = 0; i < nrow; i++){
        outfile << data[i][0];
        for (int j = 1; j < ncol; j++){
            outfile << "," << data[i][j];
        }
        outfile << endl;
    }
    outfile.close();
    std::cout << "Finish writing data into " << filename << " file." << std::endl;
}


// Function to write the current state of forest grid into a text file
void write_forest(string filename, std::vector<std::vector<int>>& data){
    ofstream outfile;
    int nrow = data.size();
    int ncol = data[0].size();
    outfile.open(filename);
    for (int i = 0; i < nrow; i++){
        outfile << data[i][0];
        for (int j = 1; j < ncol; j++){
            outfile << "," << data[i][j];
        }
        outfile << endl;
    }
    outfile.close();
}


int main(){

    //////////////////////////////////////////////////////////////////
    // 1) Creating a short simulation to understand the model rules //
    //////////////////////////////////////////////////////////////////

    int n = 20;       // Size of the forest grid (adjust as needed)
    float p = 0.5;    // Probability of a tree at each location
    int seed = 299;   // Seed for the random number generator
    omp_set_num_threads(15); // Setting the number of threads to 15
    // Generate a random forest grid of size 20 x 20, with initial probability 0.5
    std::vector<std::vector<int>> forest(n, std::vector<int>(n, 0));
    generate_forest(forest, seed, p);

    // Visualise the initial state of forest grid
    std::cout << "The Model Rules" << std::endl;
    std::cout << "Initial state of the forest : " << std::endl;
    print_forest(forest);
    // Store the forest grid into a text file to explain the model rules
    write_forest("forest", forest);

    // Create a duplicates of the forest grid to simulate a forest fire with wind direction
    std::vector<std::vector<int>> forest_down = forest;
    std::vector<std::vector<int>> forest_up = forest;
    std::vector<std::vector<int>> forest_left = forest;
    std::vector<std::vector<int>> forest_right = forest;

    // Testing the error message in check_forest function 
    std::vector<std::vector<int>> new_forest(n, std::vector<int>(n - 1, 0));
    check_forest(forest, new_forest);
    starts_to_burn(forest); 
    // Testing the error message in parallel_burning_wind function 
    parallel_burning_wind(forest, "is this wind");

    // Display the generated forest fire
    std::vector<double> results = all_burn_wind(forest, "none", false);
    std::cout << std::endl;


    /////////////////////////////////////////////////////////////////////////////////
    // 2) Creating forest fire simulations with presence of wind in all directions //
    /////////////////////////////////////////////////////////////////////////////////
    
    omp_set_num_threads(13); // Setting the number of threads to 13 
    // Run the forest fire with downward wind
    std::cout << "Forest fire with downward wind : " << std::endl;
    std::vector<double> results_down = all_burn_wind(forest_down, "down", true); // show each steps
    std::cout << std::endl;
    // Run the forest fire with leftward wind
    std::cout << "Forest fire with leftward wind : " << std::endl;
    std::vector<double> results_left = all_burn_wind(forest_left, "left", false); // change to true to show each steps
    std::cout << std::endl;
    // Run the forest fire with rightward wind
    std::cout << "Forest fire with rightward wind : " << std::endl;
    std::vector<double> results_right = all_burn_wind(forest_right, "right", false); // change to true to show each steps
    std::cout << std::endl;
    // Run the forest fire with upward wind
    std::cout << "Forest fire with upward wind : " << std::endl;
    std::vector<double> results_up = all_burn_wind(forest_up, "up", false); // change to true to show each steps
    std::cout << std::endl;

    // Store the final states of forest grids into a text file to explain the effect of wind direction on forest fire
    write_forest("forest_burn", forest);
    write_forest("forest_down", forest_down);
    write_forest("forest_left", forest_left);
    write_forest("forest_right", forest_right);
    write_forest("forest_up", forest_up);
    

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // 3) Collecting data to study the model convergence for random forest grid of size 100 x 100//
    ///////////////////////////////////////////////////////////////////////////////////////////////

    int N = 100;
    // Run forest fire on different number of random forest grids to check for model convergence
    // The data of average steps and runtime at every initial probabilities are exported into text files
    omp_set_num_threads(6); // Setting the number of threads to 6
    std::vector<std::vector<double>> converged_M1 = convergence_data_wind(1, N);
    write_file("converged_M1.csv", converged_M1);
    std::vector<std::vector<double>> converged_M5 = convergence_data_wind(5, N);
    write_file("converged_M5.csv", converged_M5);
    std::vector<std::vector<double>> converged_M10 = convergence_data_wind(10, N);
    write_file("converged_M10.csv", converged_M10);
    std::vector<std::vector<double>> converged_M25 = convergence_data_wind(25, N);
    write_file("converged_M25.csv", converged_M25);
    std::vector<std::vector<double>> converged_M50 = convergence_data_wind(50, N);
    write_file("converged_M50.csv", converged_M50);


    /////////////////////////////////////////////////////////////////////////////////////
    // 4) Collecting data to study the model convergence for 50 runs random forest grid//
    /////////////////////////////////////////////////////////////////////////////////////

    // Run forest fire on different size of random forest grids to check for model convergence
    // The data of average steps and runtime at every initial probabilities are exported into text files
    omp_set_num_threads(14); // Setting the number of threads to 14
    std::vector<std::vector<double>> converged_N5 = convergence_data_wind(50, 5);
    write_file("converged_N5.csv", converged_N5);
    std::vector<std::vector<double>> converged_N10 = convergence_data_wind(50, 10);
    write_file("converged_N10.csv", converged_N10);
    std::vector<std::vector<double>> converged_N25 = convergence_data_wind(50, 25);
    write_file("converged_N25.csv", converged_N25);
    std::vector<std::vector<double>> converged_N50 = convergence_data_wind(50, 50);
    write_file("converged_N50.csv", converged_N50);
    std::vector<std::vector<double>> converged_N100 = convergence_data_wind(50, 100);
    write_file("converged_N100.csv", converged_N100);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // 5) Collecting data to study the effect of wind direction onto average steps of forest fire//
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Run forest fire on 50 runs of random grids of size 100 x 100 with wind of different directions
    // The data of average steps and runtime at every initial probabilities are exported into text files
    omp_set_num_threads(10); // Setting the number of threads to 10
    std::vector<std::vector<double>> converged_M50_down = convergence_data_wind(50, 100, "down");
    write_file("converged_M50_down.csv", converged_M50_down);
    std::vector<std::vector<double>> converged_M50_left = convergence_data_wind(50, 100, "left");
    write_file("converged_M50_left.csv", converged_M50_left);
    std::vector<std::vector<double>> converged_M50_right = convergence_data_wind(50, 100, "right");
    write_file("converged_M50_right.csv", converged_M50_right);
    std::vector<std::vector<double>> converged_M50_up = convergence_data_wind(50, 100, "up");
    write_file("converged_M50_up.csv", converged_M50_up);


    ////////////////////////////
    // 6) Performance Analysis//
    ////////////////////////////

    // Need to study the effect of multithreading on different size of random grids
    // Fixed the initial probability to p = 0.6 and number of random grids to M = 50
    int P = 0.6;
    int M = 50;
    // Analyse the performance for random grid of size 50 x 50, 100 x 100 and 500 x 500
    int N1 = 50;
    int N2 = 100; 
    int N3 = 500;
    // Initialise vectors to store the runtime of forest fire simulations
    std::vector<std::vector<double>> walltime_1; 
    std::vector<std::vector<double>> walltime_2; 
    std::vector<std::vector<double>> walltime_3; 

    // For each sizes of random grids, run the forest fire simulations using different number of threads 1,...,16
    // For random forest grid of size 50 x 50
    for (int i = 1; i < 17; i++){ 
        double k = i;
        // Set the number of threads
        omp_set_num_threads(i);
            // Start timing
            double start = omp_get_wtime();
            // Run the forest fire simulations
            M_model_wind(M, P, N1);
            // Calculate the runtime
            double runtime = omp_get_wtime() - start;
            // Store the runtime to complete the simulations with k number of threads
            walltime_1.push_back({k, runtime / M});
    }

    // For random forest grid of size 100 x 100
    for (int i = 1; i < 17; i++){ 
        double k = i;
        // Set the number of threads
        omp_set_num_threads(i);
        // Start timing
            double start = omp_get_wtime();
            // Run the forest fire simulations
            M_model_wind(M, P, N2);
            // Calculate the runtime
            double runtime = omp_get_wtime() - start;
            // Store the runtime to complete the simulations with k number of threads
            walltime_2.push_back({k, runtime / M});
    }

    // For random forest grid of size 100 x 100
    for (int i = 1; i < 17; i++){ 
        double k = i;
        // Set the number of threads
        omp_set_num_threads(i);
        // Start timing
            double start = omp_get_wtime();
            // Run the forest fire simulations
            M_model_wind(M, P, N3);
            // Calculate the runtime
            double runtime = omp_get_wtime() - start;
            // Store the runtime to complete the simulations with k number of threads
            walltime_3.push_back({k, runtime / M});
    }
    std::cout << std::endl;

    // Export the runtime data to csv files to study the walltime and speedup of the multithreading
    write_file("walltime_1.csv", walltime_1);
    write_file("walltime_2.csv", walltime_2);
    write_file("walltime_3.csv", walltime_3);

    return EXIT_SUCCESS;
}