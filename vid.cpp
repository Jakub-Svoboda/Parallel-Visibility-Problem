/*
	Parallel visibility problem 
	Author: Jakub Svoboda
	Date 5.4.2020
	Login: xsvobo0z
	Email: xsvobo0z@stud.fit.vutbr.cz
*/

#include <math.h> 
#include <stdio.h>
#include <string>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <unistd.h>

using namespace std;

#define TAG 0
#define NEUTRAL -1.57079633


// ----------------------------- HEADERS -------------------------------
std::vector<int> parseHeights(char* param);
int getTrueSize(char* param);
std::vector<int> parseHeights(char* param, int numOfProcessors);
double op(double a, double b);
char getLeftAnswer(int32_t processID, double leftAngle, double leftMaxAngle);
char getRightAnswer(int32_t processID, double rightAngle, double rightMaxAngle);
void outputResults(int numOfProcessors, int trueSize,MPI_Status mpiStat);


// ----------------------------- CODE ----------------------------------
/**
 * @brief Get the True Size of the input - the number of heights to work with
 * 
 * @param param command line input
 * @return int Number of heights passed in parameter.
 */
int getTrueSize(char* param){
	std::vector<int> vect;
    std::stringstream ss(param);
	for (int i; ss >> i;) {
		vect.push_back(i);    
		if (ss.peek() == ',')
			ss.ignore();
	}
	return vect.size();
}

/**
 * @brief Parses the command line parameter into a vector of integers representing heights
 * 
 * @param param Command line parameter
 * @param numOfProcessors Number of processes alocated
 * @return std::vector<int> Vector of heights
 */
std::vector<int> parseHeights(char* param, int numOfProcessors){
	std::vector<int> vect;
    std::stringstream ss(param);

	for (int i; ss >> i;) {
		vect.push_back(i);    
		if (ss.peek() == ',')
			ss.ignore();
	}
	
	if(vect.size() < numOfProcessors*2){
		vect.resize(numOfProcessors*2);
		//cout << vect.size() << "  " << numOfProcessors << endl;
	}
	return vect;
}

/**
 * @brief Operation for prescan. In this case - MAX
 * 
 * @param a Left operand
 * @param b Right Operand
 * @return double result
 */
double op(double a, double b){
	return max(a,b);
}

/**
 * @brief Convert the left angle and max angle to visibility.
 * 
 * @param processID ID of calling process
 * @param leftAngle 
 * @param leftMaxAngle 
 * @return char _, v, u
 */
char getLeftAnswer(int32_t processID, double leftAngle, double leftMaxAngle){
	if (processID == 0){
		return '_';
	}else if (leftAngle > leftMaxAngle){
		return 'v';
	}else{
		return 'u';
	}
}

/**
 * @brief Convert the right angle and max angle to visibility.
 * 
 * @param processID ID of calling process
 * @param rightAngle 
 * @param rightMaxAngle 
 * @return char _, v, u
 */
char getRightAnswer(int32_t processID, double rightAngle, double rightMaxAngle){
	if (rightAngle > rightMaxAngle){
		return 'v';
	}else{
		return 'u';
	}
}

/**
 * @brief Process zero prints the vector in the correct format and order.
 * 
 * @param numOfProcessors Number of processes
 * @param trueSize True size of the passed input
 * @param mpiStat 
 */
void outputResults(int numOfProcessors, int trueSize,MPI_Status mpiStat){
	std::vector<char> results;
	for(int i=0; i<numOfProcessors;i++){
		results.push_back(0);
		results.push_back(0);
		MPI_Recv(&results[2*i], 1, MPI_CHAR, i, TAG, MPI_COMM_WORLD, &mpiStat);
		MPI_Recv(&results[2*i+1], 1, MPI_CHAR, i, TAG, MPI_COMM_WORLD, &mpiStat);
	}
	for(int i = 0; i<trueSize; i++) 
		cout << results[i] << " " ;
		
	cout << endl;
}


int main(int argc, char* argv[]){
	int32_t numOfProcessors = 0;							//zero processors by default
	int32_t processID = 0;									//ID of this process

	MPI_Init(&argc, &argv);									//initialize MPI
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcessors); 		//get the number of availible processes
	MPI_Comm_rank(MPI_COMM_WORLD ,&processID); 				//get the ID of this process

	MPI_Status mpiStat;

	if (numOfProcessors == 0){
		cerr << "The number of processors cannot be zero." << endl;
		exit(1);
	}

	std::vector<int> heights;
	int leftHeight, rightHeight, observerHeight, trueSize, n;

	//Parse heights from argument
	if (processID == 0){
		heights = parseHeights(argv[1], numOfProcessors);
		trueSize = getTrueSize(argv[1]);
		n = heights.size();
	}

	//Broadcast TrueHeight
	MPI_Bcast(&trueSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//Broadcast n (padded)
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Send left heights
	if (processID == 0) {	
		for (int i = 0; i < numOfProcessors; i++){
			MPI_Send(&heights[i*2], 1, MPI_INT, i, TAG, MPI_COMM_WORLD); 
		}
		observerHeight = heights[0];
	}
	MPI_Recv(&leftHeight, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &mpiStat);

	//Send right heights
	if (processID == 0) {	
		for (int i = 0; i < numOfProcessors; i++){
			MPI_Send(&heights[i*2+1], 1, MPI_INT, i, TAG, MPI_COMM_WORLD); 
		}
	}
	MPI_Recv(&rightHeight, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &mpiStat);
	

	//Send observerHeight
	if (processID == 0) {	
		for (int i = 0; i < numOfProcessors; i++){
			MPI_Send(&observerHeight, 1, MPI_INT, i, TAG, MPI_COMM_WORLD); 
		}
	}	
	MPI_Recv(&observerHeight, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &mpiStat);


	//Calculate angle
	double leftAngle, rightAngle;
	leftAngle = atan((leftHeight - observerHeight) / (double)(processID*2));
	rightAngle = atan((rightHeight - observerHeight) / (double)(processID*2+1));

	if (processID == 0) leftAngle = NEUTRAL;		//fix div 0

	double leftMaxAngle = leftAngle;
	double rightMaxAngle = rightAngle;

	double received = rightMaxAngle;


	// ------- UPSWEEP -------
	rightMaxAngle = op(leftMaxAngle, received);								
	for (int layer=0; layer < log2(numOfProcessors*2); layer++){
		int computingIds = pow(2,layer);		//1,2,4,8...
		int recIds = pow(2,layer+1);
		if ((processID+1) % computingIds == 0){		//if I am the layer that is working now:
			if((processID+1) % recIds == 0){			//I am receiving
				MPI_Recv(&received, 1, MPI_DOUBLE, processID-((int)pow(2,layer)), TAG, MPI_COMM_WORLD, &mpiStat);
				rightMaxAngle = op(received, rightMaxAngle);
			}else if (processID+1 < n/2){
				MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID+((int)pow(2,layer)), TAG, MPI_COMM_WORLD);		//send to parent
			}
		}
	}


	// ------- CLEAR -------
	if(processID+1 == n/2){		//If I am the last process:
		rightMaxAngle = NEUTRAL;
	}


	// ------- DOWNSWEEP -------
	for (int layer=log2(numOfProcessors); layer > 0; layer--){
		int computingIds = pow(2,layer);		//....8,4,2,1
		int sendIds = pow(2,layer-1);
		if ((processID+1) % sendIds == 0 && ((processID+1) % computingIds != 0)){		//If I am the left child process below current level		
			MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID+sendIds, TAG, MPI_COMM_WORLD);
			MPI_Recv(&rightMaxAngle, 1, MPI_DOUBLE, processID+sendIds, TAG, MPI_COMM_WORLD, &mpiStat);	
		}
		if ((processID+1) % computingIds == 0){		//if I am the layer that is working now (receiving first round):
			double received;
			MPI_Recv(&received, 1, MPI_DOUBLE, processID-sendIds, TAG, MPI_COMM_WORLD, &mpiStat);
			MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID-sendIds, TAG, MPI_COMM_WORLD);		//send back left
			rightMaxAngle = op(received, rightMaxAngle);
		}
	}
	double tmp = leftMaxAngle;
	leftMaxAngle = rightMaxAngle;
	rightMaxAngle = op(leftMaxAngle, tmp);

	
	// ------- OUTPUT -------
	char leftAnswer = getLeftAnswer(processID, leftAngle, leftMaxAngle);
	char rightAnswer = getRightAnswer(processID, rightAngle, rightMaxAngle);

	// every process sends the visibility char to process 0
	MPI_Send(&leftAnswer, 1, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
	MPI_Send(&rightAnswer, 1, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);

	if(processID == 0){		//Processor 0 will receive results and print them
		outputResults(numOfProcessors, trueSize, mpiStat);	
	}

	MPI_Finalize(); 
	return 0;
}