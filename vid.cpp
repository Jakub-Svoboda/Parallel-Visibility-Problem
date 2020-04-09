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
//#define NEUTRAL -1.57079633
#define NEUTRAL 0 //TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOooooooOooooooooooooooooooooooooooooooooooooooo

// ----------------------------- HEADERS -------------------------------
int32_t receiveNumber(int32_t processID, MPI_Status mpiStat);
std::vector<int> parseHeights(char* param);


// ----------------------------- CODE ----------------------------------

int32_t receiveNumber(int32_t processID, MPI_Status mpiStat){
	int32_t number;		
	MPI_Recv(&number, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &mpiStat); 
	return number;
}

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

double op(double a, double b){
	return a+b;
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


	switch(processID){
		case 0:
			leftMaxAngle = 3;
			rightMaxAngle = 1;
			break;
		case 1:
			leftMaxAngle = 7;
			rightMaxAngle = 0;
			break;
		case 2:
			leftMaxAngle = 4;
			rightMaxAngle = 1;
			break;
		case 3:
			leftMaxAngle = 6;
			rightMaxAngle = 3;
			break;
		default:
			leftMaxAngle = 1;
			rightMaxAngle = 1;	
	}



	double received = rightMaxAngle;

	// UPSWEEP:
	rightMaxAngle = op(leftMaxAngle, received);								
	for (int layer=0; layer < log2(numOfProcessors*2); layer++){
		int computingIds = pow(2,layer);		//1,2,4,8...
		int recIds = pow(2,layer+1);
		if ((processID+1) % computingIds == 0){		//if I am the layer that is working now:
			
			if((processID+1) % recIds == 0){			//I am receiving
				//cout << processID << " I have: " << leftMaxAngle << " and " << rightMaxAngle << endl;
				MPI_Recv(&received, 1, MPI_DOUBLE, processID-((int)pow(2,layer)), TAG, MPI_COMM_WORLD, &mpiStat);
				//cout << processID << " in layer "<< layer <<  " Receiving from " << processID-((int)pow(2,layer)) << " " << received << endl;
				rightMaxAngle = op(received, rightMaxAngle);

			}else if (processID+1 < n/2){
				//cout << processID << " in layer "<< layer <<  " SENDING TO " << processID+((int)pow(2,layer)) << " " << rightMaxAngle << endl;
				MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID+((int)pow(2,layer)), TAG, MPI_COMM_WORLD);		//send to parent
			}
		}
	}



	// CLEAR:
	if(processID+1 == n/2){		//If I am the last process:
		rightMaxAngle = NEUTRAL;//TODDDDDDDDDDDDDDDDDDDDDOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
	}



	// DOWNSWEEP:
	for (int layer=log2(numOfProcessors); layer > 0; layer--){
		int computingIds = pow(2,layer);		//....8,4,2,1
		int sendIds = pow(2,layer-1);
		if ((processID+1) % sendIds == 0 && ((processID+1) % computingIds != 0)){		//If I am the left child process below current level		
			MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID+sendIds, TAG, MPI_COMM_WORLD);
			//cout << processID << " (child) in layer "<< layer <<  " SENDING TO " << processID+sendIds << " " << rightMaxAngle << endl;
			MPI_Recv(&rightMaxAngle, 1, MPI_DOUBLE, processID+sendIds, TAG, MPI_COMM_WORLD, &mpiStat);	
			//cout << processID << " (child) in layer "<< layer <<  " Receiving from " << processID+sendIds << " " << rightMaxAngle << endl;
		}
		if ((processID+1) % computingIds == 0){		//if I am the layer that is working now (receiving first round):
			double received;
			MPI_Recv(&received, 1, MPI_DOUBLE, processID-sendIds, TAG, MPI_COMM_WORLD, &mpiStat);
			//cout << processID << " (parent) in layer "<< layer <<  " Receiving from " << processID-sendIds << " " << received << endl;
			MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, processID-sendIds, TAG, MPI_COMM_WORLD);		//send back left
			//cout << processID << " (parent) in layer "<< layer <<  " SENDING TO " << processID-sendIds << " " << rightMaxAngle << endl;
			rightMaxAngle = op(received, rightMaxAngle);
		}
	}
	double tmp = leftMaxAngle;
	leftMaxAngle = rightMaxAngle;
	rightMaxAngle = op(leftMaxAngle, tmp);

	

	//cout << processID << "  " <<  leftMaxAngle << " " << rightMaxAngle << endl;
	
	MPI_Send(&leftMaxAngle, 1, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
	MPI_Send(&rightMaxAngle, 1, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);


	std::vector<double> results;
	if(processID == 0){
		for(int i=0; i<numOfProcessors;i++){
			results.push_back(0);
			results.push_back(0);
			MPI_Recv(&results[2*i], 1, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, &mpiStat);
			MPI_Recv(&results[2*i+1], 1, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, &mpiStat);
		}
		for(int i = 0; i<trueSize; i++) 
			cout << results[i] << " " ;
			
		cout << endl;
	}



	MPI_Finalize(); 
	return 0;
}