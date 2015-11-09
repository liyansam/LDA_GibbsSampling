#include <iostream>
#include "LDA.h"
using namespace std;

int main()
{
	int K, N, M, L;
	double alpha, beta;
	char parameterNode[20];
	char *setfile;
	char *filename;
	char *outputfile;
	outputfile = new char[20];
	setfile = new char[30];
	filename = new char[20];
	cout << ("Input the settings(filename):");
	cin >> setfile;
	ifstream fin(setfile);
	
	fin >> parameterNode;
	sscanf_s(parameterNode, "topic=%d", &K);
	fin >> parameterNode;
	sscanf_s(parameterNode, "word=%d", &N);
	fin >> parameterNode;
	sscanf_s(parameterNode, "doc=%d", &M);
	fin >> parameterNode;
	sscanf_s(parameterNode, "iteration=%d", &L);
	fin >> parameterNode;
	sscanf_s(parameterNode, "alpha=%lf", &alpha);
	fin >> parameterNode;
	sscanf_s(parameterNode, "beta=%lf", &beta);
	fin >> filename;
	filename = filename + 5;
	printf("parameters loaded in.\n");

	LDA myLDA(K, N, M, L, alpha, beta);
	myLDA.loadData(filename);
	myLDA.GibbsSampling();
	cout << "Input a filename for output:";
	cin >> outputfile;
	myLDA.printResult(outputfile);
	cout << "Done.";

	return 0;
}