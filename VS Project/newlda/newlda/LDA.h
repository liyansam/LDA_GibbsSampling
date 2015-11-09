#ifndef LDA_H
#define LDA_H
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
using namespace std;

class LDA
{
private:
	int K, N, M, L;
	double alpha, beta;
	int **nDocWord;
	int *cntDocWord;
	int **nDocTopic, *sumDocTopic, **nTopicTerm, *sumTopicTerm;
	double **phiTopicTerm, **thetaDocTopic;
	int **kDocWord, **tDocWord;
	double *newDist;
public:
	LDA(int myk, int myn, int mym, int myl, double myalpha, double mybeta);
	int sample(double *distribution, int size);
	void loadData(char* filename);
	void calcNewDist(int t, int m);
	void GibbsSampling();
	void printResult(char* outputfile);
};

#endif