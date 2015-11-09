#include "LDA.h"

LDA::LDA(int myk, int myn, int mym, int myl, double myalpha, double mybeta): K(myk), N(myn), M(mym), L(myl), alpha(myalpha), beta(mybeta)
{
	//initialisation - allocate space and zero clearing
	sumDocTopic = new int[M];
	sumTopicTerm = new int[K];
	cntDocWord = new int[M];
	newDist = new double[K];

	nDocWord = new int*[M];
	for (int i = 0; i < M; i++)
		nDocWord[i] = new int[N];

	nDocTopic = new int*[M];
	for(int i=0; i<M; i++)
		nDocTopic[i] = new int[K];

	nTopicTerm = new int*[K];
	for(int i=0; i<K; i++)
		nTopicTerm[i] = new int[N];

	kDocWord = new int*[M];
	for(int i=0; i<M; i++)
		kDocWord[i] = new int[N];

	tDocWord = new int*[M];
	for(int i=0; i<M; i++)
		tDocWord[i] = new int[N];

	thetaDocTopic = new double*[M];
	for (int i = 0; i < M; i++)
		thetaDocTopic[i] = new double[K];

	phiTopicTerm = new double*[K];
	for (int i = 0; i < K; i++)
		phiTopicTerm[i] = new double[N];
	
	for (int i = 0; i < M; i++)
		sumDocTopic[i] = 0;

	for (int i = 0; i < K; i++)
		sumTopicTerm[i] = 0;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < K; j++)
			nDocTopic[i][j] = 0;

	for (int i = 0; i < K; i++)
		for (int j = 0; j < N; j++)
			nTopicTerm[i][j] = 0;

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			nDocWord[i][j] = 0;

	printf("Space allocated.\n");
}

int LDA::sample(double *distribution, int size)
{
	double randN = (double)(rand() % 999) * 0.001f + 0.001; 
	double sum = 0;
	int i;
	for(i=0; i<size; i++)
	{
		sum = sum + distribution[i];
		if(randN <= sum)
			return i;
	}
	return i;
}

void LDA::loadData(char* filename)
{
	ifstream in(filename);
	char node[10];
	int n, i = 0;
	int wordNum;
	int temp;
	
	while (in >> n)
	{
		for (int j = 0; j<n; j++)
		{
			in >> node;
			sscanf_s(node, "%d:%d", &wordNum, &temp);
			nDocWord[i][wordNum] = temp;
		}
		i++;
	}
	printf("Data loaded.\n");
}

void LDA::calcNewDist(int t, int m)
{
	double sum = 0;
	for(int i=0; i<K; i++)
	{
		newDist[i] = (nTopicTerm[i][t] + beta) / (sumTopicTerm[i] + beta * N) * (nDocTopic[m][i] + alpha);
		sum = sum + newDist[i];
	}
	for (int i=0; i<K; i++)
		newDist[i] = newDist[i] / sum;
}

void LDA::GibbsSampling()
{
	//initialisation - initialise the counts and sums
	double *mult;
	mult = new double[K];
	for(int i=0; i<K; i++)
		mult[i] = 1.0/K;

	int n = 0;
	for(int m=0; m<M; m++)
	{
		for(int t=0; t<N; t++)
		{
			if(nDocWord[m][t]>0)
			{
				for(int i=0, k; i<nDocWord[m][t]; i++)
				{
					k = sample(mult, K);
					nDocTopic[m][k]++;
					sumDocTopic[m]++;
					nTopicTerm[k][t]++;
					sumTopicTerm[k]++;
					kDocWord[m][n] = k;
					tDocWord[m][n] = t;
					n++;
				}
			}		
		}
		cntDocWord[m] = n;
		n = 0;
	}
	printf("Counts and sums initialised.\n");

	//Gibbs Sampling over burn-in period and sampling period
	int iteration = 0;
	while(true)
	{
		printf("iteration %d\n", iteration);

		for(int m=0; m<M; m++)
		{
			for (int n = 0; n < cntDocWord[m]; n++)
			{
				//for the current k and t, decrement counts and sums
				int k, t;
				k = kDocWord[m][n];
				t = tDocWord[m][n];
				nDocTopic[m][k]--;
				nTopicTerm[k][t]--;
				sumTopicTerm[k]--;
				//sample from new distribution of topics
				this->calcNewDist(t, m);
				k = this->sample(newDist, K);
				//for the new topic index, increment counts and sums
				nDocTopic[m][k]++;
				nTopicTerm[k][t]++;
				sumTopicTerm[k]++;
				kDocWord[m][n] = k;
				tDocWord[m][n] = t;
			}
		}
		//check convergence and readout parameters after L iterations
		if(iteration >= L)
		{
			for(int m=0; m<M; m++)
				for(int k=0; k<K; k++)
					thetaDocTopic[m][k] = (nDocTopic[m][k] + alpha) / (sumDocTopic[m] + alpha * K);
			for(int k=0; k<K; k++)
				for(int t=0; t<N; t++)
					phiTopicTerm[k][t] = (nTopicTerm[k][t] + beta) / (sumTopicTerm[k] + beta * N);
			break;
		}
		iteration++;
	}
}

void LDA::printResult(char* outputfile)
{
	ofstream out(outputfile);
	for (int k = 0; k<K; k++)
	{
		for (int j = 0; j<N; j++)
			out << phiTopicTerm[k][j] << " ";
		out << endl;
	}
	out.close();
}