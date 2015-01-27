//# include <iostream>
//using namespace std;
//# include <fstream>
# include <math.h>
//# include <iomanip>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <R.h>
# include <Rmath.h>
#define  Min(a,b)     (a<b)?a:b
#define  Max(a,b)     (a>b)?a:b




extern "C" {
	int rc(int n);
	int rc2(int n, int del);
	int **SLHD(int m, int t, int k);
	void distmatrix(int **A, int n, int k, double *d, int s);
	void avgdist(int n, int k, double *d, double *avgdist_cur, int s);
	double combavgdist(int n, int k, double *d, double *avgdist_cur, int s);
	void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s);
	void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old);
	void update_avgdist(int n, int k, double *d, double *avgdist_old, double *avgdist_cur, int s);
	double update_combavgdistI(int n, int k, double *d, double *avgdist_old, double *avgdist_cur, int s);

	void Ddistmatrix(double **A, int n, int k, double *d, int s);
	void Dupdate_distmatrix(double **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s);

	
void Ddistmatrix(double **A, int n, int k, double *d, int s) // To compute the interpoint distance matrix
{
	const int dim = (int)(n*(n - 1)*0.5);
	for (int i = 0; i<dim; i++)
	{
		d[i] = 0;
	}
	int count = 0;

	for (int k1 = 0; k1<(n - 1); k1++)
	{
		for (int k2 = (k1 + 1); k2<n; k2++)
		{
			for (int k3 = 0; k3<k; k3++)
			{
				d[count] += s*log(fabs(double(*(*(A + k3) + k1) - *(*(A + k3) + k2))));
			}

			count++;
		}
	}

}





void Dupdate_distmatrix(double **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s) // To update the interpoint distance matrix
{
	int row1 = Min(selrow1, selrow2);
	int row2 = Max(selrow1, selrow2);
	 
	int position1, position2;

	if (row1>0){
		for (int h = 0; h<row1; h++) //h<row1<row2
		{
			 
			position1 = (int)row1 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			d_old[position1] = d[position1];
			d_old[position2] = d[position2];
			d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
			d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

		}
	}

	for (int h = (row1 + 1); h<row2; h++) //row1<h<row2
	{
		 
		position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
		position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
		d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

	}

	if (row2<(n - 1)){
		for (int h = (row2 + 1); h<n; h++) //row1<row2<h
		{
			 
			position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
			position2 = (int)h + 1 - pow((double)(row2 + 1), 2)*0.5 + (n - 0.5)*(row2 + 1) - n - 1;
			d_old[position1] = d[position1];
			d_old[position2] = d[position2];
			d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
			d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

		}
	}

}





void MaxProLHD(int *mRow, int *Col, int *nstarts, int *IterMax, int *Total_Iter, int *design, double *measure, double *temp0, int *ntotalI, int *svalue)
{
		
	GetRNGstate();
    const int k=*Col;                 // define k: number of factors  <= change here
	const int m=*mRow;               // define m: number of runs in each slice <= change here
	const int t=1;                // define t: number of slices    <= always equals 1 here for LHD
	const int n=m*t;
	const int nsearch=*nstarts;             // number of different random starts
	//const int p=1;
	const double tfac=0.95;
	const int Imax=*IterMax;            // maximum number of tries without improving xbest at each temperature for Simulated Annealing 
	int nImax=Min(5*n*(n-1)*k,Imax);   
	int max_itotal=*Total_Iter;               // maximum number of iterations for each random start

	int s=*svalue;    // s value in the distance definition
	
	
	int max_itotal1;
	max_itotal1 = max_itotal;



	//int ntotal = 0;
	double t0;
	double xcrit;
	double critbest;
	double crittry;
	int itotal;
	double temp;
	int **xbest;
	int **xtry;
	int **x;

	// distance matrix
	int dim = (int)(n*(n - 1)*0.5);
	double *d;
	d = new double[dim];

	double *d_old;
	d_old = new double[dim];
	for (int i = 0; i<dim; i++)
	{
		d_old[i] = 0;
	}


	double *avgdist_cur = new double;
	*avgdist_cur = 0;
	double *avgdist_old = new double;
	*avgdist_old = 0;





	// Initialize xtry, xbest, and x;
	xtry = new int*[k];
	xbest = new int*[k];
	x = new int*[k];
	for (int i = 0; i<k; i++)
	{
		xbest[i] = new int[n];
		xtry[i] = new int[n];
		x[i] = new int[n];
	}

	//////initialized the best design ////////////////
	
	xbest = SLHD(m, t, k);
	distmatrix(xbest, n, k, d, s);
	critbest = combavgdist(n, k, d, avgdist_cur, s);


	/////calculate starting temperature.
	/*
	int www = n - 1;
	double avgd = pow((double)(www), (-2 * p));
	double avgdred = pow((double)(www), (2 - 2 * p))*pow((double)(www - 1), -2);
	//double avgdred = dim*pow((double)(www), (2-2*p))*pow((double)(www-1),-2);
	double delta0 = avgdred - avgd;
	t0 = -delta0*pow((double)(log(0.99)), (-1));
	*temp0 = t0;
	*/
	t0=*temp0;

	////Loop for different random starts//////////
	itotal = 0;
	for (int isearch = 1; isearch <= nsearch; isearch++)
	{
		//////initial design ////////////////
		
		x = SLHD(m, t, k);
		for (int n2 = 0; n2<k; n2++)
		{
			for (int n1 = 0; n1<n; n1++)
			{
				*(*(xtry + n2) + n1) = *(*(x + n2) + n1);
			}
		}
		distmatrix(xtry, n, k, d, s);
		xcrit = combavgdist(n, k, d, avgdist_cur, s);
		crittry = xcrit;

		///////////initialize tempertures and counts///////////////////
		temp = t0;
		int ichange = 1;
		////////// variable temperature loop ////////////////////////
		while (ichange == 1)
		{
			ichange = 0;
			//// constant temperature loop /////////////////////
			int ipert = 1;
			while (ipert<nImax)
			{
				if (itotal>max_itotal1) break;
				itotal = itotal + 1;
				//////// switch to be tried is elements 
				/////// change two component in a column ////////// 
				int ind;

				int tran1;
				int tran2;
				
				ind = rc(k - 1);			
				tran1 = rc(m);				
				tran2 = rc2(m, tran1);
				/////////perturb x to xtry////////////////////
				*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran1);
				*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran2);
				//////////////////////////////////////////////
				update_distmatrix(xtry, n, ind + 1, tran1, tran2, d, d_old, s);
				crittry = update_combavgdistI(n, k, d, avgdist_old, avgdist_cur, s);


				////// is xtry better than xbest? //////////////////////////////
				if (crittry<critbest)
				{
					//////////yes: replace x, xbest by xtry ; set iterp=1;ichange=1////////////////////////
					ichange = 1;
					for (int nn2 = 0; nn2<k; nn2++)
					{
						for (int nn1 = 0; nn1<n; nn1++)
						{
							*(*(xbest + nn2) + nn1) = *(*(xtry + nn2) + nn1);
						}
					}
					*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
					*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
					critbest = crittry;
					ipert = 1;
					xcrit = crittry;
				}
				else
				{
					//////////No:, increase ipert by 1. is xtry better than x?
					ipert = ipert + 1;

					if (crittry<xcrit)
					{
						////// xtry is better than x; replace x by xtry ///////////
						*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
						*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
						ichange = 1;

						xcrit = crittry;
						//////////////////////////////////////////////////////////
					}
					else
					{
						///////// xtry is worst than x////////////////////////////
						double delta1 = crittry - xcrit;
						double prob = exp(-delta1*pow((double)(temp), (-1)));
						//GetRNGstate();
						double q = unif_rand();
						//PutRNGstate();
						if (prob >= q)
						{///// replce x by xtry by prob///////////
							*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
							*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
							// ichange=1;
							xcrit = crittry;
						}///////////////////////////////////
						else
						{///// reset x try to x for the next pertubation		
							*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran1);
							*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran2);
							revert_distmatrix(n, tran1, tran2, d, d_old);
							*avgdist_cur = *avgdist_old;

						}////////////////////////////////////////// 
						//////////////////////////////////////////////////////////
					}
				}
			}
			//// end of constant temperature loop ////////////
			temp = temp*tfac;
		}
		///////// End of variable temperature loop///////////////////


	}
	/////end of search loop////////////////////////////



	
	/// Output the best design for the current loop
        for(int ii=0;ii<n;ii++)
		{
			for(int jj=0;jj<k;jj++)
			{
				*(design+ii*k+jj)=*(*(xbest+jj)+ii);
			}
		}

	PutRNGstate();	
	
	*measure=critbest;
	*ntotalI = itotal;
					

	// delete xtry, xbest, and x;
	for (int i = 0; i<k; i++)
	{
		delete[] xbest[i];
		delete[] xtry[i];
		delete[] x[i];

	}

	delete[]xbest;
	delete[]xtry;
	delete[]x;

	delete[]d;
	delete[]d_old;

	delete avgdist_cur;
	delete avgdist_old;
}



	void MaxProImprove(int *mRow, int *Col, int *localopm, double *initialdesign, int *nstarts, int *IterMax, int *Total_Iter, double *design, double *measure, double *temp0, int *ntotalI, int *svalue)
{
	GetRNGstate();	
    const int k=*Col;                 // define k: number of factors  <= change here
	const int m=*mRow;               // define m: number of runs in each slice <= change here
	const int t=1;                // define t: number of slices    <= always equals 1 here for LHD
	const int n=m*t;
	const int nsearch=*nstarts;             // number of different random starts
	//const int p=1;
	const double tfac=0.95;
	const int Imax=*IterMax;            // maximum number of tries without improving xbest at each temperature for Simulated Annealing 
	int nImax=Min(5*n*(n-1)*k,Imax);   
	int max_itotal=*Total_Iter;               // maximum number of iterations for each random start
	
	int s=*svalue;    // s value in the distance definition
	
	int max_itotal1;
	max_itotal1 = max_itotal;



	//int ntotal = 0;
	double t0;
	double xcrit;
	double critbest;
	double crittry;
	int itotal;
	double temp;
	double **xbest;
	double **xtry;
	double **x;

	// distance matrix
	int dim = (int)(n*(n - 1)*0.5);
	double *d;
	d = new double[dim];

	double *d_old;
	d_old = new double[dim];
	for (int i = 0; i<dim; i++)
	{
		d_old[i] = 0;
	}


	double *avgdist_cur = new double;
	*avgdist_cur = 0;
	double *avgdist_old = new double;
	*avgdist_old = 0;





	// Initialize xtry, xbest, and x;
	xtry = new double*[k];
	xbest = new double*[k];
	x = new double*[k];
	for (int i = 0; i<k; i++)
	{
		xbest[i] = new double[n];
		xtry[i] = new double[n];
		x[i] = new double[n];
	}

	//////initial design as the best design ////////////////
	for (int ii = 0; ii<n; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(xbest + jj) + ii) = *(initialdesign + ii*k + jj);
		}
	}
	

	Ddistmatrix(xbest, n, k, d, s);
	critbest = combavgdist(n, k, d, avgdist_cur,s);


	/////calculate starting temperature.
	/*
	double www = *width;
	double avgd = pow((double)(www), (-2 * p));
	//double avgdred = pow((double)(www), (2 - 2 * p))*pow((double)(www - www/(n-1)), -2);
	double avgdred = pow((double)(www - www/(n-1)), (-2 * p));
	double delta0 = dim*(avgdred - avgd);
	t0 = -delta0*pow((double)(log(0.99)), (-1));
	*temp0 = t0;
	*/
	t0=*temp0;

	////Loop for different random starts//////////
	itotal = 0;
	for (int isearch = 1; isearch <= nsearch; isearch++)
	{
		//////initial design ////////////////
		for (int ii = 0; ii<n; ii++)
		{
			for (int jj = 0; jj<k; jj++)
			{
				*(*(x + jj) + ii) = *(initialdesign + ii*k + jj);
				*(*(xtry + jj) + ii) = *(initialdesign + ii*k + jj);
			}
		}

		Ddistmatrix(xtry, n, k, d, s);
		xcrit = combavgdist(n, k, d, avgdist_cur, s);
		crittry = xcrit;

		///////////initialize tempertures and counts///////////////////
		temp = t0;
		int ichange = 1;
		////////// variable temperature loop ////////////////////////
		while (ichange == 1)
		{
			ichange = 0;
			//// constant temperature loop /////////////////////
			int ipert = 1;
			while (ipert<nImax)
			{
				if (itotal>max_itotal1) break;
				itotal = itotal + 1;
				//////// switch to be tried is elements 
				/////// change two component in a column ////////// 
				int ind;

				int tran1;
				int tran2;
				ind = rc(k - 1);

				tran1 = rc(m);
				tran2 = rc2(m, tran1);
				/////////perturb x to xtry////////////////////
				*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran1);
				*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran2);
				//////////////////////////////////////////////
				Dupdate_distmatrix(xtry, n, ind + 1, tran1, tran2, d, d_old,s);
				crittry = update_combavgdistI(n, k, d,  avgdist_old, avgdist_cur, s);


				////// is xtry better than xbest? //////////////////////////////
				if (crittry<critbest)
				{
					//////////yes: replace x, xbest by xtry ; set iterp=1;ichange=1////////////////////////
					ichange = 1;
					for (int nn2 = 0; nn2<k; nn2++)
					{
						for (int nn1 = 0; nn1<n; nn1++)
						{
							*(*(xbest + nn2) + nn1) = *(*(xtry + nn2) + nn1);
						}
					}
					*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
					*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
					critbest = crittry;
					ipert = 1;
					xcrit = crittry;
				}
				else
				{
					//////////No:, increase ipert by 1. is xtry better than x?
					ipert = ipert + 1;

					if (crittry<xcrit)
					{
						////// xtry is better than x; replace x by xtry ///////////
						*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
						*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
						ichange = 1;

						xcrit = crittry;
						//////////////////////////////////////////////////////////
					}
					else
					{
						///////// xtry is worst than x////////////////////////////
						if (*localopm == 0){
							double delta1 = crittry - xcrit;
							double prob = exp(-delta1*pow((double)(temp), (-1)));
							//GetRNGstate();
							double q = unif_rand();
							//PutRNGstate();
							if (prob >= q)
							{///// replce x by xtry by prob///////////
								*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
								*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
								// ichange=1;
								xcrit = crittry;
							}///////////////////////////////////
							else
							{///// reset x try to x for the next pertubation		
								*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran1);
								*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran2);
								revert_distmatrix(n, tran1, tran2, d, d_old);
								*avgdist_cur = *avgdist_old;

							}////////////////////////////////////////// 
						}
						if (*localopm == 1){
							///// reset x try to x for the next pertubation		
							*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran1);
							*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran2);
							revert_distmatrix(n, tran1, tran2, d, d_old);
							*avgdist_cur = *avgdist_old;

						}

						//////////////////////////////////////////////////////////
					}
				}
			}
			//// end of constant temperature loop ////////////
			if (*localopm == 1) break;
			temp = temp*tfac;
		}
		///////// End of variable temperature loop///////////////////


	}
	/////end of search loop////////////////////////////



	
	/// Output the best design for the current loop
        for(int ii=0;ii<n;ii++)
		{
			for(int jj=0;jj<k;jj++)
			{
				*(design+ii*k+jj)=*(*(xbest+jj)+ii);
			}
		}

		
	
	*measure=critbest;
	*ntotalI = itotal;
					
	PutRNGstate();
	// delete xtry, xbest, and x;
	for (int i = 0; i<k; i++)
	{
		delete[] xbest[i];
		delete[] xtry[i];
		delete[] x[i];

	}

	delete[]xbest;
	delete[]xtry;
	delete[]x;

	delete[]d;
	delete[]d_old;

	delete avgdist_cur;
	delete avgdist_old;
}






int **SLHD(int m, int t, int k)
{
	int te;
	int **SLHD;
	int *r;
	int n = m*t;
	r = new int[m];
	SLHD = new int*[k];
	for (int iii = 0; iii<k; iii++)
	{
		SLHD[iii] = new int[n];
	}
	for (int js = 0; js<t; js++)  // first dimention is 1,2,3,.....
	{
		for (int j2 = 0; j2<m; j2++){
			*(*(SLHD + 0) + m*js + j2) = (j2 + 1);
		}
	}

	for (int j1 = 0; j1<(k - 1); j1++)
	{
		for (int jss = 0; jss<t; jss++){

			for (int cc = 0; cc<m; cc++)
			{
				*(r + cc) = cc + 1;
			}
			for (int c = 0; c<m; c++)
			{
				te = rc(m - c);
				*(*(SLHD + j1 + 1) + jss*m + c) = *(r + te);

				for (int c1 = 0; c1<(m - c - 1); c1++)
				{
					if (c1 >= te)
					{
						*(r + c1) = *(r + c1 + 1);
					}
				}
			}

		}
	}

	int **SLHDS;
	SLHDS = new int*[k];
	for (int iii = 0; iii<k; iii++)
	{
		SLHDS[iii] = new int[n];
	}

	int xsubs;
	for (int j3 = 0; j3<k; j3++){
		for (int j5 = 0; j5<m; j5++){
			xsubs = j5*t + 1;
			for (int j4 = 0; j4<n; j4++){
				if (*(*(SLHD + j3) + j4) == (j5 + 1)){
					*(*(SLHDS + j3) + j4) = xsubs;
					xsubs++;
				}
			}
		}
	}


	return(SLHDS);
	for (int iiii = 0; iiii<k; iiii++)
	{
		delete[] SLHD[iiii];
		delete[] SLHDS[iiii];
	}
	delete[] SLHD;
	delete[] SLHDS;
	delete[] r;
}




int rc2(int n, int del)
{
	int rctwo;

	rctwo = rc(n - 1);
	if (rctwo >= del)  rctwo++;

	return(rctwo);
}

int rc(int n) // choose randomly from 0 to (n-1)
{
	int r;
	double u;
	//GetRNGstate();
	u = unif_rand();
	//PutRNGstate();
	r = (int)(n*u);
	return(r);
}




void distmatrix(int **A, int n, int k, double *d, int s) // To compute the interpoint distance matrix
{
	const int dim = (int)(n*(n - 1)*0.5);
	for (int i = 0; i<dim; i++)
	{
		d[i] = 0;
	}
	int count = 0;

	for (int k1 = 0; k1<(n - 1); k1++)
	{
		for (int k2 = (k1 + 1); k2<n; k2++)
		{
			for (int k3 = 0; k3<k; k3++)
			{
				d[count] += s*log(fabs(double(*(*(A + k3) + k1) - *(*(A + k3) + k2))));
			}

			count++;
		}
	}

}



void avgdist(int n, int k, double *d, double *avgdist_cur, int s) // To compute the average reciprocal interpoint distance
{
	const int dim = (int)(n*(n - 1)*0.5);
	

	double Dmin;
	Dmin=d[1];
	for (int i=0; i<dim; i++)
	{
		if(d[i]<Dmin) Dmin=d[i];
	}

	double avgdist = 0;

	for (int i = 0; i<dim; i++)
	{
		avgdist += exp(double(Dmin-d[i]));
	}

	avgdist = log(double(avgdist)) - Dmin;

	avgdist = exp(double((avgdist-log(double(dim)))*pow((double)(k*s),-1)));

	*avgdist_cur = avgdist;

}




double combavgdist(int n, int k, double *d, double *avgdist_cur, int s) // To compute the combined average reciprocal interpoint distance
{
	double combavgdist = 0;
	avgdist(n, k, d, avgdist_cur, s);
	combavgdist = *avgdist_cur;
	return(combavgdist);


}






void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s) // To update the interpoint distance matrix
{
	int row1 = Min(selrow1, selrow2);
	int row2 = Max(selrow1, selrow2);
	 
	int position1, position2;

	if (row1>0){
		for (int h = 0; h<row1; h++) //h<row1<row2
		{
			 
			position1 = (int)row1 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			d_old[position1] = d[position1];
			d_old[position2] = d[position2];
			d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
			d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

		}
	}

	for (int h = (row1 + 1); h<row2; h++) //row1<h<row2
	{
		 
		position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
		position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
		d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

	}

	if (row2<(n - 1)){
		for (int h = (row2 + 1); h<n; h++) //row1<row2<h
		{
			 
			position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
			position2 = (int)h + 1 - pow((double)(row2 + 1), 2)*0.5 + (n - 0.5)*(row2 + 1) - n - 1;
			d_old[position1] = d[position1];
			d_old[position2] = d[position2];
			d[position1] = d[position1] + s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h))));
			d[position2] = d[position2] + s*log(fabs(double(*(*(A + col) + row2) - *(*(A + col) + h)))) - s*log(fabs(double(*(*(A + col) + row1) - *(*(A + col) + h))));

		}
	}

}




void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old) // To revert the interpoint distance matrix
{
	int row1 = Min(selrow1, selrow2);
	int row2 = Max(selrow1, selrow2);
	int position1, position2;

	if (row1>0){
		for (int h = 0; h<row1; h++) //h<row1<row2
		{
			position1 = (int)row1 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
			d[position1] = d_old[position1];
			d[position2] = d_old[position2];

		}
	}

	for (int h = (row1 + 1); h<row2; h++) //row1<h<row2
	{
		position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
		position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
		d[position1] = d_old[position1];
		d[position2] = d_old[position2];
	}

	if (row2<(n - 1)){
		for (int h = (row2 + 1); h<n; h++) //row1<row2<h
		{
			position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
			position2 = (int)h + 1 - pow((double)(row2 + 1), 2)*0.5 + (n - 0.5)*(row2 + 1) - n - 1;
			d[position1] = d_old[position1];
			d[position2] = d_old[position2];
		}
	}

}






void update_avgdist(int n, int k,  double *d, double *avgdist_old, double *avgdist_cur, int s) // To update the average reciprocal interpoint distance
{
	*avgdist_old = *avgdist_cur;

	const int dim = (int)(n*(n - 1)*0.5);
	
	double Dmin;
	Dmin=d[1];
	for (int i=0; i<dim; i++)
	{
		if(d[i]<Dmin) Dmin=d[i];
	}

	double avgdist = 0;

	for (int i = 0; i<dim; i++)
	{
		avgdist += exp(double(Dmin-d[i]));
	}

	avgdist = log(double(avgdist)) - Dmin;

	avgdist = exp(double((avgdist-log(double(dim)))*pow((double)(k*s),-1)));

	*avgdist_cur = avgdist;

}




double update_combavgdistI(int n, int k, double *d, double *avgdist_old, double *avgdist_cur, int s)
{
	 
	double combavgdist = 0;

	update_avgdist(n, k, d, avgdist_old, avgdist_cur, s);

	combavgdist = *avgdist_cur;

	return(combavgdist);


}





}
