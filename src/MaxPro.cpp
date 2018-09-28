
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <R.h>
# include <Rmath.h>
#define  Min(a,b)     (a<b)?a:b
#define  Max(a,b)     (a>b)?a:b




extern "C" {
	int **SLHD(int m, int t, int k);
	int rc(int n);
	int rc2(int n, int del);
	void distmatrix(int **A, int n, int k, double *d, int s);
	void avgdist(int n, int k, double *d, double *avgdist_cur, int s);
	double combavgdist(int n, int k, double *d, double *avgdist_cur, int s);
	void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s);
	void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old);
	void update_avgdist(int n, int k, double *d, double *avgdist_old, double *avgdist_cur, int s);
	double update_combavgdistI(int n, int k, double *d, double *avgdist_old, double *avgdist_cur, int s);
	void Ddistmatrix(double **A, int n, int k, double *d, int s);
	void Dupdate_distmatrix(double **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s);
	void Ddistmatrix_QQ(double *lambda, double **A, int n, int k, int n_nom, double *d, int s);
	void Dupdate_distmatrix_QQ(double *lambda, double **A, int n, int k, int n_nom, int selrow1, int selrow2, double *d, double *d_old, int s);
	void Ddist_newrow_QQ(double *lambda, double **A, double *B, int n, int k, int n_nom, double *dnewrow, int s);
	void avgdist_newrow(int n, int k, double *dnewrow, double *avgdist_nr, int s);



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
	for (int js = 0; js<t; js++)  
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

int rc(int n) // random int number from 0 to (n-1)
{
	int r;
	double u;
	//GetRNGstate();
	u = unif_rand();
	//PutRNGstate();
	r = (int)(n*u);
	return(r);
}




void distmatrix(int **A, int n, int k, double *d, int s) 
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



void avgdist(int n, int k, double *d, double *avgdist_cur, int s) 
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




double combavgdist(int n, int k, double *d, double *avgdist_cur, int s) 
{
	double combavgdist = 0;
	avgdist(n, k, d, avgdist_cur, s);
	combavgdist = *avgdist_cur;
	return(combavgdist);
}



void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s) 
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



void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old) 
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



void update_avgdist(int n, int k,  double *d, double *avgdist_old, double *avgdist_cur, int s) 
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


	
void Ddistmatrix(double **A, int n, int k, double *d, int s) 
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




void Dupdate_distmatrix(double **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old, int s) 
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


	
void Ddistmatrix_QQ(double *lambda, double **A, int n, int k, int n_nom, double *d, int s) 
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
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[count] += s*log(fabs(double(*(*(A + k3) + k1) - *(*(A + k3) + k2))) + lambda[k3]);
			}

            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + k1) ==*(*(A + k3) + k2) ) {
			    			d[count] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[count] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

			count++;
		}
	}

}



void Dupdate_distmatrix_QQ(double *lambda, double **A, int n, int k, int n_nom, int selrow1, int selrow2, double *d, double *d_old, int s) 
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
            d[position1] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position1] += s*log(fabs(double(*(*(A + k3) + row1) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row1) ==*(*(A + k3) + h) ) {
			    			d[position1] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position1] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

            d[position2] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position2] += s*log(fabs(double(*(*(A + k3) + row2) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row2) ==*(*(A + k3) + h) ) {
			    			d[position2] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position2] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

		}
	}

	for (int h = (row1 + 1); h<row2; h++) //row1<h<row2
	{
		 
		position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
		position2 = (int)row2 + 1 - pow((double)(h + 1), 2)*0.5 + (n - 0.5)*(h + 1) - n - 1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position1] += s*log(fabs(double(*(*(A + k3) + row1) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row1) ==*(*(A + k3) + h) ) {
			    			d[position1] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position1] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

            d[position2] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position2] += s*log(fabs(double(*(*(A + k3) + row2) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row2) ==*(*(A + k3) + h) ) {
			    			d[position2] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position2] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

	}

	if (row2<(n - 1)){
		for (int h = (row2 + 1); h<n; h++) //row1<row2<h
		{
			 
			position1 = (int)h + 1 - pow((double)(row1 + 1), 2)*0.5 + (n - 0.5)*(row1 + 1) - n - 1;
			position2 = (int)h + 1 - pow((double)(row2 + 1), 2)*0.5 + (n - 0.5)*(row2 + 1) - n - 1;
			d_old[position1] = d[position1];
			d_old[position2] = d[position2];
			d[position1] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position1] += s*log(fabs(double(*(*(A + k3) + row1) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row1) ==*(*(A + k3) + h) ) {
			    			d[position1] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position1] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }

            d[position2] = 0;
			for (int k3 = 0; k3<(k-n_nom); k3++)
			{
				d[position2] += s*log(fabs(double(*(*(A + k3) + row2) - *(*(A + k3) + h))) + lambda[k3]);
			}
            if (n_nom>0){
            	for (int k3 = (k-n_nom); k3<k; k3++)
			    {
			    	if ( *(*(A + k3) + row2) ==*(*(A + k3) + h) ) {
			    			d[position2] += s*log(lambda[k3]);
			    	}
			    	else {
			    		d[position2] += s*log(1.0 + lambda[k3]);
			    	}
			    }
            }


		}
	}

}



void Ddist_newrow_QQ(double *lambda, double **A, double *B, int n, int k, int n_nom, double *dnewrow, int s) 
{
	for (int i = 0; i<n; i++)
	{
		dnewrow[i] = 0;
	}
	int count = 0;

	for (int k1 = 0; k1<n; k1++)
	{
				
		for (int k3 = 0; k3<(k-n_nom); k3++)
		{
			dnewrow[count] += s*log(fabs(double(*(*(A + k3) + k1) - B[k3] )) + lambda[k3]);
		}

		if (n_nom>0){
			for (int k3 = (k-n_nom); k3<k; k3++)
			{
				if ( *(*(A + k3) + k1) == B[k3] ) {
					dnewrow[count] += s*log(lambda[k3]);
				}
				else {
					dnewrow[count] += s*log(1.0 + lambda[k3]);
				}
			}
		}

		count++;

	}

}


void avgdist_newrow(int n, int k, double *dnewrow, double *avgdist_nr, int s) 
{
	 
	double Dmin;
	Dmin=dnewrow[1];
	for (int i=0; i<n; i++)
	{
		if(dnewrow[i]<Dmin) Dmin=dnewrow[i];
	}

	double avgdist = 0;

	for (int i = 0; i<n; i++)
	{
		avgdist += exp(double(Dmin-dnewrow[i]));
	}

	avgdist = log(double(avgdist)) - Dmin;

	avgdist = exp(double((avgdist-log(double(n)))*pow((double)(k*s),-1)));

	*avgdist_nr = avgdist;

}




	void MaxProQQ(double *lambda, int *ncol_nom, int *mRow, int *Col, int *localopm, double *initialdesign, int *nstarts, int *IterMax, int *Total_Iter, double *design, double *measure, double *temp0, int *ntotalI, int *svalue)
{
	GetRNGstate();	
    const int k=*Col;                 
	const int m=*mRow;               
	const int t=1;                
	const int n=m*t;
	const int nsearch=*nstarts;             
	const double tfac=0.95;
	const int Imax=*IterMax;             
	int nImax=Min(5*n*(n-1)*k,Imax);   
	int max_itotal=*Total_Iter;               
	
	int s=*svalue;    
	
	int max_itotal1;
	max_itotal1 = max_itotal;


    const int n_nom=*ncol_nom; 

    const int colc=k-n_nom;  

	
	double t0;
	double xcrit;
	double critbest;
	double crittry;
	int itotal;
	double temp;
	double **xbest;
	double **xtry;
	double **x;

	
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


	xtry = new double*[k];
	xbest = new double*[k];
	x = new double*[k];
	for (int i = 0; i<k; i++)
	{
		xbest[i] = new double[n];
		xtry[i] = new double[n];
		x[i] = new double[n];
	}

	
	for (int ii = 0; ii<n; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(xbest + jj) + ii) = *(initialdesign + ii*k + jj);
		}
	}
	

	Ddistmatrix_QQ(lambda, xbest, n, k, n_nom, d, s);
	critbest = combavgdist(n, k, d, avgdist_cur,s);

	t0=*temp0;

	
	itotal = 0;
	for (int isearch = 1; isearch <= nsearch; isearch++)
	{
		
		for (int ii = 0; ii<n; ii++)
		{
			for (int jj = 0; jj<k; jj++)
			{
				*(*(x + jj) + ii) = *(initialdesign + ii*k + jj);
				*(*(xtry + jj) + ii) = *(initialdesign + ii*k + jj);
			}
		}

		Ddistmatrix_QQ(lambda, xtry, n, k, n_nom, d, s);
		xcrit = combavgdist(n, k, d, avgdist_cur, s);
		crittry = xcrit;

		
		temp = t0;
		int ichange = 1;
		
		while (ichange == 1)
		{
			ichange = 0;
			
			int ipert = 1;
			while (ipert<nImax)
			{
				if (itotal>max_itotal1) break;
				itotal = itotal + 1;
				
				int ind;

				int tran1;
				int tran2;
				ind = rc(colc);

				tran1 = rc(m);
				tran2 = rc2(m, tran1);
							
				*(*(xtry + ind ) + tran2) = *(*(x + ind ) + tran1);
				*(*(xtry + ind ) + tran1) = *(*(x + ind ) + tran2);
			    
				Dupdate_distmatrix_QQ(lambda, xtry, n, k, n_nom, tran1, tran2, d, d_old,s);
				crittry = update_combavgdistI(n, k, d,  avgdist_old, avgdist_cur, s);

				if (crittry<critbest)
				{
					
					ichange = 1;
					for (int nn2 = 0; nn2<k; nn2++)
					{
						for (int nn1 = 0; nn1<n; nn1++)
						{
							*(*(xbest + nn2) + nn1) = *(*(xtry + nn2) + nn1);
						}
					}
					
					*(*(x + ind ) + tran1) = *(*(xtry + ind ) + tran1);
					*(*(x + ind ) + tran2) = *(*(xtry + ind ) + tran2);
				    
					critbest = crittry;
					ipert = 1;
					xcrit = crittry;
				}
				else
				{
					
					ipert = ipert + 1;

					if (crittry<xcrit)
					{
						
						*(*(x + ind ) + tran1) = *(*(xtry + ind ) + tran1);
						*(*(x + ind ) + tran2) = *(*(xtry + ind ) + tran2);
					    
						ichange = 1;

						xcrit = crittry;
						
					}
					else
					{
						
						if (*localopm == 0){
							double delta1 = crittry - xcrit;
							double prob = exp(-delta1*pow((double)(temp), (-1)));
							//GetRNGstate();
							double q = unif_rand();
							//PutRNGstate();
							if (prob >= q)
							{
								
								*(*(x + ind ) + tran1) = *(*(xtry + ind ) + tran1);
								*(*(x + ind ) + tran2) = *(*(xtry + ind ) + tran2);
							    
								// ichange=1;
								xcrit = crittry;
							}
							else
							{
							    	
								*(*(xtry + ind ) + tran1) = *(*(x + ind ) + tran1);
								*(*(xtry + ind ) + tran2) = *(*(x + ind ) + tran2);
							    
								revert_distmatrix(n, tran1, tran2, d, d_old);
								*avgdist_cur = *avgdist_old;

							} 
						}
						if (*localopm == 1){
								
								
							*(*(xtry + ind ) + tran1) = *(*(x + ind ) + tran1);
							*(*(xtry + ind ) + tran2) = *(*(x + ind ) + tran2);
						    
							revert_distmatrix(n, tran1, tran2, d, d_old);
							*avgdist_cur = *avgdist_old;

						}

						
					}
				}
			}
			
			if (*localopm == 1) break;
			temp = temp*tfac;
		}
		


	}
	
	
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



void MaxProMeasure(double *lambda, double *design, int *ncol_nom, int *Col, int *n, int *svalue, double *measure)
{

	const int k=*Col;                   
    const int n_nom=*ncol_nom;
    const int nTotal=*n;
    const int s=*svalue;

    double **D;
	D = new double*[k];
	for (int i = 0; i<k; i++)
	{
		D[i] = new double[nTotal];
	}

	for (int ii = 0; ii<nTotal; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(D + jj) + ii) = *(design + ii*k + jj);			 
		}
	}

	int dim = (int)(nTotal*(nTotal - 1)*0.5);
	double *d_all;
	d_all = new double[dim];
	for (int i = 0; i<dim; i++)
	{
		d_all[i] = 0;
	}


	double *avgdist_all = new double;
	*avgdist_all = 0;

	Ddistmatrix_QQ(lambda, D, nTotal, k, n_nom, d_all, s);
	avgdist(nTotal, k, d_all, avgdist_all, s);

	*measure=*avgdist_all;

	for (int i = 0; i<k; i++)
	{
		delete[] D[i];
	}

	delete[]D;

	delete[]d_all;
	
	delete avgdist_all;

}



void MaxProAugment(double *lambda, int *ncol_nom, int *Col, int *nExisting, int *nNewrow, int *nCandidate, double *existingdesign, double *canddesign, double *design, int *svalue, double *measure, int *measurecheck)
{
	
	const int k=*Col;                   
    const int n_nom=*ncol_nom; // number of columns that are nominal factors (as the rightmost columns in the design matrix)
    //const int colc=k-n_nom; // number of columns for continuous and discrete numeric factors
    
    int nExist=*nExisting;
    const int nNew=*nNewrow;
    const int nCand=*nCandidate;
    int nTotal= nExist+nNew;

	int s=*svalue;    
	
    double **CandD;
	CandD = new double*[k];
	for (int i = 0; i<k; i++)
	{
		CandD[i] = new double[nCand];
	}

	for (int ii = 0; ii<nCand; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(CandD + jj) + ii) = *(canddesign + ii*k + jj);			 
		}
	}

    double **ExistD;
	ExistD = new double*[k];
	for (int i = 0; i<k; i++)
	{
		ExistD[i] = new double[nTotal];
	}

	for (int ii = 0; ii<nExist; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(ExistD + jj) + ii) = *(existingdesign + ii*k + jj);			 
		}
	}

	for (int ii = nExist; ii<nTotal; ii++)
	{
		for (int jj = 0; jj<k; jj++)
		{
			*(*(ExistD + jj) + ii) = 0;			 
		}
	}


	double *d_new;
	d_new = new double[nTotal];
	for (int i = 0; i<nTotal; i++)
	{
		d_new[i] = 0;
	}

	int dim = (int)(nTotal*(nTotal - 1)*0.5);
	double *d_all;
	d_all = new double[dim];
	for (int i = 0; i<dim; i++)
	{
		d_all[i] = 0;
	}


	double *avgdist_all = new double;
	*avgdist_all = 0;
	double *avgdist_newr = new double;
	*avgdist_newr = 0;

	double *xrow;
	xrow = new double[k];
	for (int j = 0; j<k; j++)
	{
		xrow[j] = 0;
	}



    int Min_i=0;
	double Min_avgdist=1e+308;
	int viorule=0;
	int flag=0;

	for(int niter=0; niter<nNew; niter++)
	{

		Min_i=0;
		Min_avgdist=1e+308;
		viorule=0;
		flag=0;

		for (int i=0; i<nCand; i++)
		{

			viorule=0;

			for (int j = 0; j<k; j++)
			{
				xrow[j] = *(*(CandD + j) + i);

				if(lambda[j]==0)
				{
					for(int ccc=0; ccc<nExist; ccc++)
					{
						if(xrow[j]== *(*(ExistD + j) + ccc) ) 
						{
							viorule=1;
							break;
						}
					}
				}

				if(viorule==1) break;
			}

			if(viorule==1)
			{
				continue;
			}



			Ddist_newrow_QQ(lambda, ExistD, xrow, nExist, k, n_nom, d_new, s);
			avgdist_newrow(nExist, k, d_new, avgdist_newr, s);

			if(*avgdist_newr<Min_avgdist)
			{
				Min_avgdist=*avgdist_newr;
				Min_i=i;
				flag=1;
			}


		}

		if(flag==0)
		{
			break;	
		}
		else
		{

			for (int j = 0; j<k; j++)
			{
				*(*(ExistD + j) + nExist) = *(*(CandD + j) + Min_i);		 
			}

			nExist=nExist+1;

		}

	}


//Output the final result

	for(int ii=0;ii<nTotal;ii++)
	{
		for(int jj=0;jj<k;jj++)
		{
			*(design+ii*k+jj)=*(*(ExistD+jj)+ii);
		}
	}


	if(flag==1)
	{
		Ddistmatrix_QQ(lambda, ExistD, nTotal, k, n_nom, d_all, s);
		avgdist(nTotal, k, d_all, avgdist_all, s);

		*measure=*avgdist_all;

		*measurecheck=flag;
	}
	else
	{
		*measure=0;

		*measurecheck=flag;
	}
	
	for (int i = 0; i<k; i++)
	{
		delete[] ExistD[i];
		delete[] CandD[i];
	}

	delete[]ExistD;
	delete[]CandD;


	delete[]d_new;
	delete[]d_all;
	delete[]xrow;

	delete avgdist_all;
	delete avgdist_newr;
}



void MaxProLHD(int *mRow, int *Col, int *nstarts, int *IterMax, int *Total_Iter, int *design, double *measure, double *temp0, int *ntotalI, int *svalue)
{
		
	GetRNGstate();
    const int k=*Col;                 
	const int m=*mRow;               
	const int t=1;                
	const int n=m*t;
	const int nsearch=*nstarts;             
	//const int p=1;
	const double tfac=0.95;
	const int Imax=*IterMax;             
	int nImax=Min(5*n*(n-1)*k,Imax);   
	int max_itotal=*Total_Iter;               

	int s=*svalue;    
	
	
	int max_itotal1;
	max_itotal1 = max_itotal;

	double t0;
	double xcrit;
	double critbest;
	double crittry;
	int itotal;
	double temp;
	int **xbest;
	int **xtry;
	int **x;

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


	xtry = new int*[k];
	xbest = new int*[k];
	x = new int*[k];
	for (int i = 0; i<k; i++)
	{
		xbest[i] = new int[n];
		xtry[i] = new int[n];
		x[i] = new int[n];
	}

	
	xbest = SLHD(m, t, k);
	distmatrix(xbest, n, k, d, s);
	critbest = combavgdist(n, k, d, avgdist_cur, s);

	t0=*temp0;

	itotal = 0;
	for (int isearch = 1; isearch <= nsearch; isearch++)
	{
		
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

		temp = t0;
		int ichange = 1;

		while (ichange == 1)
		{
			ichange = 0;

			int ipert = 1;
			while (ipert<nImax)
			{
				if (itotal>max_itotal1) break;
				itotal = itotal + 1;

				int ind;

				int tran1;
				int tran2;
				
				ind = rc(k - 1);			
				tran1 = rc(m);				
				tran2 = rc2(m, tran1);

				*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran1);
				*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran2);

				update_distmatrix(xtry, n, ind + 1, tran1, tran2, d, d_old, s);
				crittry = update_combavgdistI(n, k, d, avgdist_old, avgdist_cur, s);

				if (crittry<critbest)
				{
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
					
					ipert = ipert + 1;

					if (crittry<xcrit)
					{
						
						*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
						*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
						ichange = 1;

						xcrit = crittry;
						
					}
					else
					{
						
						double delta1 = crittry - xcrit;
						double prob = exp(-delta1*pow((double)(temp), (-1)));
						//GetRNGstate();
						double q = unif_rand();
						//PutRNGstate();
						if (prob >= q)
						{
							*(*(x + ind + 1) + tran1) = *(*(xtry + ind + 1) + tran1);
							*(*(x + ind + 1) + tran2) = *(*(xtry + ind + 1) + tran2);
							
							xcrit = crittry;
						}
						else
						{
							*(*(xtry + ind + 1) + tran1) = *(*(x + ind + 1) + tran1);
							*(*(xtry + ind + 1) + tran2) = *(*(x + ind + 1) + tran2);
							revert_distmatrix(n, tran1, tran2, d, d_old);
							*avgdist_cur = *avgdist_old;

						}
						
					}
				}
			}
			temp = temp*tfac;
		}


	}
	

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





}
