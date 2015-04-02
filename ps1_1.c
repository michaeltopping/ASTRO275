#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


double RAtoRad(float RAh, float RAm, float RAs);
double dectoRad(float dech, float decm, float decs, char pm);
float slope(float *x, float *y, int N);
float intercept(float *x, float *y, float m, int N);
float regress(float *radii, float *vels, float *RAs, float *decs,float Hmin, float Xmin, float Ymin, float Zmin, float Hmax, float Xmax, float Ymax, float Zmax, float dH, float dX, float dY, float dZ, int N);
float det3d(float mat[3][3]);


//MAIN
int main()
{

	//Input the data
	int N = 24;
	float radii[N];
	float vels[N];

	float RAarr[24];
	float decArr[24];

	//open the file for reading
	FILE *inputFile;
	inputFile = fopen("Hub_1929.dat", "r");
	
	float r, v, m, M; //data parameters
	char ID[5]; 
	float RAh, RAm, RAs;
	float dech, decm, decs;
	char h, min, s; //dummy characters
	char pm; //+ or -
	
	//this will be the matrices that hold the data
	//they form the equation y = B a
	//where y is the velocities, a is H and the 3 velocity componenents
	//and be is the data
	float B[24][4];
	float y[24];
	float a[4];
	
	//read in the data and store it as our variables
	for (int ii=0; ii<N; ii++) {
		fscanf(inputFile, "%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", &r, &v, &m, &M, ID, &RAh, &h, &RAm, &min, &RAs, &s, &pm, &dech, &h, &decm, &min, &decs, &s);
		radii[ii] = r;
		vels[ii] = v;
		RAarr[ii] = RAtoRad(RAh, RAm, RAs);
		decArr[ii] = dectoRad(dech, decm, decs, pm);
		B[ii][0] = r;
		B[ii][1] = cos(RAarr[ii])*cos(decArr[ii]);
		B[ii][2] = sin(RAarr[ii])*cos(decArr[ii]);
		B[ii][3] = sin(decArr[ii]);
		y[ii] = v;

	}
	fclose(inputFile);  //close the file

	//print B matrix
// 	for (int row = 0; row < 24; row++) {
// 		for (int column = 0; column < 4; column++) {
// 			printf("%5.2f  ", B[row][column]);
// 		}
// 		printf("\n");
// 	}
	
	
	
	//find transpose of B:
	printf("Calculating transpose of B\n");
	float Bt[4][24];  //B transpose
	for (int row = 0; row < 24; row++) {
		for (int column = 0; column < 4; column++) {
			Bt[column][row] = B[row][column];
		}
	}
	
	
	//multiply Bt and B:
	printf("Multiplying B_transpose by B\n");
	float BtB[4][4]; //B transpose * B
	for (int column=0; column < 4; column++) {
		for (int row=0; row < 4; row++) {
			for (int obj=0; obj < 24; obj++) {
				BtB[row][column] += Bt[row][obj]*B[obj][column];
			}
		}
	}
	
	
	
	//now we find the inverse of BtB
	//we first need the matrix of minors:
	float MoM[4][4];
	float RR[3][3]; //matrix after removing one row and one column
	
	//loop through each element in the matrix of minors:
	printf("Calculating the matrix of minors\n");
	for (int row = 0; row < 4; row++) { //for each row in BtB
		for (int column = 0; column < 4; column++) { //for each column in BtB

			//these variables will check if the next row/column has been removed
			//leaving a 3x3 matrix
			int nextBtBrow = 0;
			int nextBtBcolumn = 0;
			
			//loop through the rows of the new reduced matrix
			for (int rowRR = 0; rowRR < 3; rowRR++) {
				nextBtBcolumn = 0;
				if (nextBtBrow == row) {
					nextBtBrow++;
				}
				//loop through the columns of the reduced matrix
				for (int columnRR = 0; columnRR < 3; columnRR++) {
					if (nextBtBcolumn == column) {
						nextBtBcolumn++;
					} 
					//assign the values of the new reduced matrix
					RR[rowRR][columnRR] = BtB[nextBtBrow][nextBtBcolumn];
					nextBtBcolumn++;
				}
				nextBtBrow++;
			}
		
			//calculate the values for the matrix of minors as the determinant of
			//the reduced matrix
			MoM[row][column] = det3d(RR);
			
		}
	}

	
	//get the matrix of cofactors
	//for odd (row# + column#) the values are multiplied by -1
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			MoM[row][column] = pow(-1, (row+column))*MoM[row][column];
		}
	}
	
	
	//transpose the matrix of minors
	float MoMt[4][24];  //find the transpose of the cofactored matrix of minors
	for (int ii = 0; ii < 24; ii++) {
		for (int jj = 0; jj < 4; jj++) {
			MoMt[jj][ii] = MoM[ii][jj];
		}
	}
	
	
	//find the determinant of the original matrix BtB
	float BtBdet;
	printf("Calculating the determinant of BtB\n");
	//loop through columns of the B transpose * B matrix
	for (int column = 0; column < 4; column++) {

		//for determining which row/column will be removed
		int nextBtBrow = 0;
		int nextBtBcolumn = 0;
		
		//loop through the rows of the reduced matrix
		for (int rowRR = 0; rowRR < 3; rowRR++) {
			nextBtBcolumn = 0;
			if (nextBtBrow == 0) {
				nextBtBrow++;
			}
			for (int columnRR = 0; columnRR < 3; columnRR++) {
				if (nextBtBcolumn == column) {
					nextBtBcolumn++;
				} 
				RR[rowRR][columnRR] = BtB[nextBtBrow][nextBtBcolumn];
				nextBtBcolumn++;
			}
			nextBtBrow++;
		}
	
			
		//calculate the determinant of the BtB matrix
		BtBdet += pow(-1,column)*BtB[0][column]*det3d(RR);
	}
	
	//then divide transpose matrix of minors by the determinant:
	float BtBinv[4][4];
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			BtBinv[row][column] = MoMt[row][column]/BtBdet;

		}
	}
	
	

	//multiply B transpose *B with B:
	printf("Multiplying B_transposeB inverse by B\n");
	float BtBinvBt[4][24];
	for (int kk=0; kk < 4; kk++) {
		for (int jj=0; jj < 4; jj++) {
			for (int ii=0; ii < 24; ii++) {
				BtBinvBt[kk][ii] += BtBinv[kk][jj]*Bt[jj][ii];
			}
		}
	}
	
	
	//multiply by y, the velocities of the galaxies to get a vector of the answers!
	float tot;
	char vars[5] = "HXYZ\0";
	for (int param = 0; param < 4; param++) {
		tot = 0;
		for (int obj = 0; obj < 24; obj++) {
			tot += BtBinvBt[param][obj]*y[obj];
		}
		printf("%c = %f\n", vars[param], tot);
	}
	
	return 0;
}


//inputs a 3x3 matrix and returns the determinant
//just hardcoded
float det3d(float mat[3][3])
{
	float det;
		det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2]) - mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2]) + mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
	return det;
}


//convert RA h, min, sec to radians
double RAtoRad(float RAh, float RAm, float RAs)
{
	float RAdeci;
	RAdeci = RAh * 15 + RAm * (15.0 / 60.0) + RAs * (15.0 / 3600);
	return (RAdeci * M_PI / 180.0);

}

//convert declination deg, min, sec to radians
double dectoRad(float dech, float decm, float decs, char pm)
{
	float decdeci;
	
	decdeci = abs(dech) + decm  / 60.0 + decs  / 3600.0;
	if (pm == '-') {
		decdeci = -decdeci;
	}
	
	return (decdeci * M_PI / 180.0);
}

