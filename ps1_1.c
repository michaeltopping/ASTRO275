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

	FILE *inputFile;
	
	inputFile = fopen("Hub_1929.dat", "r");
	
	float r, v, m, M;
	char ID[5];
	float RAh, RAm, RAs;
	float dech, decm, decs;
	char h, min, s;
	char pm;
	char dec[11];
	float B[24][4];
	float y[24];
	float a[4];
	
	for (int ii=0; ii<N; ii++) {
		fscanf(inputFile, "%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", &r, &v, &m, &M, ID, &RAh, &h, &RAm, &min, &RAs, &s, &pm, &dech, &h, &decm, &min, &decs, &s);
		//printf("%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", r, v, m, M, ID, RAh, h, RAm, min, RAs, s, pm, dech, h, decm, min, decs, s);
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
	fclose(inputFile);


	for (int row = 0; row < 24; row++) {
		for (int column = 0; column < 4; column++) {
			printf("%5.2f  ", B[row][column]);
		}
		printf("\n");
	}
	
	
	
	//find transpose of B:
	printf("Calculating transpose of B\n");
	float Bt[4][24];
	for (int row = 0; row < 24; row++) {
		for (int column = 0; column < 4; column++) {
			Bt[column][row] = B[row][column];
		}
	}
	
	
	printf("Bt:\n");

	
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 24; column++) {
			printf("%5.2f  ", Bt[row][column]);
		}
		printf("\n");
	}
	
	
	//multiply Bt and B:
	printf("Multiplying B_transpose by B\n");
	float BtB[4][4];
	for (int column=0; column < 4; column++) {
		for (int row=0; row < 4; row++) {
			for (int obj=0; obj < 24; obj++) {
				BtB[row][column] += Bt[row][obj]*B[obj][column];
			}
		}
	}
	
	
	
	printf("BtB:\n");
	
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			printf("%f     ", BtB[column][row]);
		}
		printf("\n");
	}
	
	
	//now we find the inverse of BtB
	//we first need the matrix of minors:
	float MoM[4][4];
	float RR[3][3];
	int l, k;
	
	//loop through each element in the matrix of minors:
	printf("Calculating the matrix of minors\n");
	for (int row = 0; row < 4; row++) { //for each row in BtB
		for (int column = 0; column < 4; column++) { //for each column in BtB

			printf("Eliminating row: %d and column: %d\n", row, column);
			int nextBtBrow = 0;
			int nextBtBcolumn = 0;
			for (int rowRR = 0; rowRR < 3; rowRR++) {
				nextBtBcolumn = 0;
				if (nextBtBrow == row) {
					nextBtBrow++;
				}
				for (int columnRR = 0; columnRR < 3; columnRR++) {
					if (nextBtBcolumn == column) {
						nextBtBcolumn++;
					} 
					//printf("Accessing row: %d, column: %d\n", nextBtBrow, nextBtBcolumn);
					RR[rowRR][columnRR] = BtB[nextBtBrow][nextBtBcolumn];
					nextBtBcolumn++;
				}
				nextBtBrow++;
			}
		
				
			printf("Reduced Matrix:\n");
	
			for (int row = 0; row < 3; row++) {
				for (int column = 0; column < 3; column++) {
					printf("%f     ", RR[row][column]);
				}
				printf("\n");
			}
			
			printf("\n");
			MoM[row][column] = det3d(RR);
			
		}
	}

	
	printf("MoM:\n");
	
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			printf("%f     ", MoM[column][row]);
		}
		printf("\n");
	}
	
	
	//get the matrix of cofactors
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			MoM[row][column] = pow(-1, (row+column))*MoM[row][column];
		}
	}
	
	
	//transpose the matrix of minors
	float MoMt[4][24];
	for (int ii = 0; ii < 24; ii++) {
		for (int jj = 0; jj < 4; jj++) {
			MoMt[jj][ii] = MoM[ii][jj];
		}
	}
	
	//find the determinant of the original matrix BtB
	float BtBdet;
	printf("Calculating the determinant of BtB\n");
	for (int column = 0; column < 4; column++) {

		int nextBtBrow = 0;
		int nextBtBcolumn = 0;
		for (int rowRR = 0; rowRR < 3; rowRR++) {
			nextBtBcolumn = 0;
			if (nextBtBrow == 0) {
				nextBtBrow++;
			}
			for (int columnRR = 0; columnRR < 3; columnRR++) {
				if (nextBtBcolumn == column) {
					nextBtBcolumn++;
				} 
				//printf("Accessing row: %d, column: %d\n", nextBtBrow, nextBtBcolumn);
				RR[rowRR][columnRR] = BtB[nextBtBrow][nextBtBcolumn];
				nextBtBcolumn++;
			}
			nextBtBrow++;
		}
	
			
		printf("Reduced Matrix:\n");

		for (int row = 0; row < 3; row++) {
			for (int column = 0; column < 3; column++) {
				printf("%f     ", RR[row][column]);
			}
			printf("\n");
		}
		
		printf("\n");
		
		
		BtBdet += pow(-1,column)*BtB[0][column]*det3d(RR);
	}
	printf("%f\n", BtBdet);
	
	//then divide by the inverse:
	float BtBinv[4][4];
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			BtBinv[row][column] = MoMt[row][column]/BtBdet;

		}
	}
	
	printf("BtBinv:\n");
	
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			printf("%f     ", BtBinv[column][row]);
		}
		printf("\n");
	}
	
	
	
	
	//multiply B and B:
	printf("Multiplying B_transposeB inverse by B\n");
	float BtBinvBt[4][24];
	for (int kk=0; kk < 4; kk++) {
		for (int jj=0; jj < 4; jj++) {
			for (int ii=0; ii < 24; ii++) {
				BtBinvBt[kk][ii] += BtBinv[kk][jj]*Bt[jj][ii];
			}
		}
	}
	
	
	//multiply by y
	float tot;
	for (int param = 0; param < 4; param++) {
		tot = 0;
		for (int obj = 0; obj < 24; obj++) {
			tot += BtBinvBt[param][obj]*y[obj];
		}
		printf("%f\n", tot);
	}
	
	return 0;
}


float det3d(float mat[3][3])
{
	float det;
		det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2]) - mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2]) + mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
	return det;
}



double RAtoRad(float RAh, float RAm, float RAs)
{
	float RAdeci;
	RAdeci = RAh * 15 + RAm * (15.0 / 60.0) + RAs * (15.0 / 3600);
	return (RAdeci * M_PI / 180.0);

}

double dectoRad(float dech, float decm, float decs, char pm)
{
	float decdeci;
	
	decdeci = abs(dech) + decm  / 60.0 + decs  / 3600.0;
	if (pm == '-') {
		decdeci = -decdeci;
	}
	
	return (decdeci * M_PI / 180.0);
}



float regress(float *radii, float *vels, float *RAs, float *decs,float Hmin, float Xmin, float Ymin, float Zmin, float Hmax, float Xmax, float Ymax, float Zmax, float dH, float dX, float dY, float dZ, int N)
{
	float Rsqr = 0;
	float minRes= DBL_MAX;
	float Hlow = 0;
	float Xlow = 0;
	float Ylow = 0;
	float Zlow = 0;
	
	//loop through all possible values
	//H ~ 465
	//X ~ -70
	//Y ~ 240
	//Z ~ -200
	for(float H = Hmin; H < Hmax; H += dH) {
		//printf("H is now: %f\n", H);
		for(float X = Xmin; X < 60; X += dX) {
			//printf("X is now: %f\n", X);

			for(float Y = Ymin; Y < Ymax; Y += dY) {
				//printf("Y is now: %f\n", Y);

				for(float Z = Zmin; Z < -190; Z += dZ) {
					Rsqr = 0;
					for(int ii = 0; ii < N; ii++) {

						Rsqr += pow((vels[ii] - H*radii[ii] - X*cos(RAs[ii])*cos(decs[ii]) - Y*sin(RAs[ii])*cos(decs[ii]) - Z*sin(decs[ii])),2);
					
					}
					if(Rsqr < minRes) {
						minRes = Rsqr;
						Hlow = H;
						Xlow = X;
						Ylow = Y;
						Zlow = Z;
					}
				}
			}
		}
	}
	printf("minRes: %f with H, X, Y, Z: %f %f %f %f\n",minRes, Hlow, Xlow, Ylow, Zlow);
	regress(radii, vels, RAs, decs, Hlow-2*dH, Xlow-2*dX, Ylow-2*dY, Zlow-2*dZ, Hlow+2*dH, Xlow+2*dX, Ylow+2*dY, Zlow+2*dZ, dH/10., dX/10., dY/10., dZ/10., N);

	return 0.0;
}


