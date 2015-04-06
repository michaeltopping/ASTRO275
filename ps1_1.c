#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

double RAtoRad(float RAh, float RAm, float RAs);
double dectoRad(float dech, float decm, float decs, char pm);
float slope(float *x, float *y, int N);
float intercept(float *x, float *y, float m, int N);
float regress(float *radii, float *vels, float *RAs, float *decs,float Hmin, float Xmin, float Ymin, float Zmin, float Hmax, float Xmax, float Ymax, float Zmax, float dH, float dX, float dY, float dZ, int N);
float det3d(float mat[3][3]);
float det4d(float **matrix);
float **transpose( int rows, int columns, float **matrix);
void free_board(float **matrix, int rows);
float **invert( int size, float **matrix);
float **oddIndexNegative( int size, float **matrix);
float **multiply(int iRows, int iColumns, int jColumns, float **matrix1, float **matrix2);
void findFit(char filename[], int Nobj, int Nparam);
float findChiSqr(float fit[], float v[], float **B, int Nobj, int Nparam);
float chiSqr(float fit[], float v[], float **B, float vErrs[], int Nobj, int Nparam);
float getVelError(float dErr, float H, float D);


//MAIN
int main()
{
	//problem 1
	findFit("Hub_1929.dat", 24, 4);
	findFit("RPK_SN1a.dat", 20, 4);
	return 0;
}


void findFit(char filename[], int Nobj, int Nparam)
{
	//Input the data
	int N = Nobj;
	int rows = Nobj;
	int columns = Nparam;

	//open the file for reading
	FILE *inputFile;
	inputFile = fopen(filename, "r");
	
	float r, v, m, M; //data parameters
	char ID[5];   //to read in the object ID
	float RAh, RAm, RAs;  //RA vars
	float dech, decm, decs;  //dec vars
	char h, min, s; //dummy character variables
	char pm; //+ or -
	
	//this will be the matrices that hold the data
	//they form the equation y = B a
	//where y is the velocities, a is H and the 3 velocity componenents
	//and be is the data
	float y[Nobj];
	float b[Nobj];
	float a[Nparam];
	float vErrs[Nobj];

	//create array of rows that data will be stored in
    float **B = (float **)calloc(rows,  sizeof(float *)); 

	//columns
    for (int row = 0; row < rows; row++) {
        B[row] = (float *)calloc(columns, sizeof(float));
    }
    
    //check which data file we are looking at
    printf("%d\n", strcmp(filename, "Hub_1929.dat"));
	if (strcmp(filename, "Hub_1929.dat") == 0) {
		printf("Reading in Hub_1929.dat\n");

		//read in the data and store it as our variables
		for (int ii=0; ii<N; ii++) {
			//parse the data
			fscanf(inputFile, "%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", &r, &v, &m, &M, ID, &RAh, &h, &RAm, &min, &RAs, &s, &pm, &dech, &h, &decm, &min, &decs, &s);

			B[ii][0] = r;
			B[ii][1] = cos(RAtoRad(RAh, RAm, RAs))*cos(dectoRad(dech, decm, decs, pm));
			B[ii][2] = sin(RAtoRad(RAh, RAm, RAs))*cos(dectoRad(dech, decm, decs, pm));
			B[ii][3] = sin(dectoRad(dech, decm, decs, pm));
			y[ii] = v;

		}
		fclose(inputFile);  //close the file
	
	} else if (strcmp(filename, "RPK_SN1a.dat") == 0) {
		printf("Reading in RPK_SN1a.dat\n");
		float logv;
		float distMod;
		float sigmai;
		float delta;
		float Av;
		float Avg;
		float H = 63.87; //had to run this program once to find this
		//read in the data and store it as our variables
		for (int ii=0; ii<N; ii++) {
			//parse the data
			fscanf(inputFile, "%s %f %f %f %f %f %f %f %c %f %c %f %c %c %f %c %f %c %f %c", ID, &logv, &distMod, &sigmai, &delta, &Av, &Avg, &RAh, &h, &RAm, &min, &RAs, &s, &pm, &dech, &h, &decm, &min, &decs, &s);

			B[ii][0] = pow(10, (distMod/5)-5);
			B[ii][1] = cos(RAtoRad(RAh, RAm, RAs))*cos(dectoRad(dech, decm, decs, pm));
			B[ii][2] = sin(RAtoRad(RAh, RAm, RAs))*cos(dectoRad(dech, decm, decs, pm));
			B[ii][3] = sin(dectoRad(dech, decm, decs, pm));
			
			y[ii] = pow(10,logv);
			vErrs[ii] = getVelError(sigmai, H, B[ii][0]);
			printf("velocity is: %f +/- %f\n", y[ii], vErrs[ii]);

		}

		fclose(inputFile);  //close the file		

	} else {
		printf("Failed to read in data\n");
	}
		
	//get the matrix X, which is B_i / sigma_i
	//create array of rows that data will be stored in
    float **X = (float **)calloc(rows,  sizeof(float *)); 

	//columns
    for (int row = 0; row < rows; row++) {
        X[row] = (float *)calloc(columns, sizeof(float));
    }
    
	for (int row = 0; row < rows; row++) {
		for (int param = 0; param < columns; param++) {
			X[row][param] = B[row][param]/vErrs[row];
		}
	}
	
	//find transpose of B:
	printf("Calculating transpose of B\n");
 	float **Bt = transpose(rows, columns, B);
 	float **Xt = transpose(rows, columns, X);

	//multiply Bt and B:
	printf("Multiplying B_transpose by B\n");
	float **BtB = multiply(columns, rows, columns, Bt, B);
	float **XtX = multiply(columns, rows, columns, Xt, X);

	//find the inverse of BtB
	printf("Calculating inverse\n");
	float **BtBinverse = invert(columns, BtB);
	float **XtXinverse = invert(columns, XtX);
	
	//float **mult = multiply(columns, columns, rows, BtBinverse, Bt);
	//multiply B transpose *B with B:
	printf("Multiplying B_transposeB inverse by B\n");
	float **BtBinvBt = multiply(columns, columns, rows, BtBinverse, Bt);
	float **XtXinvXt = multiply(columns, columns, rows, XtXinverse, Xt);

	//multiply by y, the velocities of the galaxies to get a vector of the answers!
	float tot;
	float fit[4];
	char vars[5] = "HXYZ\0";
	for (int param = 0; param < 4; param++) {
		tot = 0;
		for (int obj = 0; obj < 24; obj++) {
			if (strcmp(filename, "Hub_1929.dat") == 0) {
				tot += BtBinvBt[param][obj]*y[obj];
			}
			if (strcmp(filename, "RPK_SN1a.dat") == 0) {
				tot += XtXinvXt[param][obj]*y[obj]/vErrs[obj];
			}
		}
		fit[param] = tot;
		printf("%c = %f\n", vars[param], tot);
	}
	
	float sigma = findChiSqr(fit, y, B, Nobj, Nparam);
	float ChiSqr = chiSqr(fit, y, B, vErrs, Nobj, Nparam);
	printf("Sigma = %f\n", sigma);
	printf("Chi^2 = %f\n", ChiSqr);
	
	float errors[4];
	
	for (int param = 0; param < Nparam; param++) {
		if (strcmp(filename, "Hub_1929.dat") == 0) {
			errors[param] = sqrt(BtBinverse[param][param])*sigma;
		}
		if (strcmp(filename, "RPK_SN1a.dat") == 0) {
			errors[param] = sqrt(XtXinverse[param][param]);
		}
		printf("%c = %f +/- %f\n", vars[param], fit[param], errors[param]);
	}
	
	//free memory
	free_board(Bt, columns); 
	free_board(BtB, columns);
	free_board(BtBinverse, columns);
	free_board(B, rows);
	free_board(BtBinvBt, columns);
	
}
//this finds sigma, assuming chi squared is one per degree of freedom
float findChiSqr(float fit[], float v[], float **B, int Nobj, int Nparam)
{
	float total = 0;
	
	for (int obj = 0; obj < Nobj; obj++) {
		total += pow((v[obj] - fit[0]*B[obj][0] - fit[1]*B[obj][1] - fit[2]*B[obj][2] - fit[3]*B[obj][3]),2);
	}
	
	int Ndof = (Nobj - Nparam - 1);
	
	return sqrt(total/Ndof);
}

//get chiSqr
float chiSqr(float fit[], float v[], float **B, float vErrs[], int Nobj, int Nparam)
{
	float total = 0;
	
	for (int obj = 0; obj < Nobj; obj++) {
		total += pow((v[obj] - fit[0]*B[obj][0] - fit[1]*B[obj][1] - fit[2]*B[obj][2] - fit[3]*B[obj][3]) / vErrs[obj],2);
	}
		
	return (total);
}


//find the transpose of a matrix
float **transpose( int rows, int columns, float **matrix)
{
	//rows
    float **Bt = (float **)calloc(columns, sizeof(float *)); 

	//columns
    for (int column = 0; column < columns; column++) {
        Bt[column] = (float *)calloc(rows, sizeof(float));
    }
    //loop through each row and column, set them to columns and rows
	for (int row = 0; row < columns; row++) {
		for (int column = 0; column < rows; column++) {
			Bt[row][column] = matrix[column][row];
		}
	}
	return Bt;
}


//multiply two matrices
float **multiply(int iRows, int iColumns, int jColumns, float **matrix1, float **matrix2)
{
	//rows
    float **BtB = (float **)calloc(iRows, sizeof(float *)); 

	//columns
    for (int row = 0; row < iRows; row++) {
        BtB[row] = (float *)calloc(jColumns, sizeof(float));
    }
	
	//formula to calculate the matrix multiplication
	for (int row=0; row < iRows; row++) {
		for (int column=0; column < jColumns; column++) {
			for (int obj=0; obj < iColumns; obj++) {
				BtB[row][column] += matrix1[row][obj]*matrix2[obj][column];
			}
		}
	}
	return BtB;
}


//find the inverse of a matrix
float **invert( int size, float **matrix)
{
	//will need the determinant
	float det = det4d(matrix);
	
	//allocate memory
	float element;
	//rows
    float **inverse = (float **)calloc(size, sizeof(float *)); 

	//columns
    for (int row = 0; row < size; row++) {
        inverse[row] = (float *)calloc(size, sizeof(float));
    }

	//now we find the inverse of BtB
	//we first need the matrix of minors:
	float RR[3][3]; //matrix after removing one row and one column
	
	//loop through each element in the matrix of minors:
	printf("Calculating the matrix of minors\n");
	for (int row = 0; row < size; row++) { //for each row in BtB
		for (int column = 0; column < size; column++) { //for each column in BtB

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
					RR[rowRR][columnRR] = matrix[nextBtBrow][nextBtBcolumn];
					nextBtBcolumn++;
				}
				nextBtBrow++;
			}
		
			//calculate the values for the matrix of minors as the determinant of
			//the reduced matrix
			element = det3d(RR) / det;
			//printf("Inside invert: inverted element is: %f\n", element);

			inverse[row][column] = element;
		}
	}

	//then transpose the inverse matrix
	return transpose(size, size, oddIndexNegative( size, inverse));
}


//make the elements with odd (row+column) the opposite sign
float **oddIndexNegative( int size, float **matrix)
{
	//float Bt[4][24];  //B transpose
	//rows
    float **resMat = (float **)calloc(size, sizeof(float *)); 

	//columns
    for (int row = 0; row < size; row++) {
        resMat[row] = (float *)calloc(size, sizeof(float));
    }

	//float BtB[4][4]; //B transpose * B
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			resMat[row][column] = pow(-1, (row+column))*matrix[row][column];
		}
	}	

	return resMat;
}


//inputs a 3x3 matrix and returns the determinant
//just hardcoded
float det3d(float mat[3][3])
{
	float det;
		det = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2]) - mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2]) + mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
	return det;
}

//determinant of a 4d matrix
float det4d(float **matrix)
{
	//find the determinant of the original matrix BtB
	float det;
	float RR[3][3];
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
				RR[rowRR][columnRR] = matrix[nextBtBrow][nextBtBcolumn];
				nextBtBcolumn++;
			}
			nextBtBrow++;
		}
	
			
		//calculate the determinant of the BtB matrix
		det += pow(-1,column)*matrix[0][column]*det3d(RR);
	}
	
	return det;

}

//converts error from distance to velocity
float getVelError(float dErr, float H, float D)
{
	float error;
	
	error = abs(H * D * log(10.0) * dErr / 5);
	
	return error;
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

// free dynamically allocated memory
void free_board(float **matrix, int rows) 
{
    int row;
    for (row = 0; row < rows; row++) {
         free(matrix[row]);
    }
    free(matrix);
 }