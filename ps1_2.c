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
float det4d(float **matrix);
float **transpose( int rows, int columns, float **matrix);
//float **multiply( int rows, int columns, float **matrix1, float **matrix2);
void free_board(float **matrix, int rows);
float **invert( int size, float **matrix);
float **oddIndexNegative( int size, float **matrix);
float **multiply(int iRows, int iColumns, int jColumns, float **matrix1, float **matrix2);


//MAIN
int main()
{

	//Input the data
	int N = 24;
	int rows = 24;
	int columns = 4;
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
	float y[24];
	float a[4];
	
	//create array that data will be stored in
    float **B = (float **)malloc(rows * sizeof(float *)); 

	printf("Dynamically allocated array of column pointers\n");
	//columns
    for (int row = 0; row < rows; row++) {
        B[row] = (float *)malloc(columns * sizeof(float));
    }
	
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


	
	
	//find transpose of B:
	printf("Calculating transpose of B\n");
 	float **Bt = transpose(rows, columns, B);

	//multiply Bt and B:
	printf("Multiplying B_transpose by B\n");
	float **BtB = multiply(rows, columns, columns, B, Bt);
	
	//find the inverse of BtB
	printf("Calculating inverse\n");
	float **BtBinverse = invert(columns, BtB);
	
	//float **mult = multiply(columns, columns, rows, BtBinverse, Bt);
	//multiply B transpose *B with B:
	printf("Multiplying B_transposeB inverse by B\n");
	float BtBinvBt[4][24];
	for (int kk=0; kk < 4; kk++) {
		for (int jj=0; jj < 4; jj++) {
			for (int ii=0; ii < 24; ii++) {
				BtBinvBt[kk][ii] += BtBinverse[kk][jj]*Bt[jj][ii];
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
	
	
	//free memory
	free_board(Bt, columns); 
	free_board(BtB, columns);
	free_board(BtBinverse, columns);
	free_board(B, rows);
	
	
	return 0;
}

//find the transpose of a matrix
float **transpose( int rows, int columns, float **matrix)
{
	printf("Transpose Here\n");
	//rows
    float **Bt = (float **)malloc(columns * sizeof(float *)); 

	printf("Dynamically allocated array of column pointers\n");
	//columns
    for (int column = 0; column < columns; column++) {
        Bt[column] = (float *)malloc(rows * sizeof(float));
    }
    printf("Dynamically allocated arrays of row pointers\n");
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
    float **BtB = (float **)malloc(iColumns * sizeof(float *)); 

	printf("Dynamically allocated array of column pointers\n");
	//columns
    for (int row = 0; row < iRows; row++) {
        BtB[row] = (float *)malloc(iColumns * sizeof(float));
    }
    printf("Dynamically allocated arrays of row pointers\n");

	//float BtB[4][4]; //B transpose * B
	for (int column=0; column < jColumns; column++) {
		for (int row=0; row < iColumns; row++) {
			for (int obj=0; obj < iRows; obj++) {
				BtB[row][column] += matrix2[row][obj]*matrix1[obj][column];
			}
		}
	}
	

	return BtB;
}



//find the inverse of a matrix
float **invert( int size, float **matrix)
{

	float det = det4d(matrix);

	float element;
	printf("Transpose Here\n");
	//float Bt[4][24];  //B transpose
	//rows
    float **inverse = (float **)malloc(size * sizeof(float *)); 

	printf("Dynamically allocated array of column pointers\n");
	//columns
    for (int row = 0; row < size; row++) {
        inverse[row] = (float *)malloc(size * sizeof(float));
    }
    printf("Dynamically allocated arrays of row pointers\n");



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



	//inverse = oddIndexNegative( size, inverse);

	return transpose(size, size, oddIndexNegative( size, inverse));
}




//multiply two matrices
float **oddIndexNegative( int size, float **matrix)
{
	printf("Transpose Here\n");
	//float Bt[4][24];  //B transpose
	//rows
    float **resMat = (float **)malloc(size * sizeof(float *)); 

	printf("Dynamically allocated array of column pointers\n");
	//columns
    for (int row = 0; row < size; row++) {
        resMat[row] = (float *)malloc(size * sizeof(float));
    }
    printf("Dynamically allocated arrays of row pointers\n");

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


float det4d(float **matrix)
{
	//find the determinant of the original matrix BtB
	float det;
	float RR[3][3];
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
				RR[rowRR][columnRR] = matrix[nextBtBrow][nextBtBcolumn];
				nextBtBcolumn++;
			}
			nextBtBrow++;
		}
	
			
		//calculate the determinant of the BtB matrix
		det += pow(-1,column)*matrix[0][column]*det3d(RR);
	}
	
	printf("BtB determinant: %f", det);
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


// free dynamically allocated memory
void free_board(float **matrix, int rows) 
{
    int row;

    // first free each row
    for (row = 0; row < rows; row++) {
         free(matrix[row]);
    }

    // Eventually free the memory of the pointers to the rows
    free(matrix);
 }