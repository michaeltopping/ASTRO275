#include <stdio.h>
#include <stdlib.h>

struct Angles{
	float theta;
	float phi;
};

struct Angles sexagesimal(int RAh, int RAm, float RAs, int dech, int decm, float decs);
float avg(float *arr, int size);
float slope(float *x, float *y, int N);
float intercept(float *x, float *y, float m, int N);



//MAIN
int main()
{

	//Input the data
	float radii[] = {0.032, 0.034, 0.214, 0.263, 0.275, 0.275, 0.45, 0.5, 0.5, 0.63, 0.8, 0.9, 0.9, 0.9, 0.9, 1.0, 1.1, 1.1, 1.4, 1.7, 2.0, 2.0, 2.0, 2.0};
	float vels[] = {170, 290, -130, -70, -185, -220, 200, 290, 270, 200, 300, -30, 650, 150, 500 ,920, 450, 500, 500, 960, 500, 850, 800, 1090};
	float mags[] = {1.5, 0.5, 9.0};
	float absMags[] = {-16.0, -17.2, -12.7};
	int N = 24;

	struct Angles angles = sexagesimal(0, 51, 0, -73, 6, 0);

	float m = slope(radii, vels, N);
	float b = intercept(radii, vels, m, N);
	printf("Slope of: %f\n", m);
	printf("Intercept of: %f\n", b);
	return 0;
}



//gets RA and dec and returns angles to the object
struct Angles sexagesimal(int RAh, int RAm, float RAs, int dech, int decm, float decs) 
{
	float RAdeci;
	float decdeci;
	
	RAdeci = RAh * 15 + RAm * (15.0 / 60.0) + RAs * (15.0 / 3600);
	decdeci = abs(dech) + decm  / 60.0 + decs  / 3600.0;
	if (dech < 0) {
		decdeci = -decdeci;
	}
	
	struct Angles angles;
	angles.theta = decdeci;
	angles.phi = RAdeci;
	
	return angles;
}



//calculate the average of a list of numbers
float avg(float *arr, int size)
{
	float total = 0;
	for (int ii = 0; ii < size; ii++) {
		total += arr[ii];
	}
	float average = total/size;
	return average;
}


//calculate the slope from two arrays
float slope(float *x, float *y, int N)
{
	//create an array of x*y
	float xy[N];
	for (int ii = 0; ii < N; ii ++) {
		xy[ii] = x[ii]*y[ii];
	}

	//create an array of x**2
	float xx[N];
	for (int ii = 0; ii < N; ii ++) {
		xx[ii] = x[ii]*x[ii];
	}

	//this gives the slope m
	float m = (avg(xy, N) - avg(x,N)*avg(y,N)) / (avg(xx,N) - avg(x,N)*avg(x,N));
	
	return m;
}


//calculate the y-intercept from two arrays and the least squares regression slope
float intercept(float *x, float *y, float m, int N)
{
	float b = avg(y,N) - m*avg(x,N);
	return b;
}