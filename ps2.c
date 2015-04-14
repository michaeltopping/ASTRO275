#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float sum(float arr[], int N);
float *multiply(float arr1[], float arr2[], int N);
float f(float halfLifeI, float N, float P, float t);
float findRoot(float t0, float maxErr, int maxIter, float halfLifeI, float N, float P, float dt);
float df(float halfLifeI, float N, float P, float t, float dt);



int main()
{

	int iN = 9;

	float Rb[9] = {0.6843, 0.3895, 0.2243, 0.7609, 0.0227, 0.4666, 0.6631, 0.4606, 0.8904};
	float Sr[9] = {0.7937, 0.7744, 0.7615, 0.7989, 0.7472, 0.7788, 0.7912, 0.7771, 0.8074};
	
	//calculate the slope:
	float ssAge;
	ssAge = (47.5/log(2))*log((sum(multiply(Rb, Sr, iN), iN) - sum(Rb, iN)*sum(Sr, iN)/iN) / 
		   (sum(multiply(Rb, Rb, iN), iN) - pow(sum(Rb, iN), 2)/iN)+1);
 	printf("Age of meteorite: %.3f Gyr\n", ssAge);
	
	
	float t0 = 10;
	float maxErr = 0.0001;
	int maxIter = 10;
	float dt = 0.01;
		
	//Thorium 232
	float halfLifeI = 14.05;
	float N = 2.2222;
	float P = 1.6;
		
	float ThAge =findRoot(t0, maxErr, maxIter, halfLifeI, N, P, dt);



	//Uranium 235
	N = .3246;
	P = 1.16;
	halfLifeI = 0.7038;
		
	float UAge =findRoot(t0, maxErr, maxIter, halfLifeI, N, P, dt);


	printf("Age of the universe (232Th): %f\n", ThAge + ssAge);
	printf("Age of the universe (235U): %f\n", UAge + ssAge);

	return 0;
}




float findRoot(float t0, float maxErr, int maxIter, float halfLifeI, float N, float P, float dt)
{
	float t1, dh;

	for (int iter = 0; iter < maxIter; iter++) {
		dh = f(halfLifeI, N, P, t0) / df(halfLifeI, N, P, t0, dt);
		t1 = t0 - dh;
        printf("At Iteration no. %3d, t = %9.6f\n", iter, t1);
		if (fabs(dh) < maxErr) {
			printf("dh: %f, maxErr: %f\n", fabs(dh), maxErr);
			printf("Maximum accuracy reached\n");
			return t1;
		}
		t0 = t1;
	}
	return 1;
}

//function
float f(float halfLifeI, float N, float P, float t)
{
	float halfLifeU = 4.47; //Gyr
	return (1 - exp(-log(2)*t / halfLifeI))/(1 - exp(-log(2)*t / halfLifeU)) - (N/P)*(halfLifeU/halfLifeI);
}

//derivative of the function
float df(float halfLifeI, float N, float P, float t, float dt)
{
	return (f(halfLifeI, N, P, t+dt) - f(halfLifeI, N, P, t)) / dt;
}

//sum an array
float sum(float arr[], int N)
{
	float total = 0;
	
	for (int ii = 0; ii < N; ii++) {
		total += arr[ii];
	}
	return total;
}

//multiply two arrays together
float *multiply(float arr1[], float arr2[], int N)
{
	float *newArr = (float*)calloc(N, sizeof(float*));
	for(int ii = 0; ii < N; ii++) {
		newArr[ii] = arr1[ii] * arr2[ii];
	}
	
	return newArr;
}