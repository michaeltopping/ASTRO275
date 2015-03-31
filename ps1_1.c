#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


double RAtoRad(float RAh, float RAm, float RAs);
double dectoRad(float dech, float decm, float decs, char pm);
float slope(float *x, float *y, int N);
float intercept(float *x, float *y, float m, int N);
float regress(float *radii, float *vels, float *RAs, float *decs, int N);


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
	
	for (int ii=0; ii<N; ii++) {
		fscanf(inputFile, "%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", &r, &v, &m, &M, ID, &RAh, &h, &RAm, &min, &RAs, &s, &pm, &dech, &h, &decm, &min, &decs, &s);
		//printf("%f %f %f %f %s %f %c %f %c %f %c %c %f %c %f %c %f %c", r, v, m, M, ID, RAh, h, RAm, min, RAs, s, pm, dech, h, decm, min, decs, s);
		radii[ii] = r;
		vels[ii] = v;
		RAarr[ii] = RAtoRad(RAh, RAm, RAs);
		decArr[ii] = dectoRad(dech, decm, decs, pm);
		
	}
	fclose(inputFile);

	//Hr + Xcosαcosδ + Ysinαcosδ + Zsinδ
	regress(radii, vels, RAarr, decArr, N);

	
	return 0;
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



float regress(float *radii, float *vels, float *RAs, float *decs, int N)
{
	float Rsqr = 0;
	float minRes= DBL_MAX;
	float Hmax = 500;
	float Xmax = 100;
	float Ymax = 275;
	float Zmax = 250;
	float Hmin = 0;
	float Xmin = 0;
	float Ymin = 0;
	float Zmin = 0;
	
	//loop through all possible values
	//H ~ 465
	//X ~ -70
	//Y ~ 240
	//Z ~ -200
	for(float H = 400.0; H < Hmax; H += 1) {
		printf("H is now: %f\n", H);
		for(float X = -Xmax; X < 0; X += 1) {
			//printf("X is now: %f\n", X);

			for(float Y = 200; Y < Ymax; Y += 1) {
				//printf("Y is now: %f\n", Y);

				for(float Z = -Zmax; Z < -150; Z += 1) {
					Rsqr = 0;
					for(int ii = 0; ii < N; ii++) {

						Rsqr += pow((vels[ii] - H*radii[ii] - X*cos(RAs[ii])*cos(decs[ii]) - Y*sin(RAs[ii])*cos(decs[ii]) - Z*sin(decs[ii])),2);
					
					}
					if(Rsqr < minRes) {
						minRes = Rsqr;
						Hmin = H;
						Xmin = X;
						Ymin = Y;
						Zmin = Z;
					}
				}
			}
		}
	}
	printf("minRes: %f with H, X, Y, Z: %f %f %f %f\n",minRes, Hmin, Xmin, Ymin, Zmin);
	return 0.0;
}


