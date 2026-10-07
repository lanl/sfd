// Version 5.4 09/25/25
// Cleaned up from Version 5.0 01/21/21

//***************************************************************************************************

#include "BoilerPlate.h"
#include "CGData.h"

double eijk[27] = {0., 0., 0., 0., 0., 1., 0., -1., 0., 0., 0., -1., 0., 0., 0., 1., 0., 0., 0., 1., 0., -1., \
	0, 0, 0, 0, 0};

double Sum_term(int j1, int m1, int j2, int m2, int j, int m);

//***************************************************************************************************

void Factorial(){
	// calculate factorials out to 3*JMAX + 1
	
	int i;
	
	Fac[0] = 1.;
	for (i = 1; i < 3*JMAX + 2; i++){
		Fac[i] = i*Fac[i-1];
	}
}

//***************************************************************************************************

double CG_calc(int j1, int m1, int j2, int m2, int j, int m){
	// Calculate the 3D Clebsch Gordan Coefficients
	// See Edmonds "Angular Momentum in Quantum Mechanics" eq 3.6.11
	// Values are zero for impossible combinations
    // Phase convention is the same as Edmonds, which is the same as in Mathematica
	
	double cg = 0.0; // value is zero if indices are out of physical bounds

	if( (j1>=0) && (j2>=0) && (j>=0) && (abs(m1)<=j1) && (abs(m2)<=j2) && (abs(m)<=j) && (m1+m2==m) && (j1+j2>=j) && (abs(j1-j2)<=j)){ // If the indices are physical, calculate the CG

		cg = sqrt((2*j+1)*Fac[j1+j2-j]*Fac[j1-j2+j]*Fac[j2-j1+j]/Fac[j1+j2+j+1])*
				sqrt(Fac[j1+m1]*Fac[j1-m1]*Fac[j2+m2]*Fac[j2-m2]*Fac[j+m]*Fac[j-m])*
				Sum_term(j1,m1,j2,m2,j,m);
	}
	return cg;
}

//***************************************************************************************************

double Phase(int n){
	// This function is the phase factor (-1)ˆn
	
	double phase;
	
	if(fmod(n,2)==0){
		phase=1.;
	}
	else{
		phase=-1.;
	}
	
	return phase;
}

//***************************************************************************************************

double Sum_term(int j1, int m1, int j2, int m2, int j, int m){
	// Summation term for Edmonds eq 3.6.11
	
	double sum_coeff;
	int z, zmin, zmax;
	
	double sum_term = 0.0;
	
	zmin = max(j2-j-m1, max(j1+m2-j, 0));
	zmax = min(j1+j2-j, min(j1-m1, j2+m2));
	
	for(z=zmin;z<zmax+1;z++){
		sum_coeff = Phase(z)/(Fac[z]*Fac[j1+j2-j-z]*Fac[j1-m1-z]*Fac[j2+m2-z]*Fac[j-j2+m1+z]*Fac[j-j1-m2+z]);
		sum_term += sum_coeff;
	}
	
	return sum_term;
}

//***************************************************************************************************

double Mag(double *vect, int length){

	int i;
	double mag;
	
	mag = 0.0;
	for (i = 0; i < length; i++){
		mag += vect[i]*vect[i];
	}
	
	return sqrt(mag);
}

//***************************************************************************************************

double Mag_c(complex<double> *Vc, int length){
	// takes advantage of Vc[i] = conj(Vc[-i])
	
	int i;
	double mag;
	
	mag = norm(Vc[0]);
	for (i = 1; i < length+1; i++){
		mag += 2.*norm(Vc[i]);
	}
	
	return sqrt(mag);
}

//***************************************************************************************************

double Trace(double *vect, int length){
// trace of a square matrix of size length
	int i;
	double tr;
	
	tr = 0.0;
	for (i = 0; i < length; i++){
		tr += vect[i*(length+1)];
	}
	
	return tr;
}

//***************************************************************************************************

double V1dotV2(double *vect1, double *vect2, int length){
	// simple dot product between two vectors of size length
	int i;
	double dot;
	
	dot = 0.0;
	for (i = 0; i < length+1; i++){
		dot += vect1[i]*vect2[i];
	}
	
	return dot;
}

//***************************************************************************************************

double V1dotV2_c(complex<double> *Vect1, complex<double> *Vect2, int J){
	// CG dot product between two vectors of momentum J
	// differs from "normal" dot product by CG coeff: (2J + 1)^-1/2
	// takes advantage of Vc[-m] = (-1)^m*conj(Vc[m])

	int m;
	complex<double> dotc, temp;
	
	// mj0 = m + J; m = 0
	dotc = CG[J][J][J][J][0][0]*Phase(J)*real(Vect1[0]*Vect2[0]);
	for (m = 1; m < J+1; m++){
		temp = CG[J][J+m][J][J-m][0][0]*Phase(J-m)*(Vect1[m]*conj(Vect2[m]) + Vect2[m]*conj(Vect1[m]));
		dotc += temp;
	}
	return real(dotc);
}

//***************************************************************************************************

void ClebschGordan(){
	// calculate Clebsch-Gordan coefficients
	// shift m indices by j so all matrix indices are non-negative
	
	FILE *ofile;
	
	int j, j1, j2, m, m1, m2;
	int mj0, mj1, mj2;
	
	Factorial();
	
	for(j = 0; j < 2*JMAX+1; j++){
		for(m = -j; m < j+1; m++){
			mj0 = m + j;
			for(j1 = 0; j1 < JMAX+1; j1++){
				for(j2 = 0; j2 < JMAX+1; j2++){
					for(m1 = -j1; m1 < j1+1; m1++){
						mj1 = m1 + j1;
						for(m2 = -j2; m2 < j2+1; m2++){
							mj2 = m2 + j2;
							CG[j1][mj1][j2][mj2][j][mj0] = CG_calc(j1, m1, j2, m2, j, m);
						}
					}
				}
			}
		}
	}
	
}

//***************************************************************************************************

void V1outV2_c(complex<double> *Vect1, int J1, complex<double> *Vect2, int J2, complex<double> *Vect, int J){
	// rank J outer product of two vector of length J1 and J2
	// only need to calculate for m >= 0; for m < 0 can use complex conjugate of |m|
	
	int m, m1, m2, mj0, mj1, mj2;
	complex<double> func1, func2, temp;
	
	for (m = 0; m < J + 1; m++){
		mj0 = m + J; // for indexing CG coeffs
		//cout << "m = " << m << endl;
		Vect[m] = 0.;
		for (m1 = -J1; m1 < J1 + 1; m1++){
			m2  = m - m1;
			mj1 = m1 + J1; // for indexing CG coeffs
			mj2 = m2 + J2; // for indexing CG coeffs
			if((m2>=-J2)&&(m2<=J2)){
				if(m1<0){
					func1 = Phase(-m1)*conj(Vect1[-m1]);
				}else{
					func1 = Vect1[m1];
				}
				if(m2<0){
					func2 = Phase(-m2)*conj(Vect2[-m2]);
				}else{
					func2 = Vect2[m2];
				}
				temp     = CG[J1][mj1][J2][mj2][J][mj0]*func1*func2;
				Vect[m] += temp;
			}
		}
	}
}

//***************************************************************************************************

