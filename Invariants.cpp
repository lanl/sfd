// Version 5.4 10/31/25

//***************************************************************************************************
// 080918: rescales the 7/31 version, adds radial nodes, general clean-up
// 102918: cleaned up P4 scaling and redefined the internal rotations
// 103125: Cleaning up for release

#include "BoilerPlate.h"
double V1dotV2_c(complex<double> *vect1, complex<double> *vect2, int length);
void V1outV2_c(complex<double> *Vect1, int j1, complex<double> *Vect2, int j2, complex<double> *Vect, int j);

double tiny = 1.e-15;
double Rtrt3 = sqrt(sqrt(3.));
double Rtrt5 = sqrt(sqrt(5.));
double Rtrt7 = sqrt(sqrt(7.));
double srt3  = sqrt(3.);

//***************************************************************************************************

void P1_Invar(complex<double> *P1, double *P1_Inv,
			  complex<double> *Z, complex<double> *X, complex<double> *Y, double *O1_Inv){
	// Finds all the invariants of rank 1 tensor (vector) and 3 orientation factors

    double magi;
    
	P1_Inv[0] = sqrt(V1dotV2_c(P1, P1, 1)); // CG normalization; r^1 units
    magi = 1./max(tiny,P1_Inv[0]);

	// orientations with respect to external frame
    O1_Inv[0] = Rtrt3*V1dotV2_c(P1, Z, 1)*magi;
	O1_Inv[1] = Rtrt3*V1dotV2_c(P1, X, 1)*magi;
	O1_Inv[2] = Rtrt3*V1dotV2_c(P1, Y, 1)*magi;
}

//***************************************************************************************************

void P2_Invar(complex<double> *P2, double *P0, double *P2_Inv,
			  complex<double> *Z2, complex<double> *Y2, complex<double> *X2, double *O2_Inv){
	// Finds the three invariants of a symmtric rank 2 tensor and the 3 orientation factors

    double magi, norm2;
    complex<double> P2_2[3];
 
	// The first invariant is the norm of the order 2 tensor: net non-sphericity of distribution
	norm2 = sqrt(V1dotV2_c(P2, P2, 2));
	P2_Inv[0] = norm2;               // r^2 units
	magi = 1./max(tiny,norm2); // r^-2 units

	// next invariant is the skew metric of the order 2 tensor
	// First make the rank 2 outer product of the rank 2 tensor with itself
	V1outV2_c(P2, 2, P2, 2, P2_2, 2);
	// And dot that back into the original rank 2 tensor: skewness
	// Normalize by (P2•P2) would generate skewness * r2: magnitude of skewness in distance
	P2_Inv[1] = V1dotV2_c(P2,P2_2,2)*magi*magi; // r^2 units

	// The third invariant is the trace scalar defined in the conversion to spherical tensors
	// ****** THIS IS IN r^2 and is set in T4.cpp*********//
	P2_Inv[2] = P0[1]; // r^2 units

	// orientations with respect to external frame: approach based on the second moments
    // And the orientations vs the X2 unit tensors
    O2_Inv[0] = Rtrt5*V1dotV2_c(P2, Z2, 2)*magi;
    O2_Inv[1] = Rtrt5*V1dotV2_c(P2, X2, 2)*magi;
    O2_Inv[2] = Rtrt5*V1dotV2_c(P2, Y2, 2)*magi;

}

//***************************************************************************************************

void P3_Invar(complex<double> *P3a, complex<double> *P1a, double *P3_Inv, complex<double> *P_1,
			  complex<double> *Z3, complex<double> *X3, complex<double> *Y3, double *O3_Inv){
	// Finds the 5 invariants of symmetric rank 3 tensor and the 2 internal and 3 external orientation values

	// Temporary products
    complex<double> P1a_2[3], P3a_2[3], P3a_22[3], P3_33a2[3], P3a_4[5], P3a_44[5];
    double mag3i, norm112, norm3, norm31, norm332, norm333, norm334;
    
	// First invariant is the norm of the order 3 tensor
    //norm3 = P3sc0*V1dotV2_c(P3a, P3a, 3);
	norm3 = sqrt(V1dotV2_c(P3a, P3a, 3));
	P3_Inv[0] = norm3;           // r^3
	mag3i = 1./max(tiny,norm3);  // r^-3

	// Measure the skew of the rank 3 tensor
	// Make rank 2 outer product and make the norm
    V1outV2_c(P3a, 3, P3a, 3, P3a_2, 2); // r^6
    norm332 = sqrt(V1dotV2_c(P3a_2, P3a_2, 2)); // r^6
    P3_Inv[1] = norm332*mag3i;  // r^3
	
	// Make rank 2 outer product again and dot back into original for skewness
	V1outV2_c(P3a_2, 2, P3a_2, 2, P3a_22, 2); // r^12
 	P3_Inv[2] = V1dotV2_c(P3a_2, P3a_22, 2)*mag3i/max(norm332*norm332,tiny); // r^3

    // Make rank 4 outer product and find its skewness
    V1outV2_c(P3a, 3, P3a, 3, P3a_4, 4); // r^6
    norm334 = sqrt(V1dotV2_c(P3a_4, P3a_4, 4)); // r^6
    V1outV2_c(P3a_4, 4, P3a_4, 4, P3a_44, 4); // r^12
    P3_Inv[3] = V1dotV2_c(P3a_4, P3a_44, 4)*mag3i/max(norm334*norm334,tiny); // r^3

	// Fourth invariant is the vector "norm" of the order 1 vectors
	norm31 = sqrt(V1dotV2_c(P1a, P1a, 1));
	P3_Inv[4] = norm31;

	// Now the orientations between the primary 33 matrix and 31 vector using rank 2 matrices
	// make rank 2 outer products of the primary vectors
    V1outV2_c(P1a, 1, P1a, 1, P1a_2, 2); // nu(1,1)2 r^6
    norm112 = sqrt(V1dotV2_c(P1a_2, P1a_2, 2)); // r^6
    V1outV2_c(P3a_2, 2, P3a_2, 2, P3_33a2, 2); // nu(3,3)2 x nu(3,3)2 r^12
    norm333 = sqrt(V1dotV2_c(P3_33a2, P3_33a2, 2)); // r^12

    // make contractions of these for orientations and scale to norm3*norm31
    P3_Inv[5] = V1dotV2_c(P3a_2, P1a_2, 2)*sqrt(norm3*norm31)/max(norm332*norm112,tiny); // r^3
    P3_Inv[6] = V1dotV2_c(P3_33a2, P1a_2, 2)*sqrt(norm3*norm31)/max(norm333*norm112,tiny); // r^3
    
    // And the orientations vs the X3 unit tensors
    O3_Inv[0] = Rtrt7*V1dotV2_c(P3a, Z3, 3)*mag3i;
    O3_Inv[1] = Rtrt7*V1dotV2_c(P3a, X3, 3)*mag3i;
    O3_Inv[2] = Rtrt7*V1dotV2_c(P3a, Y3, 3)*mag3i;

}

//***************************************************************************************************

void P4_Invar(complex<double> *P4a, complex<double> *P2a, double *P0, complex<double> *P2, double *P4_Inv,
			  complex<double> *Z4, complex<double> *X4, complex<double> *Y4, double *O4_Inv){
	// Finds all the invariants of the symmetric rank 4 tensor for positions: 15 DOF

    // Temporary outer products
    complex<double> P2_2[3], P4_2[3], P4_22[3], P4_224[5], P4_4[5], P4_44[5], P4_4d4[5];
    double norm4, mag4i, norm2, magi, norm23, norm44, norm42, mag2i;

    // First invariant is the norm of the order 4 tensors: the rank 0 contraction
    // Traditional norms would be this term *sqrt9 and then taking the sqrt
    norm4 = sqrt(V1dotV2_c(P4a, P4a, 4)); // r^4
    P4_Inv[0] = norm4;
    mag4i = 1./max(norm4,tiny);

    // Rank 2 outer product of the rank 4 tensor with itself
    // Norm and skewness of rank 2 outer product
    V1outV2_c(P4a, 4, P4a, 4, P4_2, 2); // r^8
    norm2 = sqrt(V1dotV2_c(P4_2,P4_2,2)); // r^8
    magi  = 1./max(tiny,norm2);
    P4_Inv[1] = norm2*mag4i;              // r^4
    V1outV2_c(P4_2, 2, P4_2, 2, P4_22, 2); // r^16
    norm23 = V1dotV2_c(P4_2,P4_22,2); // r^24
    P4_Inv[2] = norm23*mag4i*magi*magi; // r^4

    // Rank 4 contractions: make skewness metric of initial tensor
    V1outV2_c(P4a, 4, P4a, 4, P4_4, 4); // r^8 4x4
    V1outV2_c(P4a, 4, P4_4, 4, P4_44, 4); // r^12 4x4x4

    // Rank 4 contractions: make 4th order contraction of initial tensor
    V1outV2_c(P4_4, 4, P4_4, 4, P4_4d4, 4); // r^16 (4x4)x(4x4)
    norm44 = max(tiny,V1dotV2_c(P4_4,P4_4,4)); // r^16
    P4_Inv[3] = V1dotV2_c(P4_4,P4_4d4,4)*mag4i/norm44;  // r^4 (4x4).((4x4)x(4x4))

    // 5th order contraction
    V1outV2_c(P4_2, 2, P4_2, 2, P4_224, 4); // r^16
    P4_Inv[4] = V1dotV2_c(P4a,P4_224,4)*mag4i*magi*sqrt(magi); // r^4

    // 7th order contraction
    P4_Inv[5] = V1dotV2_c(P4_44,P4_4d4,4)*mag4i*mag4i/norm44; // r^4 (4x4x4).((4x4)x(4x4))

    // next simple invariant is the norm of the rank 2 subspace tensor
    norm42 = sqrt(V1dotV2_c(P2a, P2a, 2)); // r^4
    P4_Inv[6] = norm42;

    // Dot that back into the original rank 2 tensor and normalize: skewness
    V1outV2_c(P2a, 2, P2a, 2, P2_2, 2);
    mag2i = 1./max(norm42,tiny);
    P4_Inv[7] = V1dotV2_c(P2a,P2_2,2)*mag2i*mag2i;  // r^4

    // Last simple invariant is scalar defined in the conversion to spherical tensors
    P4_Inv[8] = P0[2];

    // Now the orientations between the primary 44 matrix and 42 matrix using rank 2 matrices
    // make contractions of these and normalize
    P4_Inv[9] = V1dotV2_c(P4_2, P2a, 2)*sqrt(norm4*norm42)/(norm2*norm42);// r^12 * r^4 / (r^8 * r^4)
    P4_Inv[10] = V1dotV2_c(P4_22, P2a, 2)*sqrt(norm4*norm42)/(norm2*norm2*norm42);// r^20 * r^4 / (r^16 * r^4)
    P4_Inv[11] = V1dotV2_c(P4_2, P2_2, 2)*sqrt(norm4*norm42)/(norm2*norm42*norm42);// r^16 * r^4 / (r^8 * r^8)

    // And the orientations vs the X4 unit tensors: based on absolute distortion
    O4_Inv[0] = srt3*V1dotV2_c(P4a, Z4, 4)*mag4i;
    O4_Inv[1] = srt3*V1dotV2_c(P4a, X4, 4)*mag4i;
    O4_Inv[2] = srt3*V1dotV2_c(P4a, Y4, 4)*mag4i;
}

//***************************************************************************************************

