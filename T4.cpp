// Version 5.4 07/25/25

//***************************************************************************************************
// 080918: Using standard definition of trace, and added radial term to P2_2_ST
// 092523: Normalization factors made consistent with Colin Adam's thesis
// 072525: Cleaning up for release

#include "BoilerPlate.h"
#include "Constants.h"

//***************************************************************************************************

void P1_2_ST(double *P_1, complex<double> *T1){
	// convert position into complex ang mom vector
	// just creating positive half, get negatives by complex conjugate
	// Should scale as Xlm's for radial dependence: take N(X00) = 1
	// N(Y1m) = sqrt(3); N(R11) = sqrt(2/3)  ==> sqrt(2)
	
	T1[0] = P_1[2]*sqrt2; // z component
	T1[1] = -(P_1[0] + I*P_1[1]); // +1 = -(x + i*y)/sqrt(2)
}

//***************************************************************************************************

void P2_2_ST(double *P_2, complex<double> *P2, double *P0){
	// takes general symmetric 3x3 position matrix,
	// converts them into order 2 and 0 component, respectively, spherical tensors (complex vectors)
	// assumes that matrix is passed as a vector, so index accordingly
	// return with spherical harmonic normalization and phase, defined using Jerphagnon's approach
	// Should scale as Xlm's for radial dependence: take N(X00) = 1
	// N(Y2m) = sqrt(15/2); N(R22) = sqrt(4/15); N(Y2m)*N(R22) = sqrt(2)
	// N(Y00) = 1; N(R20) = sqrt(2/3); N(Y00)*N(X20) = sqrt(2/3)
	
	P0[1] = (P_2[0] + P_2[4] + P_2[8])*sqrt2o3; // J*i*i 
	// redefine symmetric matrix as angular momentum vector: J
	P2[0] =   (P_2[8] + P_2[8] - P_2[0] - P_2[4])*sqrt3i;           // 2*zz - xx - yy
	P2[1] = -((P_2[6] + P_2[2]) + I*(P_2[7] + P_2[5]))*sqrt2i; // -(xz + zx) - i*(yz + zy)
	P2[2] =  ((P_2[0] - P_2[4]) + I*(P_2[3] + P_2[1]))*sqrt2i; //  (xx - yy) + i*(xy + yx)
}

//***************************************************************************************************

void P3_2_ST(double *T_3, complex<double> *P3a, complex<double> *P1a){
	// takes 3x3x3 symmetric tensor, for positions or 3rd derivative of scalar field,
    // all three indices can be interchanged (rst)
	// converts them into one 1, and one 3 rank spherical tensors
	// assumes that matrix is passed as a vector
	// Use spherical harmonic normalization factors and phase
	// Should scale as Xlm's for radial dependence: take N(X00) = 1
	// N(Y3m) = sqrt(5*7/2); N(R33) = sqrt(8/3*5*7); N(Y3m)*N(R33) = sqrt(4/3)
	// N(Y1m) = sqrt(3); N(R31) = sqrt(4/15); N(Y1m)*N(R31) = sqrt(4/5)
	
	// define the vector trace
	// (0,1): R1m*r2  N(Y1m) = sqrt(3); N(R31) = sqrt(4/15); N(Y1m)*N(R20) = sqrt(4/5)
	P1a[0] =  (T_3[18] + T_3[22] + T_3[26])*sqrt4o5;             //  -(zxx + zyy + zzz) = -zrr
	P1a[1] = -(T_3[0]  + T_3[4]  + T_3[8]                       //   (xxx + xyy + xzz) =  xrr
			   + I*(T_3[9]  + T_3[13] + T_3[17]))*sqrt2o5;            //  i(yxx + yyy + yzz) = iyrr
 
	// define traceless symmetric matrix components
	// (2,3) R3m: N(Y3m) = sqrt(35/2); N(R33) = sqrt(8/3*5*7); N(Y3m)*N(R33) = sqrt(4/3)
	P3a[0] =  (T_3[26] + T_3[26] - T_3[2]  - T_3[6]  - T_3[14] - T_3[16] - T_3[18] - T_3[22])*sqrt2o15;
	// 2*zzz - xxz - xzx - yyz - yzy - zxx - zyy
	P3a[1] = ((T_3[4]  + T_3[10] + T_3[12] + T_3[0] + T_3[0] + T_3[0]) //    xyy + yxy + yyx + 3*xxx
			  - four*(T_3[8]  + T_3[20] + T_3[24])                    //  -(xzz + zxz + zzx)
			  + I*((T_3[1]  + T_3[3]  + T_3[9] + T_3[13] + T_3[13] + T_3[13]) //  i(xxy + xyx + yxx + 3*yyy)
				   - four*(T_3[17] + T_3[23] + T_3[25])))*athird*sqrt10i;                  // -i(yzz + zyz + zzy)
	P3a[2] = (T_3[2]  + T_3[6]  + T_3[18] - T_3[14] - T_3[16] - T_3[22]           // xxz + xzx + zxx - (yyz + yzy + zzy)
         + I*(T_3[5]  + T_3[7]  + T_3[11] + T_3[15] + T_3[19] + T_3[21]))*athird; // xyz + xzy + yxz + yzx + zxy + zyx
	P3a[3] = (T_3[4]  + T_3[10] + T_3[12] - T_3[0]                     //    xyy + yxy + yyx - xxx
         - I*(T_3[1]  + T_3[3]  + T_3[9]  - T_3[13]))*sqrt6i;          // -i(xxy + xyx + yxx - yyy)
}

//***************************************************************************************************

void P4_2_ST(double *T_4, complex<double> *P4a, complex<double> *P2a, double *P0){
	// takes symmetric 3x3x3x3 matrix, separates into spherical tensor components
	// which consists of scalar and 2nd and 4th order traceless tensors
	// assumes that matrix is passed as a vector
	// Should scale as Xlm's for radial dependence: take N(X00) = 1

	// define the trace in the conventional sense: R40
	// R40: N(Y00) = 1; N(R40) = sqrt(2/15); N(Y00)*N(R40) = sqrt(2/15)
	P0[2] = (T_4[0] + T_4[4] + T_4[8] + T_4[36] + T_4[40] + T_4[44] + T_4[72] + T_4[76] + T_4[80])*sqrt2o15;
          // xxxx + xxyy + xxzz + yyxx + yyyy + yyzz + zzxx + zzyy + zzzz
	// Alternative scalar
	// (2,1,0): V00*r2 : N(Y00) = 1; N(R40) = sqrt(2/15); N(Y00)*N(R40) = sqrt(2/15)
	//P0[2] = ((T_4[0]  + T_4[0]  - T_4[4]  - T_4[8]                                // 2xxxx - xxyy - xxzz
			  //+ T_4[40] + T_4[40] - T_4[36] - T_4[44]                               // 2yyyy - yyxx - yyzz
			  //+ T_4[80] + T_4[80] - T_4[72] - T_4[76])*sqrt45i                      // 2zzzz - zzxx - zzyy
			 //+ (T_4[10] + T_4[12] + T_4[20] + T_4[24] + T_4[28] + T_4[30]           //  xyxy + xyyx + xzxz + xzzx + yxxy + yxyx
				//+  T_4[50] + T_4[52] + T_4[56] + T_4[60] + T_4[68] + T_4[70])*sqrt20i*sqrt2o15); //  yzyz + yzzy + zxxz + zxzx + zyyz + zyzy
	
	// no vector components
	// define 2nd order components: R42m
	// R42m: N(Y2m) = sqrt(15/2); N(R42) = sqrt(8/(3*5*7)); N(Y2m)*N(R42) = sqrt(4/7)
	P2a[0] = -((T_4[0]  + T_4[4]  + T_4[8]  + T_4[36] + T_4[40] + T_4[44])   //    xxxx + xxyy + xxzz + yyxx + yyyy + yyzz
			 - two*(T_4[72] + T_4[76] + T_4[80]))*sqrt2o21;                      // -2(zzxx + zzyy + zzzz)
	P2a[1] = -((T_4[18] + T_4[22] + T_4[26] + T_4[54] + T_4[58] + T_4[62]            //    xzxx + xzyy + xzzz + zxxx + zxyy + zxzz
		   + I*(T_4[45] + T_4[49] + T_4[53] + T_4[63] + T_4[67] + T_4[71])))*sqrt7i; //  i(yzxx + yzyy + yzzz + zyxx + zyyy + zyzz)
	P2a[2] =-((-T_4[0]  - T_4[4]  - T_4[8]  + T_4[40] + T_4[36] + T_4[44])   //    yyxx - xxxx + yyyy - xxyy + yyzz - xxzz
		  + I*(-T_4[9]  - T_4[13] - T_4[17] - T_4[27] - T_4[31] - T_4[35]))*sqrt7i;  // -i(xyxx + xyyy + xyzz + yxxx + yxyy + yxzz)
	// Alternative matrix
	// (2,3,2) V2m *i*i : N(Y2m) = sqrt(15/2); N(R42) = sqrt(8/(3*5*7)); N(Y2m)*N(R42) = sqrt(4/7)
/*	P2a[0] = (-((T_4[0] + T_4[40] + T_4[56] + T_4[60] + T_4[68] + T_4[70] + T_4[72] + T_4[76])*sqrt3o70
				// (xxxx + yyyy + zxxz + zxzx + zyyz + zyzy + zzxx + zzyy)
				+ (T_4[4] + T_4[10] + T_4[12] + T_4[28] + T_4[30] + T_4[36])*sqrt210i  // +(xxyy + xyxy + xyyx + yxxy + yxyx + yyxx)
				- (T_4[8] + T_4[20] + T_4[24] + T_4[44] + T_4[50] + T_4[52])*sqrt8o105 // -(xxzz + xzxz + xzzx + yyzz + yzyz + yzzy)
				- T_4[80]*sqrt6o35))*sqrt4o7;                                                   // -zzzz
	P2a[1] = (-((T_4[32] + T_4[34] + T_4[38] + T_4[42] + T_4[46] + T_4[48])*sqrt5o252  //  (yxyz + yxzy + yyxz + yyzx + yzxy + yzyx)
				+ (T_4[2]  + T_4[6]  + T_4[18] + T_4[62] + T_4[74] + T_4[78])*sqrt16o315 // +(xxxz + xxzx + xzxx + zxzz + zzxz + zzzx)
				- (T_4[14] + T_4[16] + T_4[22] + T_4[58] + T_4[64] + T_4[66])*sqrt315i   // -(xyyz + xyzy + xzyy + zxyy + zyxy + zyyx)
				- (T_4[26] + T_4[54])*sqrt35i)                                            // -(xzzz + zxxx)
			  - I*((T_4[5]  + T_4[7]  + T_4[11] + T_4[15] + T_4[19] + T_4[21])*sqrt5o252  // i(xxyz + xxzy + xyxz + xyzx + xzxy + xzyx)
				   + (T_4[41] + T_4[43] + T_4[49] + T_4[71] + T_4[77] + T_4[79])*sqrt16o315 //+i(yyyz + yyzy + yzyy + zyzz + zzyz + zzzy)
				   - (T_4[29] + T_4[33] + T_4[45] + T_4[55] + T_4[57] + T_4[63])*sqrt315i   //-i(yxxz + yxzx + yzxx + zxxy + zxyx + zyxx)
				   - (T_4[53] + T_4[67])*sqrt35i))*sqrt4o7;                                          //-i(yzzz + zyyy)
	P2a[2] = (-((T_4[8]  + T_4[20] + T_4[24] - T_4[44] - T_4[50] - T_4[52])*sqrt315i   //  (xxzz + xzxz + xzzx - yyzz - yzyz - yzzy)
				+ (-T_4[56] - T_4[60] + T_4[68] + T_4[70] - T_4[72] + T_4[76])*sqrt5o252  // (-zxxz - zxzx + zyyz + zyzy - zzxx + zzyy)
				+ ( T_4[4]  + T_4[10] + T_4[12] - T_4[28] - T_4[30] - T_4[36])*sqrt7o180  //  (xxyy + xyxy + xyyx - yxxy - yxyx - yyxx)
				+ (-T_4[0]  + T_4[40])*sqrt9o140)                                          //  (yyyy - xxxx)
			  - I*((T_4[17] + T_4[23] + T_4[25] + T_4[35] + T_4[47] + T_4[51])*sqrt315i   // i(xyzz + xzyz + xzzy + yxzz + yzxz + yzzx)
				   -   ( T_4[59] + T_4[61] + T_4[65] + T_4[69] + T_4[73] + T_4[75])*sqrt5o252  //-i(zxyz + zxzy + zyxz + zyzx + zzxy + zzyx)
				   -   ( T_4[1]  + T_4[3]  + T_4[9]  + T_4[31] + T_4[37] + T_4[39])*sqrt16o315 //-i(xxxy + xxyx + xyxx + yxyy + yyxy + yyyx)
				   +   ( T_4[13] + T_4[27])*sqrt35i))*sqrt4o7;*/
	// no 3rd order components
	// and then the 4th order component: R44m
	// R44m: N(Y4m) = sqrt(5*7*9/8); N(R44) = sqrt(16/(3*5*7*9)); N(Y4m)*N(R44) = sqrt(2/3)
	P4a[0] =  ((T_4[4]  + T_4[10] + T_4[12] + T_4[28] + T_4[30] + T_4[36])    //   xxyy + xyxy + xyyx + yxxy + yxyx + yyxx
       + three*(T_4[0]  + T_4[40])                                            // 3(xxxx + yyyy)
       + four*(-T_4[8]  - T_4[20] - T_4[24] - T_4[44] - T_4[50] - T_4[52]     // -(xxzz + xzxz + xzzx + yyzz + yzyz + yzzy)
			  - T_4[56] - T_4[60] - T_4[68] - T_4[70] - T_4[72] - T_4[76])    // -(zxxz + zxzx + zyyz + zyzy + zzxx + zzyy)
        + eight*T_4[80])*sqrt420i;                                            //  4zzzz
	P4a[1] =  ((T_4[14] + T_4[16] + T_4[22] + T_4[32] + T_4[34] + T_4[38]     //   xyyz + xyzy + xzyy + yxyz + yxzy + yyxz
			  + T_4[42] + T_4[46] + T_4[48] + T_4[58] + T_4[64] + T_4[66])    //   yyzx + yzxy + yzyx + zxyy + zyxy + zyxy
			 + three*(T_4[2]  + T_4[6]  + T_4[18] + T_4[54])                  // 3(xxxz + xxzx + xzxx + zxxx)
		    + four*(-T_4[26] - T_4[62] - T_4[74] - T_4[78])                   // -(xzzz + zxzz + zzxz + zzzx)
		  + I*((T_4[5]  + T_4[7]  + T_4[11] + T_4[15] + T_4[19] + T_4[21]     // i(xxyz + xxzy + xyxz + xyzx + xzxy + xzyx)
		 	  + T_4[29] + T_4[33] + T_4[45] + T_4[55] + T_4[57] + T_4[63])    // i(yxxz + yxzx + yzxx + zxxy + zxyx + zyxx)
			 + three*(T_4[41] + T_4[43] + T_4[49] + T_4[67])                  // i(yyyz + yyzy + yzyy + zyyy)
		    + four*(-T_4[53] - T_4[71] - T_4[77] - T_4[79])))*sqrt336i;       //-i(yzzz + zyzz + zzyz + zzzy)
	P4a[2] = (two*(-T_4[0]  + T_4[8]  + T_4[20] + T_4[24] + T_4[40] - T_4[44] //  -xxxx + xxzz + xzxz + xzzx + yyyy - yyzz
			  - T_4[50] - T_4[52] + T_4[56] + T_4[60] - T_4[68] - T_4[70]     //  -yzyz - yzzy + zxxz + zxzx - zyyz - zyzy
			  + T_4[72] - T_4[76])                                            //   zzxx - zzyy
		 + I*((-T_4[1]  - T_4[3]  - T_4[9]  - T_4[13]                         //-i(xxxy + xxyx + xyxx + xyyy)
			  - T_4[27] - T_4[31] - T_4[37] - T_4[39])                        //-i(yxxx + yxyy + yyxy + yyyx)
			 + two*(T_4[17] + T_4[23] + T_4[25] + T_4[35] + T_4[47] + T_4[51] // i(xyzz + xzyz + xzzy + yxzz + yzxz + yzzx)
			  + T_4[59] + T_4[61] + T_4[65] + T_4[69] + T_4[73] + T_4[75])))*sqrt42i;// i(zxyz + zxzy + zyxz + zyzx + zzxy + zzyx)
	P4a[3] = ((-T_4[2]  - T_4[6]  + T_4[14] + T_4[16] - T_4[18] + T_4[22]   //  -xxxz - xxzx + xyyz + xyzy - xzxx + xzyy
			  + T_4[32] + T_4[34] + T_4[38] + T_4[42] + T_4[46] + T_4[48]   //   yxyz + yxzy + yyxz + yyzx + yzxy + yzyx
			  - T_4[54] + T_4[58] + T_4[64] + T_4[66])                      //  -zxxx + zxyy + zyxy + zyyx
		  + I*(-T_4[5]  - T_4[7]  - T_4[11] - T_4[15] - T_4[19] - T_4[21]   //-i(xxyz + xxzy + xyxz + xyzx + xzxy + xzyx)
			  - T_4[29] - T_4[33] + T_4[41] + T_4[43] - T_4[45] + T_4[49]   //-i(yxxz + yxzx - yyyz - yyzy + yzxx - yzyy)
			  - T_4[55] - T_4[57] - T_4[63] + T_4[67]))*sqrt48i;            //-i(zxxy + zxyx + zyxx - zyyy)
	P4a[4] = (( T_4[0]  - T_4[4]  - T_4[10] - T_4[12] - T_4[28] - T_4[30] - T_4[36] + T_4[40])
                                                                            //   xxxx - xxyy - xyxy - xyyx - yxxy - yxyx - yyxx + yyyy
		   + I*(T_4[1]  + T_4[3]  + T_4[9]  - T_4[13] + T_4[27] - T_4[31] - T_4[37] - T_4[39]))*sqrt24i;
                                                                            // i(xxxy + xxyx + xyxx - xyyy + yxxx - yxyy - yyxy - yyyx)
}


