// Version 5.0 02/21/21
// Version 5.4 11/10/25 Clean up for release

// Global dims
#define NPATOM 12 // coordinates per atom: position, velocity, force,
#define NAID 5    // how many ID flags per atom
#define LMAX 100  // max number of atoms per cell; tied to radius

// Functions
void Factorial();
double CG_calc(int j1, int m1, int j2, int m2, int j, int m);
bool CG_check(int j1, int m1, int j2, int m2, int j, int m);
double Phase(int n);
double Sum_term(int j1, int m1, int j2, int m2, int j, int m);
void ClebschGordan();
double Mag(double vect[], int length);
double Mag_c(complex<double> vect[], int length);
double Trace(double vect[], int length);
void P1_2_ST(double P_1[], complex<double> T1[]);
double V1dotV2(double vect1[], double vect2[], int length);
double V1dotV2_c(complex<double> vect1[], complex<double> vect2[], int length);
void P1_Invar(complex<double> P1[], double P1_Inv[],
			  complex<double> Z[], complex<double> X[], complex<double> Y[], double O1_Inv[]);
void P2_2_ST(double P_2[], complex<double> P2[], double P0[]);
void P2_Invar(complex<double> P2[], double P0[], double P2_Inv[],
			  complex<double> Z2[], complex<double> Y2[], complex<double> X2[], double O2_Inv[]);
void P3_2_ST(double P_3[], complex<double> P3a[], complex<double> P1[]);

void P3_Invar(complex<double> P3a[], complex<double> P1a[], double P3_Inv[], complex<double> P_1[],
	complex<double> Z3[], complex<double> X3[], complex<double> Y3[], double O3_Inv[]);

void P4_2_ST(double T_4[], complex<double> P4a[], complex<double> P2a[], double P0[]);

void P4_Invar(complex<double> P4a[], complex<double> P2a[], double P0[], complex<double> P2[], double P4_Inv[],
			  complex<double> Z4[], complex<double> X4[], complex<double> Y4[], double O4_Inv[]);

void V1outV2_c(complex<double> Vect1[], int j1, complex<double> Vect2[], int j2, complex<double> Vect[], int j);

// Other Global Variables

char ParFileName[20] = "params.txt";
char XYZFileName[128] = "fin.txt";
char BINFileName[20] = "binfile";
char BoxFileName[20] = "binfile.parm";
char MomFileName[20] = "Moment.txt";
char VelFileName[20] = "Velocities.txt";
char ForFileName[20] = "Forces.txt";
char CorFileName[20] = "Correlate.txt";
char AtomFileName[20] = "Atoms.txt";
char CellFileName[20] = "Cells.txt";
char DistFileName[20] = "Dists.txt";

bool XPER, YPER, ZPER;  // Logicals for boundary conditions
char xdum[3],ydum[3],zdum[3]; //Strings for boundary conditions
int timestep;// Store the timestep
bool BINYES;            //Logical for reading binary file
bool REFWEIGHT;         //Logical for weighting using reference configuration
double XLO, XHI, YLO, YHI, ZLO, ZHI;  // box boundaries
double RXLO, RXHI, RYLO, RYHI, RZLO, RZHI;  // box boundaries
double XLEN, YLEN, ZLEN;  // box lengths
double RXLEN, RYLEN, RZLEN;  // box lengths
double radius;  // radius of neighborhood
double sigma; // Gaussian coefficient
double XDIM, YDIM, ZDIM;
int XCELLS, YCELLS, ZCELLS;
double RXDIM, RYDIM, RZDIM;
int RXCELLS, RYCELLS, RZCELLS;
int DSTYLE;

/* variables controlling the DMA of arrays, assigned in main, GetCoords */
int natom;
int ncells;
int rncells;

//for reading SPaSM binary files
#define BUFSIZE       1024   /* Number of atoms to read at a time */

int dx[13] = { 1, 1, 0,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};
int dy[13] = { 0, 1, 1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
int dz[13] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
