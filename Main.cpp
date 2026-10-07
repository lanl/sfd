//***************************************************************************************************
//
// Code for the calculation of Strain Functional Descriptors (SFD)
// Version 5.4 09/25/25
// The mathematical basis for the SFD is described in:
// Kober, E. M., Tavenner, J. P., Adams, C. M. & Mathew, N.
// "Strain functionals: A complete and symmetry-adapted set of descriptors to
// characterize atomistic configurations" https://arxiv.org/abs/2402.04191v1 (2024).
// (Revisions and full publication pending)
//
// This version calculates the SFD through 4th order for atomic positions
// The input files are primarily expected to be in a LAMMPS dump format (text file)
// The capability to read SPaSM binary files is also included
//
// Original code was written in c by E. M. Kober
// Code modifications and dynamic memory allocation by Nithin Mathew

// Within this code, 'the reference', from which this entire approach is taken, is
// "3D Moment Forms: Their Construction and Application to Object Identification
// and Positioning"
// by Chong-Huah Lo and Hon-Son Don
// from IEEE Transactions on Pattern Analysis and Machine Intelligence,
// Vol 2, No. 10, October 1989

// REVISIONS:

// Main Function
// Moments code to read straight from LAMMPS file
// Added DSTYLE index to params.txt to identify dumpstyle
// 072525: cleaning up for release version

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <sys/stat.h>
#include <float.h>
#include <cmath>
#include <complex>

#ifdef OMP
#include <omp.h>
#endif

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h> // for read and write functions

#include "NFSmoment.h"
#include "MainData.h"

/* Sav Data structure */
// typedef struct { double x; double y; double z;} Vector;
// typedef struct { int type; Vector r; Vector s;
// } SaveParticle;


//***************************************************************************************************
int getSavFileNAtoms ()
{
    
// For reading a SPaSM binary file
    
  struct stat StatResults;
  int FileSize;

  if(stat(BINFileName, &StatResults) == 0)
    FileSize = StatResults.st_size;
  else{
    cout << "Error: Unable to access atomic coordinate file "
         << BINFileName << endl;
    return 2;
  }
  if(FileSize % (3*sizeof(SaveParticle)) != 0){
    cout << "Warning: Input file " << BINFileName
         << " does not appear to contain an integral number of atoms."
         << endl;
    return 1;
  }
  natom = FileSize / (sizeof(SaveParticle));
  cout << "Estimated number of atoms = " << natom << endl;
	return 0;
}

//***************************************************************************************************
void readSavFile(double **tmpCoords, int **tmpAtomID)
{
    
// For reading a SPaSM binary file
    
    SaveParticle *pt;
    SaveParticle *databuf;

    int savFileDescriptor;
    int i, n, np, fmax;
    int xshift, xtrunc, yshift, ytrunc, zshift, ztrunc;
  	double tx, ty, tz;
  	double tvx, tvy, tvz;
	double abit, big;

    abit = 1.e-7;
    big  = 1.e7;

  // number of particles read.
  	np = 0;

  // allocate memory to temporry buffer that is used to read sav file.
  	fmax = BUFSIZE * sizeof(SaveParticle);
  	databuf = (SaveParticle *) malloc(fmax);

  // open sav file.
  	cout << "Reading sav file " << BINFileName << endl;
  	if ((savFileDescriptor = open(BINFileName,O_RDONLY)) < 0) {
    		cout << "error, can't open Sav file " << BINFileName << endl;
    		exit (1);
  	}

  	do {
    		n=read(savFileDescriptor, databuf,fmax);

    		if (n>0)
      			n /=sizeof(SaveParticle);
    		else if (n < 0) {
      			cout << "error reading file, is SaveParticle correct?" << endl;
      			exit(1);
    		}

    		pt = databuf;

    		for (i = 0, pt = databuf; i < n; i++, pt++) {

       			tmpAtomID[np][0] = np;

                // for now make all the atoms have type 0
      			tmpAtomID[np][1] = 0;

      
      			tmpCoords[np][0] = pt->r.x;
      			tmpCoords[np][1] = pt->r.y;
      			tmpCoords[np][2] = pt->r.z;
      			
      			np++;
		}
	} while (n == BUFSIZE);

  	cout << "np read from Sav file = " << np <<endl;
	cout << "box lengths from binfile.parm file" << endl;

  	cout << "XLO = " << XLO << " XHI = " << XHI << endl;
  	cout << "YLO = " << YLO << " YHI = " << YHI << endl;
  	cout << "ZLO = " << ZLO << " ZHI = " << ZHI << endl;
    //      reset to 0 - max
  	for(i = 0; i < natom; i++){
  		tmpCoords[i][0] -= XLO;
		tmpCoords[i][1] -= YLO;
		tmpCoords[i][2] -= ZLO;
	}

        XHI -= XLO;
        YHI -= YLO;
        ZHI -= ZLO;
        XLO = 0.;
        YLO = 0.;
        ZLO = 0.;

        XPER = true;
        YPER = false;
        ZPER = true;

  	XLEN = abs(XHI - XLO);
  	YLEN = abs(YHI - YLO);
  	ZLEN = abs(ZHI - ZLO);
        XCELLS  = max(2,int(0.7*XLEN/radius));
        YCELLS  = max(2,int(0.7*YLEN/radius));
        ZCELLS  = max(2,int(0.7*ZLEN/radius));
        ncells   = XCELLS*YCELLS*ZCELLS;
        XDIM    = XLEN/XCELLS;
        YDIM    = YLEN/YCELLS;
        ZDIM    = ZLEN/ZCELLS;
        cout << "XCELLS = " << XCELLS << " YCELLS " << YCELLS << " ZCELLS = " << ZCELLS << endl;
        cout << "NCELLS = " << ncells << endl;
	
	cout << "Done with binary file" << endl;
  	close(savFileDescriptor);

	xshift = 0;
        xtrunc = 0;
        yshift = 0;
        ytrunc = 0;
        zshift = 0;
        ztrunc = 0;
        for(i = 0; i < natom; i++){

                if(tmpCoords[i][0] < XLO) {
                        if (XPER) {
                                tmpCoords[i][0] += XLEN;
                                xshift += 1;
                        }
                        else {
                                tmpCoords[i][0] = XLO;
                                xtrunc += 1;
                        }
                }
                if(tmpCoords[i][0] >= XHI) {
                        if (XPER) {
                                tmpCoords[i][0] -= XLEN;
                                xshift += 1;
                        }
                        else {
                                tmpCoords[i][0] = XHI - abit;
                                xtrunc += 1;
                        }
                }
               if(tmpCoords[i][1] < YLO) {
                        if (YPER) {
                                tmpCoords[i][1] += YLEN;
                                yshift += 1;
                        }
                        else {
                                tmpCoords[i][1] = YLO;
                                ytrunc += 1;
                        }
                }
                if(tmpCoords[i][1] >= YHI) {
                        if (YPER) {
                                tmpCoords[i][1] -= YLEN;
                                yshift += 1;
                        }
                        else {
                                tmpCoords[i][1] = YHI - abit;
                                ytrunc += 1;
                        }
                }
                if(tmpCoords[i][2] < ZLO) {
                        if (ZPER) {
                                tmpCoords[i][2] += ZLEN;
                                zshift += 1;
                        }

                        else {
                                tmpCoords[i][2] = ZLO;
                                ztrunc += 1;
                        }
                }
                if(tmpCoords[i][2] >= ZHI) {
                        if (ZPER) {
                                tmpCoords[i][2] -= ZLEN;
                                zshift += 1;
                        }
                        else {
                                tmpCoords[i][2] = ZHI - abit;
                                ztrunc += 1;
                        }
                }
        }

        cout << "Read in boundary conditions and shifted origin " << endl;
        cout << "XPER = " << (int)XPER << " XLO= " << XLO << " XHI = " << XHI << endl;
        cout << "YPER = " << (int)YPER << " YLO= " << YLO << " YHI = " << YHI << endl;
        cout << "ZPER = " << (int)ZPER << " ZLO= " << ZLO << " ZHI = " << ZHI << endl;
        cout << xshift << " atoms were read in that needed to be shifted for x periodicity " << endl;
        cout << xtrunc << " atoms were read in that needed to be truncated at x boundary " << endl;
        cout << yshift << " atoms were read in that needed to be shifted for y periodicity " << endl;
        cout << ytrunc << " atoms were read in that needed to be truncated at y boundary " << endl;
        cout << zshift << " atoms were read in that needed to be shifted for z periodicity " << endl;
        cout << ztrunc << " atoms were read in that needed to be truncated at z boundary " << endl;
}


//***************************************************************************************************
void GetCoords(double **Coords, double **RefCoords, int **AtomID, int **RefAtomID){
	/*
     For reading LAMMPS dump/data files in txt format; several options included
     Current version of the code can handle periodic boundaries for orthorhombic cells
     Otherwise an oversized non-periodic orthorhombic box will be defined to contain the atoms
	 Fills the array Coords with atomic coordinates loaded from the input file
	 InputFileName, assuming NumbersPerAtom double-precision floating point
	 numbers per atom, in the format
	 Assigns each atom to a local box, and saves those indices to AtomID
	 */
	double abit, big;
	int i, per, dummy;
	int xshift, xtrunc, yshift, ytrunc, zshift, ztrunc;
	char label;
	bool test;
	int junknum;
	char junk;
	
	FILE *ifile, *reffile;
	ifile = fopen(XYZFileName,"r");
	reffile = fopen(XYZFileName, "r");
	
	abit = 1.e-7;
    big  = 1.e7;
	
	cout << "Getting Coords" << endl;
	// Load atomic coordinates
    // This depends on the LAMMPS dump styles that people have used
    // Feel free to write your own
    
    if(DSTYLE==3){  // Jacob's LAMMPS style
		cout << "DSTYLE = " << DSTYLE << endl;
        fscanf(ifile,"%s %s\n", &junk,&junk);
        fscanf(ifile, "%d\n",&dummy);
        fscanf(ifile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(ifile, "%d\n",&dummy);
        fscanf(ifile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(ifile,"%lf %lf\n",&XLO,&XHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&YLO,&YHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&ZLO,&ZHI);  // xlo
        fscanf(ifile,"%s %s %s %s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%s %s\n", &junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%lf %lf\n",&RXLO,&RXHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RYLO,&RYHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RZLO,&RZHI);  // xlo
        fscanf(reffile,"%s %s %s %s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
    }

    if(DSTYLE==2){  // another LAMMPS style
        fscanf(ifile,"%s %s\n", &junk,&junk);
        fscanf(ifile, "%d\n",&dummy);
        fscanf(ifile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(ifile, "%d\n",&dummy);
        fscanf(ifile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(ifile,"%lf %lf\n",&XLO,&XHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&YLO,&YHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&ZLO,&ZHI);  // xlo
        fscanf(ifile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%s %s\n", &junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%lf %lf\n",&RXLO,&RXHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RYLO,&RYHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RZLO,&RZHI);  // xlo
        fscanf(reffile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
    }

    if(DSTYLE==1){  // another LAMMPS style
        fscanf(ifile,"%s %s\n", &junk,&junk);
        fscanf(ifile, "%d\n",&timestep);
        fscanf(ifile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(ifile, "%d\n",&dummy);
        fscanf(ifile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,xdum,ydum,zdum);
        fscanf(ifile,"%lf %lf\n",&XLO,&XHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&YLO,&YHI);  // xlo
        fscanf(ifile,"%lf %lf\n",&ZLO,&ZHI);  // xlo
        fscanf(ifile,"%s %s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%s %s\n", &junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s\n", &junk,&junk,&junk,&junk);
        fscanf(reffile, "%d\n",&dummy);
        fscanf(reffile,"%s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk);
        fscanf(reffile,"%lf %lf\n",&RXLO,&RXHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RYLO,&RYHI);  // xlo
        fscanf(reffile,"%lf %lf\n",&RZLO,&RZHI);  // xlo
        fscanf(reffile,"%s %s %s %s %s %s %s\n", &junk,&junk,&junk,&junk,&junk,&junk,&junk);
    }

    if(DSTYLE==0){
	// another LAMMPS style
        fscanf(ifile,"%s %s %s\n", &junk,&junk,&junk);
        fscanf(ifile,"%d\n",&dummy);
        fscanf(ifile,"%lf %lf %s %s\n",&XLO,&XHI,&junk,&junk);  // xlo
        fscanf(ifile,"%lf %lf %s %s\n",&YLO,&YHI,&junk,&junk);  // xlo
        fscanf(ifile,"%lf %lf %s %s\n",&ZLO,&ZHI,&junk,&junk);  // xlo
        cout << "XLO = " << XLO << " XHI = " << XHI << endl;
        cout << "YLO = " << YLO << " YHI = " << YHI << endl;
        cout << "ZLO = " << ZLO << " ZHI = " << ZHI << endl;
        fscanf(ifile,"%d %s\n",&dummy,&junk);  // xlo
        cout << "natom = " << dummy << " label = " << junk << endl;
        fscanf(ifile,"%d %s %s\n",&dummy,&junk,&junk);  // xlo
        cout << "atomtypes = " << dummy << " label = " << junk << endl;
        fscanf(ifile,"%s\n",&junk);
        cout << "label = " << junk << endl;
        fscanf(reffile,"%s %s %s\n", &junk,&junk,&junk);
        fscanf(reffile,"%d\n",&dummy);
        fscanf(reffile,"%lf %lf %s %s\n",&RXLO,&RXHI,&junk,&junk);  // xlo
        fscanf(reffile,"%lf %lf %s %s\n",&RYLO,&RYHI,&junk,&junk);  // xlo
        fscanf(reffile,"%lf %lf %s %s\n",&RZLO,&RZHI,&junk,&junk);  // xlo
        fscanf(reffile,"%d %s\n",&dummy,&junk);  // xlo
        fscanf(reffile,"%d %s %s\n",&dummy,&junk,&junk);  // xlo
        fscanf(reffile,"%s\n",&junk);
    }
    
    if((DSTYLE!=0)&&(DSTYLE!=1)&&(DSTYLE!=2)&&(DSTYLE!=3)){
        cout << "Oops, dump style not defined" << endl;
    }

    if(DSTYLE==3){
        for(i = 0; i < natom; i++){

            fscanf(ifile,"%d %d %lf %lf %lf %lf %lf\n",
				   &AtomID[i][0],&AtomID[i][1],
				   &Coords[i][0],&Coords[i][1],&Coords[i][2],
					&Coords[i][3],&Coords[i][4]);
             fscanf(reffile,"%d %d %lf %lf %lf %lf %lf\n",
				   &AtomID[i][0],&AtomID[i][1],
				   &Coords[i][0],&Coords[i][1],&Coords[i][2],
				&Coords[i][3],&Coords[i][4]);
        }
    }

    if(DSTYLE==2){
        for(i = 0; i < natom; i++){

            fscanf(ifile,"%d %lf %lf %lf\n",
                    &AtomID[i][0],
                    &Coords[i][0],&Coords[i][1],&Coords[i][2]);
            fscanf(reffile,"%d %lf %lf %lf\n",
                    &RefAtomID[i][0],
                    &RefCoords[i][0],&RefCoords[i][1],&RefCoords[i][2]);
        }
    }

    if(DSTYLE==1){
        for(i = 0; i < natom; i++){

            fscanf(ifile,"%d %d %lf %lf %lf\n",
					&AtomID[i][0],&AtomID[i][1],
					&Coords[i][0],&Coords[i][1],&Coords[i][2]);
            fscanf(reffile,"%d %d %lf %lf %lf\n",
					&RefAtomID[i][0],&RefAtomID[i][1],
					&RefCoords[i][0],&RefCoords[i][1],&RefCoords[i][2]);
        }
    }
    
    if(DSTYLE==0){
        for(i = 0; i < natom; i++){

            fscanf(ifile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					&AtomID[i][0],&AtomID[i][1],
					&Coords[i][0],&Coords[i][1],&Coords[i][2],
					&Coords[i][3],&Coords[i][4],&Coords[i][5],
					&Coords[i][6],&Coords[i][7],&Coords[i][8]);
            fscanf(reffile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					&RefAtomID[i][0],&RefAtomID[i][1],
					&RefCoords[i][0],&RefCoords[i][1],&RefCoords[i][2],
					&RefCoords[i][3],&RefCoords[i][4],&RefCoords[i][5],
					&RefCoords[i][6],&RefCoords[i][7],&RefCoords[i][8]);
        }
    }
	
    // adjustments so atoms don't sit on boundaries
	XHI -= abit;
	YHI -= abit;
	ZHI -= abit;
	cout << "XLO = " << XLO << " XHI = " << XHI << endl;
	cout << "YLO = " << YLO << " YHI = " << YHI << endl;
	cout << "ZLO = " << ZLO << " ZHI = " << ZHI << endl;
	//	reset to 0 - max

	for(i = 0; i < natom; i++){

		Coords[i][0] -= XLO;
		Coords[i][1] -= YLO;
		Coords[i][2] -= ZLO;
		Coords[i][9] = 0;
		Coords[i][10] = 0;
		Coords[i][11] = 0;
		RefCoords[i][0] -= RXLO;
		RefCoords[i][1] -= RYLO;
		RefCoords[i][2] -= RZLO;
		RefCoords[i][9] = 0;
		RefCoords[i][10] = 0;
		RefCoords[i][11] = 0;

	}
	
	XHI -= XLO;
	YHI -= YLO;
	ZHI -= ZLO;
	XLO = 0.;
	YLO = 0.;
	ZLO = 0.;
	RXHI -= RXLO;
	RYHI -= RYLO;
	RZHI -= RZLO;
	RXLO = 0.;
	RYLO = 0.;
	RZLO = 0.;
	XPER = (strcmp(xdum, "pp") == 0);
	YPER = (strcmp(ydum, "pp") == 0);
	ZPER = (strcmp(zdum, "pp") == 0);

	XLEN = XHI - XLO;
	YLEN = YHI - YLO;
	ZLEN = ZHI - ZLO;

	RXLEN = RXHI - RXLO;
	RYLEN = RYHI - RYLO;
	RZLEN = RZHI - RZLO;
	XCELLS  = max(2,int(0.7*XLEN/radius));
	YCELLS  = max(2,int(0.7*YLEN/radius));
	ZCELLS  = max(2,int(0.7*ZLEN/radius));
	ncells   = XCELLS*YCELLS*ZCELLS;
	XDIM    = XLEN/XCELLS;
	YDIM    = YLEN/YCELLS;
	ZDIM    = ZLEN/ZCELLS;
	cout << "XCELLS = " << XCELLS << " YCELLS " << YCELLS << " ZCELLS = " << ZCELLS << endl;
	cout << "NCELLS = " << ncells << endl;
	
	RXCELLS  = max(2,int(0.7*RXLEN/radius));
	RYCELLS  = max(2,int(0.7*RYLEN/radius));
	RZCELLS  = max(2,int(0.7*RZLEN/radius));
	rncells  = RXCELLS*RYCELLS*RZCELLS;
	RXDIM    = RXLEN/RXCELLS;
	RYDIM    = RYLEN/RYCELLS;
	RZDIM    = RZLEN/RZCELLS;

	
	cout << "Done with ASCII file" << endl;
	fclose(ifile);
	fclose(reffile);

	// count the number of outliers; LAMMPS cleans up before timestep, prints out after
	// don't really need to print this out, just checking when it's necessary
	// can't have atoms sitting on upper boundaries because of box numbering algorithm
	
	xshift = 0;
	xtrunc = 0;
	yshift = 0;
	ytrunc = 0;
	zshift = 0;
	ztrunc = 0;
	
	// Get the raw numbers
	for(i = 0; i < natom; i++){
			
		// Check for outliers and shift positions or bounds accordingly
		if(Coords[i][0] < XLO) {
            //cout << Coords[i][0] << endl;
			if (XPER) {
				Coords[i][0] += XLEN;
				Coords[i][9] =-1;
				xshift += 1;
			}
			else {
				Coords[i][0] = XLO;
				xtrunc += 1;
			}
		}else{
			if(Coords[i][0] >= XHI) {
                //cout << Coords[i][0] << endl;
				if (XPER) {
					Coords[i][0] -= XLEN;
					Coords[i][9] =1;
					xshift += 1;
				}
				else {
					Coords[i][0] = XHI - abit;
					xtrunc += 1;
				}
			}else{
				Coords[i][9] =0;
			}
		}
		if(Coords[i][1] < YLO) {
			if (YPER) {
				Coords[i][1] += YLEN;
				Coords[i][10] =-1;
				yshift += 1;
			}
			else {
				Coords[i][1] = YLO;
				ytrunc += 1;
			}
		}else{
			if(Coords[i][1] >= YHI) {
				if (YPER) {
					Coords[i][1] -= YLEN;
					Coords[i][10] =1;
					yshift += 1;
				}
				else {
					Coords[i][1] = YHI - abit;
					ytrunc += 1;
				}
			}else{
				Coords[i][10] =0;
			}
		}
		if(Coords[i][2] < ZLO) {
			if (ZPER) {
				Coords[i][2] += ZLEN;
				Coords[i][11] =-1;
				zshift += 1;
			}
				
			else {
				Coords[i][2] = ZLO;
				ztrunc += 1;
			}
		}else{
			if(Coords[i][2] >= ZHI) {
				if (ZPER) {
					Coords[i][2] -= ZLEN;
					Coords[i][11] =1;
					zshift += 1;
				}
				else {
					Coords[i][2] = ZHI - abit;
					ztrunc += 1;
				}
			}else{
				Coords[i][11] =0;
			}	
		}

	//do it for the reference configuration
		if(RefCoords[i][0] < RXLO) {
			if (XPER) {
				RefCoords[i][0] += RXLEN;
				RefCoords[i][9] =-1;
			}
			else {
				RefCoords[i][0] = RXLO;
			}
		}else{
			if(RefCoords[i][0] >= RXHI) {
				if (XPER) {
					RefCoords[i][0] -= RXLEN;
					RefCoords[i][9] =1;
				}
				else {
					RefCoords[i][0] = RXHI - abit;
				}
			}else{
				RefCoords[i][9] =0;
			}
		}
		if(RefCoords[i][1] < RYLO) {
			if (YPER) {
				RefCoords[i][1] += RYLEN;
				RefCoords[i][10] =-1;
			}
			else {
				RefCoords[i][1] = RYLO;
			}
		}else{
			if(RefCoords[i][1] >= RYHI) {
				if (YPER) {
					RefCoords[i][1] -= RYLEN;
					RefCoords[i][10] =1;
				}
				else {
					RefCoords[i][1] = RYHI - abit;
				}
			}else{
				RefCoords[i][10] =0;
			}
		}
		if(RefCoords[i][2] < RZLO) {
			if (ZPER) {
				RefCoords[i][2] += RZLEN;
				RefCoords[i][11] =-1;
			}
				
			else {
				RefCoords[i][2] = RZLO;
			}
		}else{
			if(RefCoords[i][2] >= RZHI) {
				if (ZPER) {
					RefCoords[i][2] -= RZLEN;
					RefCoords[i][11] =1;
				}
				else {
					RefCoords[i][2] = RZHI - abit;
				}
			}else{
				RefCoords[i][11] =0;
			}	
		}
	}
	
	cout << "Read in boundary conditions and shifted origin " << endl;
	cout << "XPER = " << (int)XPER << " XLO= " << XLO << " XHI = " << XHI << endl;
	cout << "YPER = " << (int)YPER << " YLO= " << YLO << " YHI = " << YHI << endl;
	cout << "ZPER = " << (int)ZPER << " ZLO= " << ZLO << " ZHI = " << ZHI << endl;
	cout << xshift << " atoms were read in that needed to be shifted for x periodicity " << endl;
	cout << xtrunc << " atoms were read in that needed to be truncated at x boundary " << endl;
	cout << yshift << " atoms were read in that needed to be shifted for y periodicity " << endl;
	cout << ytrunc << " atoms were read in that needed to be truncated at y boundary " << endl;
	cout << zshift << " atoms were read in that needed to be shifted for z periodicity " << endl;
	cout << ztrunc << " atoms were read in that needed to be truncated at z boundary " << endl;
	
}

//***************************************************************************************************
void GetCells(double **Coords, int **AtomID, int **CellList, int *CellPop, double lXDIM, double lYDIM, double lZDIM, int lXCELLS, int lYCELLS, int lZCELLS){
	/*
	 Assigns each atom to a local box, and saves those indices to AtomID
	 Then, make lists of atoms that are in each Cell
	 */
	
	int i, ij, ik, j, maxpop, numcell;
	bool test;
	
	cout << "entering getcells" << endl;
	numcell=lXCELLS*lYCELLS*lZCELLS;
	
	//debug
	cout << "numcells= "<<numcell<<endl;
	cout << lXDIM << " " <<lYDIM <<" " <<lZDIM <<endl;

	for(i = 0; i < numcell; i++){
		CellPop[i] = 0;
		for(j = 0; j < LMAX; j++){
			CellList[i][j] = 0;
		}
	}
	
	// now put things in boxes
	for(i = 0; i < natom; i++){
			
		AtomID[i][2] = int(Coords[i][0]/lXDIM);
		AtomID[i][3] = int(Coords[i][1]/lYDIM);
		AtomID[i][4] = int(Coords[i][2]/lZDIM);

        // This logic check should now be obsolete, but the code runs fast...
		if(AtomID[i][2]<0) cout << "ix too little: i = " << i << endl;
		if(AtomID[i][2]>lXCELLS-1) cout << "ix too big: i = " << i << endl;
		if(AtomID[i][3]<0) cout << "iy too little: i = " << i << endl;
		if(AtomID[i][3]>lYCELLS-1) cout << "iy too big: i = " << i << endl;
		if(AtomID[i][4]<0) cout << "iz too little: i = " << i << endl;
		if(AtomID[i][4]>lZCELLS-1) cout << "iz too big: i = " << i << endl;

		ij = AtomID[i][4] + AtomID[i][3]*lZCELLS + AtomID[i][2]*lZCELLS*lYCELLS;
		ik = CellPop[ij]; //no. of atoms in the cell

		CellList[ij][ik] = i; // ith atom read in for jth timestep
		CellPop[ij] += 1;
	}
	
	maxpop = 0;
	test = false;
	for(i = 0; i < numcell; i++){
		maxpop = max(maxpop,CellPop[i]);
		if(CellPop[i]>LMAX-1){
			cout << " # atoms per exceeded for cell = " << i <<  " with CellPop = " << CellPop[i] << endl;
			test = true;
		}
	}
	if(test) exit(1);
	cout << "exiting get cells: maxpop = " << maxpop << " LMAX = " << LMAX << endl;
}

//***************************************************************************************************

void GetNeighsandVectors(double ***NVecs, double **Coords, int **LNeigh, double **BDist2, int **CellList, int *NNeigh, int *CellPop){
	/*
	 Loop over atoms in each cell
	 Check self and 13 unique neighboring cells for neighboring atoms within the cutoff length
     Loop over all cells completes the neighbor list
	 */
	int i, in1, in2;
	int j, k, k2, kx, ky, kz, l1, l2, lm1, lm2, maxneigh, numcell;
	int n, nx, ny, nz;
	double r2, rad2, xij, yij, zij, rxij, ryij, rzij;
	double xs, ys, zs;
	bool doit, test;
	int o, p;

	rad2 = radius*radius;
	numcell= XCELLS*YCELLS*ZCELLS;
	

	for(i = 0; i < natom; i++){
		NNeigh[i] = 0;
		for(j = 0; j < LMAX; j++){
			LNeigh[i][j] = 0;
			BDist2[i][j] = 0.0;
			NVecs[i][j][0] = 0.0;
			NVecs[i][j][1] = 0.0;
			NVecs[i][j][2] = 0.0;
		}
	}
	
	cout << "entering GetNeighsandVectors" << endl;
	// loop over all the cells
	for(k = 0; k < numcell; k++){

	//get the cell population			
		lm1 = CellPop[k];
		kz  = k%ZCELLS;
		ky  = int((k-kz)/ZCELLS)%YCELLS;
		kx  = int((int((k-kz)/ZCELLS)-ky)/YCELLS);
			
		// loop over pairs of atoms within that cell, including self
		if(lm1>0){
			for(l1 = 0; l1 < lm1; l1++){
				// first do self
				in1 = CellList[k][l1];  // pointer to ith atom read in
				
				LNeigh[in1][NNeigh[in1]] = in1; // particular neighbor atom
				BDist2[in1][NNeigh[in1]] = 0.0;
				NVecs[in1][NNeigh[in1]][0] = 0.0;
				NVecs[in1][NNeigh[in1]][1] = 0.0;
				NVecs[in1][NNeigh[in1]][2] = 0.0;
				NNeigh[in1] += 1; // total neighbors to this atom including itself
				// now loop over neighbors
				for(l2 = l1+1; l2 < lm1; l2++){
					in2 = CellList[k][l2];
					xij = Coords[in1][0] - Coords[in2][0];
					yij = Coords[in1][1] - Coords[in2][1];
					zij = Coords[in1][2] - Coords[in2][2];
					
					r2  = pow(xij,2) + pow(yij,2) + pow(zij,2);

					if(r2 < rad2){
						LNeigh[in1][NNeigh[in1]] = in2; // particular atom number
						LNeigh[in2][NNeigh[in2]] = in1;
						BDist2[in1][NNeigh[in1]] = r2;
						BDist2[in2][NNeigh[in2]] = r2;
						NVecs[in1][NNeigh[in1]][0] = xij;
						NVecs[in2][NNeigh[in2]][0] = -xij;
						NVecs[in1][NNeigh[in1]][1] = yij;
						NVecs[in2][NNeigh[in2]][1] = -yij;
						NVecs[in1][NNeigh[in1]][2] = zij;
						NVecs[in2][NNeigh[in2]][2] = -zij;
						NNeigh[in1] += 1; // total neighbors to this atom
						NNeigh[in2] += 1; // incrementing after the fact for storage in slot 0
					}
				}
			}
		}
			
		// loop over neighboring cells
		//			exb = 0;
		for(n = 0; n < 13; n++){
			nx = kx + dx[n];
			ny = ky + dy[n];
			nz = kz + dz[n];
			xs = 0.0;
			ys = 0.0;
			zs = 0.0;
			doit = true;
			if(nx>XCELLS-1){
				if(XPER){
					nx = 0;
					xs = XLEN;
				}else{
					doit = false;
				}
			}
			if(nx<0){
				if(XPER){
					nx = XCELLS-1;
					xs = -XLEN;
				}else{
					doit = false;
				}
			}
			if(ny>YCELLS-1){
				if(YPER){
					ny = 0;
					ys = YLEN;
				}else{
					doit = false;
				}
			}
			if(ny<0){
				if(YPER){
					ny = YCELLS-1;
					ys = -YLEN;
				}else{
					doit = false;
				}
			}
			if(nz>ZCELLS-1){
				if(ZPER){
					nz = 0;
					zs = ZLEN;
				}else{
					doit = false;
				}
			}
			if(nz<0){
				if(ZPER){
					nz = ZCELLS-1;
					zs = -ZLEN;
				}else{
					doit = false;
				}
			}
			if(doit){
				k2  = nz + ny*ZCELLS + nx*YCELLS*ZCELLS;
				lm2 = CellPop[k2];
				
				// loop over atoms in primary cell
                if((lm1>0)&&(lm2>0)){
                    for(l1 = 0; l1 < lm1; l1++){
                        for(l2 = 0; l2 < lm2; l2++){
                            in1 = CellList[k][l1];  // pointer to this atom, any timestep
                            in2 = CellList[k2][l2];

                            xij = Coords[in1][0] - Coords[in2][0] - xs;
                            yij = Coords[in1][1] - Coords[in2][1] - ys;
                            zij = Coords[in1][2] - Coords[in2][2] - zs;
						
                            r2  = pow(xij,2) + pow(yij,2) + pow(zij,2);
						
                            if(r2 < rad2){
                                LNeigh[in1][NNeigh[in1]] = in2; // particular bond type
                                LNeigh[in2][NNeigh[in2]] = in1;
                                BDist2[in1][NNeigh[in1]] = r2;
                                BDist2[in2][NNeigh[in2]] = r2;
                                NVecs[in1][NNeigh[in1]][0] = xij;
                                NVecs[in2][NNeigh[in2]][0] = -xij;
                                NVecs[in1][NNeigh[in1]][1] = yij;
                                NVecs[in2][NNeigh[in2]][1] = -yij;
                                NVecs[in1][NNeigh[in1]][2] = zij;
                                NVecs[in2][NNeigh[in2]][2] = -zij;
                                NNeigh[in1] += 1; // total neighbors to this atom
                                NNeigh[in2] += 1;
                            }
						}
					}
				}
			}
		}
	}
	maxneigh = 0;
	test = false;
	for(i = 0; i < natom; i++){
		maxneigh = max(maxneigh,NNeigh[i]);
		if(NNeigh[i]>LMAX-1){
			test = true;
			cout << " # neighbors per atom exceeded for i = " << i <<  " with NNeigh = " << NNeigh[i] << endl;
		}
	}
	if(test) exit(1);
	cout << "exiting GetNeighsandVectors: maxneigh = " << maxneigh << " LMAX = " << LMAX << endl;
}

//***************************************************************************************************
void GetNeighs(double **Coords, int **LNeigh, double **BDist2, int **CellList, int *NNeigh, int *CellPop, int lXCELLS, int lYCELLS, int lZCELLS, double lXLEN, double lYLEN, double lZLEN){
	/*
	 Loop over atoms in each cell
	 Check self and 13 unique neighboring cells for neighobring atoms within cutoff length
     Loop over all cells completes the neighbor list
	 */
	int i, in1, in2;
	int j, k, k2, kx, ky, kz, l1, l2, lm1, lm2, maxneigh, numcell;
	int n, nx, ny, nz;
	double r2, rad2, xij, yij, zij, rxij, ryij, rzij;
	double xs, ys, zs;
	bool doit, test;
	int o, p;

	rad2 = radius*radius;
	numcell= lXCELLS*lYCELLS*lZCELLS;
	
	for(i = 0; i < natom; i++){
		NNeigh[i] = 0;
		for(j = 0; j < LMAX; j++){
			LNeigh[i][j] = 0;
			BDist2[i][j] = 0.0;
		}
	}
	
	cout << "entering GetNeighs" << endl;
	// loop over all the cells
	for(k = 0; k < numcell; k++){

	//get the cell population			
		lm1 = CellPop[k];
		
		kz  = k%lZCELLS;
		ky  = int((k-kz)/lZCELLS)%lYCELLS;
		kx  = int((int((k-kz)/lZCELLS)-ky)/lYCELLS);
			
		// loop over pairs of atoms within that cell, including self
		if(lm1>0){
			for(l1 = 0; l1 < lm1; l1++){
				// first do self
				in1 = CellList[k][l1];  // pointer to ith atom read in
				LNeigh[in1][NNeigh[in1]] = in1; // particular neighbor atom
				BDist2[in1][NNeigh[in1]] = 0.0;
				NNeigh[in1] += 1; // total neighbors to this atom including itself
				// now loop over neighbors
				for(l2 = l1+1; l2 < lm1; l2++){
					in2 = CellList[k][l2];
					
					xij = Coords[in1][0] - Coords[in2][0];
					yij = Coords[in1][1] - Coords[in2][1];
					zij = Coords[in1][2] - Coords[in2][2];
					
					r2  = pow(xij,2) + pow(yij,2) + pow(zij,2);

					if(r2 < rad2){
						LNeigh[in1][NNeigh[in1]] = in2; // particular atom number
						LNeigh[in2][NNeigh[in2]] = in1;
						BDist2[in1][NNeigh[in1]] = r2;
						BDist2[in2][NNeigh[in2]] = r2;
						NNeigh[in1] += 1; // total neighbors to this atom
						NNeigh[in2] += 1; // incrementing after the fact for storage in slot 0
					}
				}
			}
		}
			
		// loop over neighboring cells
		for(n = 0; n < 13; n++){
			nx = kx + dx[n];
			ny = ky + dy[n];
			nz = kz + dz[n];
			xs = 0.0;
			ys = 0.0;
			zs = 0.0;
			doit = true;
			if(nx>lXCELLS-1){
				if(XPER){
					nx = 0;
					xs = lXLEN;
				}else{
					doit = false;
				}
			}
			if(nx<0){
				if(XPER){
					nx = lXCELLS-1;
					xs = -lXLEN;
				}else{
					doit = false;
				}
			}
			if(ny>lYCELLS-1){
				if(YPER){
					ny = 0;
					ys = lYLEN;
				}else{
					doit = false;
				}
			}
			if(ny<0){
				if(YPER){
					ny = lYCELLS-1;
					ys = -lYLEN;
				}else{
					doit = false;
				}
			}
			if(nz>lZCELLS-1){
				if(ZPER){
					nz = 0;
					zs = lZLEN;
				}else{
					doit = false;
				}
			}
			if(nz<0){
				if(ZPER){
					nz = lZCELLS-1;
					zs = -lZLEN;
				}else{
					doit = false;
				}
			}
			if(doit){
				k2  = nz + ny*lZCELLS + nx*lYCELLS*lZCELLS;
				lm2 = CellPop[k2];
				
				// loop over atoms in primary cell
                if((lm1>0)&&(lm2>0)){
                    for(l1 = 0; l1 < lm1; l1++){
                        for(l2 = 0; l2 < lm2; l2++){
                            in1 = CellList[k][l1];  // pointer to this atom, any timestep
                            in2 = CellList[k2][l2];

                            xij = Coords[in1][0] - Coords[in2][0] - xs;
                            yij = Coords[in1][1] - Coords[in2][1] - ys;
                            zij = Coords[in1][2] - Coords[in2][2] - zs;
						
                            r2  = pow(xij,2) + pow(yij,2) + pow(zij,2);
						
                            if(r2 < rad2){
                                LNeigh[in1][NNeigh[in1]] = in2; // particular bond type
                                LNeigh[in2][NNeigh[in2]] = in1;
                                BDist2[in1][NNeigh[in1]] = r2;
                                BDist2[in2][NNeigh[in2]] = r2;
                                NNeigh[in1] += 1; // total neighbors to this atom
                                NNeigh[in2] += 1;
                            }
						}
					}
				}
			}
		}
	}
	
	maxneigh = 0;
	test = false;
	for(i = 0; i < natom; i++){
		maxneigh = max(maxneigh,NNeigh[i]);
		if(NNeigh[i]>LMAX-1){
			test = true;
			cout << " # neighbors per atom exceeded for i = " << i <<  " with NNeigh = " << NNeigh[i] << endl;
		}
	}
	if(test) exit(1);
	cout << "exiting GetNeighs: maxneigh = " << maxneigh << " LMAX = " << LMAX << endl;
	
	cout << "leaving GetNeighs " << endl;
}

//***************************************************************************************************
void GetVectors(double ***NVecs, double **Coords, int *NNeigh, int **LNeigh){
		
    // define the distance vectors between all neighboring atoms.
	// pass the neighbor list for the current or the reference configuration
    
	int i, j, in1, in2;
	double xij, yij, zij, xs, ys, zs;

	for(i = 0; i < natom; i++){
		for(j = 0; j < LMAX; j++){
			NVecs[i][j][0] = 0.0;
			NVecs[i][j][1] = 0.0;
			NVecs[i][j][2] = 0.0;
		}
	}

	xs=0.0;
	ys=0.0;
	zs=0.0;
		
	cout << "entering GetVectors" << endl;
	//loop over the atoms
	for(i = 0; i < natom; i++){
		
		//loop over the neighbors of the ith atom
		for(j=0; j<NNeigh[i]; j++){
			//get the id of the neighbor atom
			in1 = LNeigh[i][j];
			//calculate distance vectors
			xij = Coords[i][0] - Coords[in1][0];
			yij = Coords[i][1] - Coords[in1][1];
			zij = Coords[i][2] - Coords[in1][2];
		//do the minimum imgae convention

			if(XPER){
				if(xij <= -(0.5*XLEN)){
					xij+=XLEN;
				}else if(xij > (0.5*XLEN)){
					xij-=XLEN;
				}
			}
			if(YPER){
				if(yij <= -(0.5*YLEN)){
					yij+=YLEN;
				}else if(yij > (0.5*YLEN)){
					yij-=YLEN;
				}
			}
			if(ZPER){
				if(zij <= -(0.5*ZLEN)){
					zij+=ZLEN;
				}else if(zij > (0.5*ZLEN)){
					zij-=ZLEN;
				}
			}	

			NVecs[i][j][0] = xij;
			NVecs[in1][NNeigh[in1]][0] = -xij;
			NVecs[i][j][1] = yij;
			NVecs[in1][NNeigh[in1]][1] = -yij;
			NVecs[i][j][2] = zij;
			NVecs[in1][NNeigh[in1]][2] = -zij;
		}
	}
			
	
	cout << "exiting GetVectors" << endl;
}

//***************************************************************************************************
void CalcMoments(double ***NVecs, double ***NVels, double ***NFors, double **Coords, int **AtomID, int **LNeigh, double **BDist2, int *NNeigh){
	
	 // First, calculate Cartesian Moments for all atoms with their neighbors
	 // Then transform into complex angular momentum forms and store
	
    FILE *opfile, *pstfile;

	int i, i2, j, k, l, m, n, neighnum;
	double x_a, y_a, z_a;
	double scale, sigmi, wii, V0i;
	double Wij[LMAX];
	double pPk, pPkl, pPklm, pPklmn;

    // temporary Cartesian tensors
	double P_1[3];
	double P_2[3][3];
	double P_3[3][3][3];
    double P_4[3][3][3][3];

    // Invariant tensors
	double P0_Inv;
	double P1_Inv[1];
	double P2_Inv[3];
	double P3_Inv[7];
	double P4_Inv[12];
	double O1_Inv[3];
	double O2_Inv[3];
	double O3_Inv[4];
	double O4_Inv[5];

    // temporary vectors
	double P0[5];
    complex<double> P1[2], P2[3];
    complex<double> P3_1a[2], P3_3a[4];
	complex<double> P4_2a[3], P4_4a[5];

    double sqrt2i = sqrt(0.5);
    double rt2pi  = sqrt(6.2831853072);
    double sqrt3o2 = sqrt(1.5);
    double sqrt5 = sqrt(5.);
    double sqrt7o2 = sqrt(3.5);
    double sqrt15o8 = sqrt(15./8.);
    
    // Orientation tensors
    complex<double> Z[2] = {complex<double>(1.,0.),complex<double>(0.,0.)};
    complex<double> X[2] = {complex<double>(0.,0.),complex<double>(-sqrt2i,0.)};
    complex<double> Y[2] = {complex<double>(0.,0.),complex<double>(0.,-sqrt2i)};
    complex<double> Z2[3] = {complex<double>(1.,0.),complex<double>(0.,0.),complex<double>(0.,0.)};
    complex<double> Y2[3] = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,sqrt2i)};
    complex<double> X2[3] = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(sqrt2i,0.)};
    complex<double> Z3[4]  = {complex<double>(1.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.)};
    complex<double> X3[4]  = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(sqrt2i,0.)};
    complex<double> Y3[4]  = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,sqrt2i)};
    complex<double> Z4[5]   = {complex<double>(1.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.)};
    complex<double> X4[5]   = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(sqrt2i,0.)};
    complex<double> Y4[5]   = {complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,0.),complex<double>(0.,sqrt2i)};
    
	cout << " Entering CalcMoments " << endl;
	
	scale = 1./(2.*sigma*sigma);
	sigmi = 1./sigma;
	V0i  = 1./pow(rt2pi*sigma,3);
	
	opfile = fopen(MomFileName, "w");
	pstfile = fopen("Position_ST.txt", "w");

	if(!opfile){
		cout << "Error in WriteMoments: Cannot open moment file "
		<< MomFileName << endl;
		exit(1);
	}
	
    // writing in LAMMPS dump/data text format
    // readable by OVITO
	fprintf(opfile, "ITEM: TIMESTEP \n");
        fprintf(opfile, "%d\n", timestep);
        fprintf(opfile, "ITEM: NUMBER OF ATOMS \n");
        fprintf(opfile, "%d\n", natom);
        fprintf(opfile, "ITEM: BOX BOUNDS %s %s %s \n", xdum,ydum,zdum);
        fprintf(opfile, "%lf %lf\n", XLO, XHI);
        fprintf(opfile, "%lf %lf\n", YLO, YHI);
        fprintf(opfile, "%lf %lf\n", ZLO, ZHI);
		fprintf(opfile,"ITEM: ATOMS id x y z w P0_I0 P1_I0 "
        //fprintf(opfile,"id x y z w P0_I0 P1_I0 "
                        "P2_I0 P2_I1 P2_I2 "
                        "P3_I0 P3_I1 P3_I2 P3_I3 P3_I4 P3_I5 P3_I6 "
                        "P4_I0 P4_I1 P4_I2 P4_I3 P4_I4 P4_I5 P4_I6 P4_I7 P4_I8 P4_I9 "
						"P4_I10 P4_I11 "
                        "O1_I0 O1_I1 O1_I2 "
                        "O2_I0 O2_I1 O2_I2 "
                        "O3_I0 O3_I1 O3_I2 "
                        "O4_I0 O4_I1 O4_I2 "
                        "\n");

	fprintf(pstfile, "ITEM: TIMESTEP \n");
	fprintf(pstfile, "%d\n", timestep);
	fprintf(pstfile, "ITEM: NUMBER OF ATOMS \n");
	fprintf(pstfile, "%d\n", natom);
	fprintf(pstfile, "ITEM: BOX BOUNDS %s %s %s \n", xdum,ydum,zdum);
	fprintf(pstfile, "%lf %lf\n", XLO, XHI);
	fprintf(pstfile, "%lf %lf\n", YLO, YHI);
	fprintf(pstfile, "%lf %lf\n", ZLO, ZHI);
	fprintf(pstfile,"ITEM: ATOMS id x y z "
            "P00_0 P11_0 P11_1r P11_1i "
            "P20_0 P22_0 P22_1r P22_1i P22_2r P22_2i "
            "P31a_0 P31a_1r P31a_1i "
            "P33a_0 P33a_1r P33a_1i P33a_2r P33a_2i P33a_3r P33a_3i "
            "P40_0 P42a_0 P42a_1r P42a_1i P42a_2r P42a_2i "
			"P44a_0 P44a_1r P44a_1i P44a_2r P44a_2i P44a_3r P44a_3i P44a_4r P44a_4i"
            "\n");

	for (i = 0; i < natom; i++) {
		
		// Now that we have the positions of all the neighbors relative to the base_atom,
		// we can compute the center of position (this is equivalent to eliminating the
		// first order moment).  The difference in position between the central atom and
		// center of position will be stored as the P_1 moment.
		//
		
		P0[0] = 0.0;

		x_a  = 0.;
		y_a  = 0.;
		z_a  = 0.;

		for(j=0; j<NNeigh[i]; j++)
		{	
			Wij[j]  = exp(-BDist2[i][j]*scale);
			x_a    += Wij[j]*NVecs[i][j][0];
			y_a    += Wij[j]*NVecs[i][j][1];
			z_a    += Wij[j]*NVecs[i][j][2];
			i2      = LNeigh[i][j];
			P0[0]  += Wij[j];
		}
		wii = 1./P0[0];
		// Normalize to define weighted average position for neighborhood
		x_a *= wii;
		y_a *= wii;
		z_a *= wii;
		
		// Switch all of the coordinates s.t. the origin is now at r_a and scale the distance by sigma
		//
		for(j=0; j<NNeigh[i]; j++)
		{
			NVecs[i][j][0] = (NVecs[i][j][0] - x_a)*sigmi;  // Normalizing distances by sigma
			NVecs[i][j][1] = (NVecs[i][j][1] - y_a)*sigmi;
			NVecs[i][j][2] = (NVecs[i][j][2] - z_a)*sigmi;
		}
        
		// initialize orders just to make sure
		for (k = 0; k < 3; k++) {
			P_1[k] = 0.;  // already done
			for (l = 0; l < 3; l++) {
				P_2[k][l] = 0.;
				for (m = 0; m < 3; m++) {
					P_3[k][l][m] = 0.;
					for (n = 0; n < 3; n++) {
						P_4[k][l][m][n] = 0.;
					}
				}
			}
		}
		
		//  store scaled origin shift in P_1
		P_1[0] = x_a*sigmi;
		P_1[1] = y_a*sigmi;
		P_1[2] = z_a*sigmi;
	
		neighnum=NNeigh[i];

		for(j=0; j<neighnum; j++) //FOR-LOOP OVER NEIGHBORS
		{
			for (k = 0; k < 3; k++) {
				//P_1[k] += Wij[j]*NVecs[i][j][k]; // Don't overwrite these anymore.
				pPk = Wij[j]*NVecs[i][j][k];
				for (l = 0; l < 3; l++) {
					pPkl = pPk*NVecs[i][j][l]; // r*r really only need upper triangle for r
					P_2[k][l] += pPkl; // r*r really only need upper triangle for r
					for (m = 0; m < 3; m++) {
						pPklm = pPkl*NVecs[i][j][m]; // r*r*r
						P_3[k][l][m] += pPklm; // r*r*r
						for (n = 0; n < 3; n++) {
							pPklmn = pPklm*NVecs[i][j][n]; //
							P_4[k][l][m][n] += pPklmn; // r*r*r*r
						}
					}
				}
			}
		} //END FOR-LOOP OVER NEIGHBORS
		
		//  Normalize by total weight stored in W
		for (k = 0; k < 3; k++) {
			for (l = 0; l < 3; l++) {
				P_2[k][l] *= wii;
				for (m = 0; m < 3; m++) {
					P_3[k][l][m] *= wii;
					for (n = 0; n < 3; n++) {
						P_4[k][l][m][n] *= wii;
					}
				}
			}
		}

        // convert Cartesian tensors to ang mom vectors, spherical tensors
		// simple vectors
		P1_2_ST(&P_1[0], P1);

        // 3x3 matrix
		// For fully, antisymmetric matrix, this yields one scalar, one vector, one symmetric matrix
		// For symmetric matrix (positions), this yields one scalar and one symmetric, traceless matrix
		P2_2_ST(&P_2[0][0], P2, P0);
		P0[1] -= sqrt3o2; //*P0[0]; arrgh, that weighting was already normalized out

		// 3x3x3 tensor: for fully asymmetric form should yield
		// one scalar, three traces (1 vector), two symmetric matrices (2 vector), one symmetric tensor (3 vector)
		// For semi-symmetric: no scalar, two unique traces (1 vector), one symmetric matrix and one sym 3-tensor
		// For fully symmetric: no scalar, one unique trace (vector), no symmetric matrix and one sym 3-tensor
		P3_2_ST(&P_3[0][0][0], P3_3a, P3_1a);

		// 3x3x3x3 tensor: for fully asymmetric form should yield
        // three scalars, six vectors, six symmetric matrices, three symmetric 3-tensors, one symmetric 4-tensor
        // For semi-symmetric: no scalar, two unique traces (1 vector), one symmetric matrix and one sym 3-tensor
        // For fully symmetric: one scalar, no vectors, one symmetric matrix, no sym 3-tensor, one sym 4-tensor
		P4_2_ST(&P_4[0][0][0][0], P4_4a, P4_2a, P0);
		P0[2] += -sqrt5*(P0[1] + sqrt3o2) + sqrt15o8; // Colin's
		for (k=0; k<3; k++) {
			P4_2a[k] -= sqrt7o2*P2[k];
		}
        
		fprintf(pstfile,"%d %lf %lf %lf "
                "%lf %lf %lf %lf "
				"%lf %lf %lf %lf %lf %lf "
				"%lf %lf %lf "
                "%lf %lf %lf %lf %lf %lf %lf "
				"%lf %lf %lf %lf %lf %lf "
                "%lf %lf %lf %lf %lf "
                "%lf %lf %lf %lf\n",
				AtomID[i][0], Coords[i][0], Coords[i][1], Coords[i][2],
				P0[0], real(P1[0]), real(P1[1]), imag(P1[1]),
				P0[1], real(P2[0]),  real(P2[1]), imag(P2[1]), real(P2[2]), imag(P2[2]),
				real(P3_1a[0]), real(P3_1a[1]), imag(P3_1a[1]),
                real(P3_3a[0]), real(P3_3a[1]), imag(P3_3a[1]), real(P3_3a[2]), imag(P3_3a[2]), real(P3_3a[3]), imag(P3_3a[3]),
				P0[2], real(P4_2a[0]), real(P4_2a[1]), imag(P4_2a[1]), real(P4_2a[2]), imag(P4_2a[2]),
				real(P4_4a[0]), real(P4_4a[1]), imag(P4_4a[1]), real(P4_4a[2]), imag(P4_4a[2]),
				real(P4_4a[3]), imag(P4_4a[3]), real(P4_4a[4]), imag(P4_4a[4]));
        
//*******************now calculate the invariants*************************************/
		
		// zero order invariant of position: sum of weights => density
		//P0_Inv = P0[0]*V0i;
        P0_Inv = P0[0];

		// first order invariants: magnitude of the vectors
		P1_Invar(P1, P1_Inv, Z, X, Y, O1_Inv);
		P2_Invar(P2, P0, P2_Inv, Z2, Y2, X2, O2_Inv);
		P3_Invar(P3_3a, P3_1a, P3_Inv, P1, Z3, X3, Y3, O3_Inv);
		P4_Invar(P4_4a, P4_2a, P0, P2, P4_Inv, Z4, X4, Y4, O4_Inv);
		
		fprintf(opfile,
            "%d %lf %lf %lf "
            "%e %e %e "
            "%e %e %e "
            "%e %e %e %e %e %e %e "
            "%e %e %e %e %e %e %e %e %e %e "
            "%e %e "
            "%e %e %e "
            "%e %e %e "
            "%e %e %e "
            "%e %e %e\n",
            AtomID[i][0],Coords[i][0],Coords[i][1],Coords[i][2],
            P0[0], P0_Inv, P1_Inv[0],
            P2_Inv[0], P2_Inv[1], P2_Inv[2],
            P3_Inv[0], P3_Inv[1], P3_Inv[2], P3_Inv[3], P3_Inv[4], P3_Inv[5], P3_Inv[6],
            P4_Inv[0], P4_Inv[1], P4_Inv[2], P4_Inv[3], P4_Inv[4], P4_Inv[5], P4_Inv[6], P4_Inv[7], P4_Inv[8], P4_Inv[9],
            P4_Inv[10],P4_Inv[11],
            O1_Inv[0], O1_Inv[1], O1_Inv[2],
            O2_Inv[0], O2_Inv[1], O2_Inv[2],
            O3_Inv[0], O3_Inv[1], O3_Inv[2],
            O4_Inv[0], O4_Inv[1], O4_Inv[2]);
	}
    
	fclose(opfile);
	fclose(pstfile);

}

//***************************************************************************************************

//***************************************************************************************************

int main(int argc, char** argv){

	/*
	Main Function
     Overall structure is to read in an atomic coordinate data file and put it on an orthorhombic grid
     If data is a periodic orthorhombic structure, the periodicity can be handled
     Non-orthorhombic structures can not currently be recognized as periodic
     A cut-off distance is defined based on Gaussian weighting function << 10^-6
     Box sizes are defined based on that cutoff distance
     Looping over each box and its neighbors is used to define the atomic neighbor list for each atom
     The Cartesian moments for each atom and its neighbors is then calculated
     These are transformed to Spherical Tensors and the invariants are then calculated
     The oriented Spherical Tensors are output to the ST file
     The rotational invariant tensors (the SFD) and their orientation factors are output to the Moment file
     Both output files are written in a LAMMPS data/dump text format and can be visualized by OVITO
	*/
    
	FILE *pfile;
	FILE *boxfile;
	
	int k, l;
	int dummy;
	double cutoff;
	
	radius 	= 5.6;
	sigma  	= 1.;
	BINYES = false;
	REFWEIGHT = false;


	// The Params file defines the name of the data files
    // and the Gaussian sigma value and associated cutoff radius
    // and if a reference configuration to evaluate changes is being used.
	if(!(pfile = fopen(ParFileName,"r"))){
		cout << "No parameter file found " << ParFileName << " Using the defaults" << endl;
		cout << "Default radius = " << radius << "Default sigma = " << sigma << endl;
		exit(1);
	}else{
		cout << "Reading Parameters file " << ParFileName << endl;
		fscanf(pfile,"%127s\n",XYZFileName);
		cout << "Input File " << XYZFileName << endl;
		fscanf(pfile,"%d\n",&natom);
		cout << "Number of Atoms " << natom << endl;
		
		fscanf(pfile,"%lf\n",&radius);		
		fscanf(pfile,"%lf\n",&sigma);
		cout << "reading in radius = " << radius << " sigma = " << sigma << endl;
		fscanf(pfile,"%d\n", &dummy);
		if(dummy==1) BINYES = true;
		fscanf(pfile,"%d\n", &dummy);
		if(dummy==1) REFWEIGHT = true;
        fscanf(pfile,"%d\n", &DSTYLE);
        cout << "DSTYLE = " << DSTYLE << endl;

        
		if(REFWEIGHT){
                        cout <<"Weighting using reference configuration " <<endl;
		}

		if(BINYES){
                        cout <<"Reading coords from binary file "<< BINFileName <<endl;
                        getSavFileNAtoms();
                        cout <<"Reading box dimension from the .parm file "<< BoxFileName <<endl;
                        boxfile = fopen(BoxFileName, "r");
                        fscanf(boxfile,"%lf\n", &XLO);
                        fscanf(boxfile,"%lf\n", &YLO);
                        fscanf(boxfile,"%lf\n", &ZLO);
                        fscanf(boxfile,"%lf\n", &XHI);
                        fscanf(boxfile,"%lf\n", &YHI);
                        fscanf(boxfile,"%lf\n", &ZHI);
                        fclose(boxfile);
                }
                else{
                        cout <<"Reading coords from ASCII file "<<endl;
                }
	}

	fclose(pfile);

	cutoff = exp(-radius*radius*0.5/(sigma*sigma));
	cout << "cutoff weighting = " << cutoff << endl;
	

//allocate all the DM
	double ***DM_NVecs	=new double**[natom];
	double ***DM_NVels	=new double**[natom];
	double ***DM_NFors	=new double**[natom];
	
	double **DM_Coords;
	double **DM_RefCoords;
	//float **DM_Coords;
	double **DM_BDist2;
	int **DM_AtomID;
	int **DM_RefAtomID;
	int **DM_LNeigh;
	int **DM_CellList;
	int **DM_RefCellList;

	int *DM_NNeigh;
	int *DM_CellPop;
	int *DM_RefCellPop;
	

	DM_Coords	= new double *[natom];
	DM_RefCoords	= new double *[natom];
	DM_BDist2	= new double *[natom];
	DM_AtomID	= new int *[natom];
	DM_RefAtomID	= new int *[natom];
	DM_LNeigh	= new int *[natom];

	DM_NNeigh	= new int [natom];
	
	for (k =0 ; k <natom ; k++){
		DM_Coords[k] 	= new double [NPATOM];
		DM_RefCoords[k] = new double [NPATOM];
		DM_BDist2[k] 	= new double [LMAX];
		DM_AtomID[k] 	= new int [NAID];
		DM_RefAtomID[k] = new int [NAID];
		DM_LNeigh[k] 	= new int [LMAX];

		DM_NVecs[k]  	= new double*[LMAX];
		DM_NVels[k]  	= new double*[LMAX];
		DM_NFors[k]  	= new double*[LMAX];

		for (l =0 ; l <LMAX ; l++){
			DM_NVecs[k][l] 	= new double[3];
			DM_NVels[k][l] 	= new double[3];
			DM_NFors[k][l] 	= new double[3];
		}
		
	}

    if(BINYES){
            readSavFile(DM_Coords, DM_AtomID);
    }
    else{
            GetCoords(DM_Coords, DM_RefCoords, DM_AtomID, DM_RefAtomID);
    }

//GetCoords will populate the value of ncells and rbcells. Use that to allocate the CellPop and CellList array
	DM_CellPop	= new int [ncells];
	DM_CellList	= new int *[ncells];
	DM_RefCellPop	= new int [rncells];
	DM_RefCellList	= new int *[rncells];

	for (k=0; k<ncells; k++){
		DM_CellList[k]	= new int [LMAX];
	} 
	//for the reference configuration
	for (k=0; k<rncells; k++){
		DM_RefCellList[k]	= new int [LMAX];
	} 
//end allocate the DM

	GetCells(DM_Coords, DM_AtomID, DM_CellList, DM_CellPop, XDIM, YDIM, ZDIM, XCELLS, YCELLS, ZCELLS);

	if(REFWEIGHT){
//do it for the reference configuration
		GetCells(DM_RefCoords, DM_RefAtomID, DM_RefCellList, DM_RefCellPop, RXDIM, RYDIM, RZDIM, RXCELLS, RYCELLS, RZCELLS);
	}
	
	if(REFWEIGHT){
//get neighbors according to the reference configuration
		GetNeighs(DM_RefCoords, DM_LNeigh, DM_BDist2, DM_RefCellList, DM_NNeigh, DM_RefCellPop, RXCELLS, RYCELLS, RZCELLS, RXLEN, RYLEN, RZLEN);
//use neighbors from reference and Coords from current to creat xij vectors
		GetVectors(DM_NVecs, DM_Coords, DM_NNeigh, DM_LNeigh);
	}else{

		GetNeighsandVectors(DM_NVecs, DM_Coords, DM_LNeigh, DM_BDist2, DM_CellList, DM_NNeigh, DM_CellPop);
	}
	
   	ClebschGordan();
    
	CalcMoments(DM_NVecs, DM_NVels, DM_NFors, DM_Coords, DM_AtomID, DM_LNeigh, DM_BDist2, DM_NNeigh);
//	CalcMoments(DM_NVecs, DM_Coords, DM_AtomID, DM_LNeigh, DM_BDist2, DM_NNeigh, DM_P_0);

//free the DM

	for (k =0; k<natom ; k++){
		for (l =0 ; l <LMAX ; l++){
			if (l==0){
				delete [] DM_NVecs[k][l];
				delete [] DM_NVels[k][l];   
				delete [] DM_NFors[k][l];   

				delete [] DM_Coords[k];
				delete [] DM_RefCoords[k];
				delete [] DM_AtomID[k];
				delete [] DM_RefAtomID[k];
				delete [] DM_LNeigh[k];
				delete [] DM_BDist2[k];
				}
			else{
                                delete [] DM_NVecs[k][l];
                		delete [] DM_NVels[k][l];   
				delete [] DM_NFors[k][l];   

                        }
		}
		delete [] DM_NVecs[k];
		delete [] DM_NVels[k];
		delete [] DM_NFors[k];
	}

	delete [] DM_NVecs;
	delete [] DM_NVels;
	delete [] DM_NFors;
	
	delete [] DM_Coords;
	delete [] DM_RefCoords;
	delete [] DM_AtomID;
	delete [] DM_RefAtomID;
	delete [] DM_LNeigh;
	delete [] DM_BDist2;

	delete [] DM_NNeigh;	
	delete [] DM_CellPop;
	delete [] DM_RefCellPop;
	
	for (k=0; k<ncells; k++){
		delete [] DM_CellList[k];
	}
	
	for (k=0; k<rncells; k++){
		delete [] DM_RefCellList[k];
	}

	delete [] DM_CellList;
	delete [] DM_RefCellList;

	return 0;
}
