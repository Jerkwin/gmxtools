//! flag: -O2, -Wall
// exe:
//! inc: F:/gromacs2019.6/src
//! lib: F:/gromacs2019.6/_build/lib/libgromacs.a,F:/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/9.2.0/libgomp.a

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#define pi       3.141592653589793
#define rad2deg 57.29577951308232
#define deg2rad  0.0174532925199433

enum plotType {
	plotQTP, plotTP, plotConf
};

struct SRing {
	float         aXcoords[6];
	float         aYcoords[6];
	float         aZcoords[6];
	float         fQ;
	float         fTheta;
	float         fPhi;
	struct SRing *next;
};

struct SFrame {
	int           nNumber;
	float         fTime;
	struct SRing *rings;
	struct SFrame *next;
};

void cremer_pople(int nring, float *x, float *y, float *z,
				  float *fQ, float *ftheta, float *fphi) {

	/*
	nring: determens if it's a five- or six-membered ring
	x, y, z: a matrix containg all cartesian coordinates of relevant ring atoms

	for hexapyranose ring should pick atoms
		O5,C1,C2,C3,C4,C5
	   for agreement with Cremer & Pople Sucrose pyranoid ring coordinates
	for furanoid ring pick atoms
		O2', C5', C4', C3' C2'
	*/

	int j;
	double xnew[6], ynew[6], znew[6], zdisp[6],
		xmean, ymean, zmean, R1[3], R2[3], n[3],
		s, c, t, rNN, q3, Q, theta, phi2;

	rNN=1./nring;

	/* copy coords to arrays xnew,ynew,znew */
	for(j=0; j<nring; j++) {
		xnew[j]=x[j];
		ynew[j]=y[j];
		znew[j]=z[j];
	}

	/* find mean x,y,z coordinates */
	xmean = 0.; ymean = 0.; zmean = 0.;
	for(j=0; j<nring; j++) {
		xmean += xnew[j];
		ymean += ynew[j];
		zmean += znew[j];
	}
	xmean *= rNN; ymean *= rNN; zmean *= rNN;

	// now subtract mean -> origin at geometric mean of coordinates of ring atoms
	for(j=0; j<nring; j++) {
		xnew[j] -= xmean;
		ynew[j] -= ymean;
		znew[j] -= zmean;
	}

	/* compute R1, R2 vectors */
	R1[0]=0.; R1[1]=0.; R1[2]=0.;
	R2[0]=0.; R2[1]=0.; R2[2]=0.;
	for(j=0; j<nring; j++) {
		t=2*pi*j*rNN;
		s=sin(t);
		c=cos(t);
		R1[0] += xnew[j]*s; R1[1] += ynew[j]*s; R1[2] += znew[j]*s;
		R2[0] += xnew[j]*c; R2[1] += ynew[j]*c; R2[2] += znew[j]*c;
	}
	//printf("vector R1 = %f %f %f \n",R1[0],R1[1],R1[2]);
	//printf("vector R2 = %f %f %f \n",R2[0],R2[1],R2[2]);

	/* compute cross product */
	n[0] = R1[1]*R2[2]-R1[2]*R2[1];
	n[1] = R1[2]*R2[0]-R1[0]*R2[2];
	n[2] = R1[0]*R2[1]-R1[1]*R2[0];
	t=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0] /= t; n[1] /= t; n[2] /= t;
	//printf("n vector = %f %f %f \n",n[0],n[1],n[2]);

	/* compute z axis displacements, EQ 11*/
	for(j=0; j<nring; j++) {
		zdisp[j]=xnew[j]*n[0]+ynew[j]*n[1]+znew[j]*n[2];
		//printf("atom %d %4.4s z coord = %f z disp = %f \n",
		//			i,atom_name[k],znew[i],zdisp[i]);
	}

	/* compute total puckering amplitude, EQ 17*/
	Q=0.0;
	for(j=0; j<nring; j++) Q += zdisp[j]*zdisp[j];
	Q=sqrt(Q);
	//printf("total puckering amplitude = %5.4f \n",Q);

	/* now compute q1, phi1 using Eq 12,13 of
	   q1 should be = 0 hence no need to evaluate */
	phi2=999;
	if (nring==5 || nring==6) {
		/* now compute q2, phi2 using Eq 12,13  */
		s=0.0; c=0.0;
		for(j=0; j<nring; j++) {
			t=4*pi*j*rNN;
			s -= zdisp[j]*sin(t);
			c += zdisp[j]*cos(t);
		}
		t=sqrt(2*rNN);
		s *= t;
		c *= t;
		phi2=atan2(s, c)*rad2deg;
		if(phi2<0.) phi2 += 360;
		//printf("q2 = %8.4f phi 2 = %8.4f degrees \n", q2, phi2);
	}

	theta=999;
	if(nring==6) { /* compute q3 using eq 14 */
		q3=0.0;
		for(j=0; j<nring; j++) {
			if(j%2) q3 -= zdisp[j];
			else    q3 += zdisp[j];
		}
		q3 *= sqrt(rNN);
		theta=acos(q3/Q)*rad2deg;
		//printf("q3 = %8.4f theta = %8.4f degrees \n", q3, theta);
	}

	*fQ=Q;
	*ftheta=theta;
	*fphi=phi2;
}

int getConf(double theta, double phi, char* conf) {
		//"Berce's et al. dihedrals Cremer-Pople coordinates this work",
		//"#     t1   t2   t3    phi theta Q   th0    th1   th2",
	char txt[][100]={
		"1C4   60  -60   60    0  180  0.57 -35.26 -35.26 -35.26 ",
		"4C1  -60   60  -60    0    0  0.57  35.26  35.26  35.26 ",
		"1,4B   0   60  -60  240   90  0.76 -35.26  74.20 -35.26 ",
		"B1,4   0  -60   60   60   90  0.76  35.26 -74.20  35.26 ",
		"2,5B -60    0   60  120   90  0.76  74.20 -35.26 -35.26 ",
		"B2,5  60    0  -60  300   90  0.76 -74.20  35.26  35.26 ",
		"3,6B  60  -60    0    0   90  0.76 -35.26 -35.26  74.20 ",
		"B3,6 -60   60    0  180   90  0.76  35.26  35.26 -74.20 ",
		"1H2   45  -15    0  270  129  0.42 -42.16   9.07 -17.83 ",
		"2H1  -45   15    0   90   51  0.42  42.16  -9.07  17.83 ",
		"2H3  -60   45  -15  150   51  0.42  42.16  17.83  -9.06 ",
		"3H2   60  -45   15  330  129  0.42 -42.16 -17.83   9.06 ",
		"3H4   45  -60   45   30  129  0.42 -17.83 -42.16   9.07 ",
		"4H3  -45   60  -45  210   51  0.42  17.83  42.16  -9.07 ",
		"4H5  -15   45  -60  270   51  0.42  -9.07  42.16  17.83 ",
		"5H4   15  -45   60   90  129  0.42   9.07 -42.16 -17.83 ",
		"5H6    0  -15   45  150  129  0.42   9.07 -17.83 -42.16 ",
		"6H5    0   15  -45  330   51  0.42  -9.07  17.83  42.16 ",
		"6H1  -15    0  -15   30   51  0.42  17.83  -9.07  42.16 ",
		"1H6   15    0   15  210  129  0.42 -17.83   9.07 -42.16 ",
		"1S3  -30   60  -30  210   88  0.62   0     50.84 -50.84 ",
		"3S1   30  -60   30   30   92  0.62   0    -50.84  50.84 ",
		"5S1  -30  -30   60   90   92  0.62  50.84 -50.84   0    ",
		"1S5   30   30  -60  270   88  0.62 -50.84  50.84   0    ",
		"6S2   60  -30  -30  330   88  0.62 -50.84   0     50.84 ",
		"2S6  -60   30   30  150   92  0.62  50.84   0    -50.84 ",
		"1E    30    0    0  240  125  0.45 -35.26  17.37 -35.26 ",
		"E1   -30    0    0   60   55  0.45  35.26 -17.37  35.26 ",
		"2E   -60   30    0  120   55  0.45  46.86   0      0    ",
		"E2    60  -30    0  300  125  0.45 -46.86   0      0    ",
		"3E    60  -60   30  360  125  0.45 -35.26 -35.26  17.37 ",
		"E3   -60   60  -30  180   55  0.45  35.26  35.26 -17.37 ",
		"4E   -30   60  -60  240   55  0.45   0     46.86   0    ",
		"E4    30  -60   60   60  125  0.45   0    -46.86   0    ",
		"5E     0  -30   60  120  125  0.45  17.37 -35.26 -35.26 ",
		"E5     0   30  -60  300   55  0.45 -17.37  35.26  35.26 ",
		"6E     0    0  -30  360   55  0.45   0      0     46.86 ",
		"E6     0    0   30  180  125  0.45   0      0    -46.86 ",};

		theta *= deg2rad;
		phi   *= deg2rad;

		int idx=0;
		double dx, dy, dz, t, p, R=1, d, dmin=1E99;
		char tmp[4], str[99];
		for(int i=0; i<38; i++) {
			sscanf(txt[i], "%s %s %s %s %lf %lf", tmp, str, str, str, &p, &t);
			t *= deg2rad;
			p *= deg2rad;

			dx=R*(sin(theta)*cos(phi)-sin(t)*cos(p));
			dy=R*(sin(theta)*sin(phi)-sin(t)*sin(p));
			dz=R*(cos(theta)         -cos(t));

			d=sqrt(dx*dx+dy*dy+dz*dz);
			if(d<dmin) { dmin=d; idx=i+1; strcpy(conf, tmp); }
			//printf("%d %f %f %f %f %f %f %f\n", i, t*rad2deg, p*rad2deg, dx, dy, dz, d, dmin);
		}

		return idx;
}

void writeTimePlot(const char *fn, int nRings, char *grpnm, int plotType,
					struct SFrame *firstFrame, const gmx_output_env_t *oenv) {
	/*
	* fn       filename
	* nRings   Number atoms of the ring
	* grpnm    Name of Indexgroup
	* firstFrame pointer to Linked List
	*/

	int    i;
	char   conf[4], Title[STRLEN], Xlabel[STRLEN], Ylabel[STRLEN];
	FILE   *out;
	struct SFrame *currentFrame;
	struct SRing  *currentRing;

	switch(plotType) {
		case plotQTP:
			sprintf(Title, "Cremer-Pople Puckering Parameters for Group: %s", grpnm);
			strcpy(Xlabel, "Time (ps)");
			strcpy(Ylabel, "Q(pm), \\f{Symbol}q\\f{}(бу), \\f{Symbol}j\\f{}(бу)");
			break;
		case plotTP:
			sprintf(Title, "Cremer-Pople Puckering Polar Plot for Group %s", grpnm);
			strcpy(Xlabel, "\\f{Symbol}j\\f{}(бу)");
			strcpy(Ylabel, "\\f{Symbol}q\\f{}(бу)");
			break;
		case plotConf:
			sprintf(Title, "Canonical Ring Conformation for Group: %s", grpnm);
			strcpy(Xlabel, "Time (ps)");
			strcpy(Ylabel, "Canonical Ring Index Number and Name");
			break;
	}

	out=xvgropen(fn, Title, Xlabel, Ylabel, oenv);

	fprintf(out, "@ view 0.1, 0.1, 1.5, 0.88\n");
	fprintf(out, "@ legend on\n");
	fprintf(out, "@ legend box on\n");
	fprintf(out, "@ legend loctype view\n");
	fprintf(out, "@ legend 0.78, 0.8\n");
	fprintf(out, "@ legend length 2\n");

	switch(plotType) {
		case plotQTP:
			for(int i=0; i<nRings; i++) {
				fprintf(out, "@ s%0d legend \"Q\\s %0d\"\n", 3*i, i);
				fprintf(out, "@ s%0d legend \"\\f{Symbol}q\\f{}\\s %0d\"\n", 3*i+1, i);
				fprintf(out, "@ s%0d legend \"\\f{Symbol}j\\f{}\\s %0d\"\n", 3*i+2, i);
			}
			break;
		case plotTP:
			xvgr_line_props(out, 0, elNone, ecFrank, oenv);
			fprintf(out, "@ s0 symbol 2\n@ s0 symbol size 0.5\n@ s0 symbol fill 1\n");
			break;
		case plotConf:
			break;
	}

	currentFrame=firstFrame;
	while(currentFrame) {
		if(plotType!=plotTP) fprintf(out, "%12.6f", currentFrame->fTime);

		currentRing=currentFrame->rings;
		while(currentRing) {

			switch(plotType) {
				case plotQTP:
					fprintf(out, " %12.6f%9.3f%9.3f",
						currentRing->fQ*1000, currentRing->fTheta, currentRing->fPhi);
					break;
				case plotTP:
					fprintf(out, "%9.3f%9.3f", currentRing->fPhi, currentRing->fTheta);
					break;
				case plotConf:
					i=getConf(currentRing->fTheta, currentRing->fPhi, conf);
					fprintf(out, "%3d(%s)", i, conf);
					break;
			}

			currentRing=currentRing->next;
		}
		fprintf(out, "\n");

		currentFrame=currentFrame->next;
	}

	fclose(out);
}

void writeDistPlot(const char *fn, char *grpnm,
				  struct SFrame *firstFrame, const gmx_output_env_t *oenv) {
	/**
	* fn       filename
	* grpnm    Name of Indexgroup
	* distTheta   pointer to Array with Theta-Values
	* distPhi     pointer to Array with Phi-Values
	*/

	int    i, nValues;
	FILE   *out;
	struct SFrame *currentFrame;
	struct SRing  *currentRing;
	double distTheta[360], distPhi[360], t;

	char   Title[STRLEN];
	char   Xlabel[]="\\f{Symbol}q\\f{}(бу), \\f{Symbol}j\\f{}(бу)";
	char   Ylabel[]="Relative Probability [%]";

	for(i=0; i<360; ++i) { distTheta[i]=0.; distPhi[i]=0.; }

	/* Calculate Distribution Function :*/
	nValues=0;
	currentFrame=firstFrame;
	while(currentFrame) {
		currentRing=currentFrame->rings;
		while(currentRing) {
			nValues++;
			distTheta[(int)round(currentRing->fTheta)]++;
			distPhi[(int)round(currentRing->fPhi)]++;
			currentRing=currentRing->next;
		}
		currentFrame=currentFrame->next;
	}

	sprintf(Title, "Cremer-Pople Puckering Angles Distribution for Group %s", grpnm);
	out=xvgropen(fn, Title, Xlabel, Ylabel, oenv);

	fprintf(out, "@ view 0.1, 0.1, 1.5, 0.88\n");
	fprintf(out, "@ legend on\n");
	fprintf(out, "@ legend box on\n");
	fprintf(out, "@ legend loctype view\n");
	fprintf(out, "@ legend 0.78, 0.8\n");
	fprintf(out, "@ legend length 2\n");
	fprintf(out, "@ s0 legend \"\\f{Symbol}q\"\n");
	fprintf(out, "@ s1 legend \"\\f{Symbol}j\"\n");

	t=100./nValues;
	for(i=0; i<360; ++i)
		fprintf(out,"%4i %6.3f %6.3f\n", i, distTheta[i]*t, distPhi[i]*t);

	fclose(out);
}

void selfTest() {
  float cptheta = 0.0; /* Just to make sure */
  float cpphi   = 0.0;
  float cpQ     = 0.0;
  int   ring, ringsize, atom;
  float coords_x[6], coords_y[6], coords_z[6];
  float testCoords[3][3][6] = {
    /* Sample Molecules from Cremer & Pople, JACS, 1975, 97(6), p. 1354*/
    { /* Furanoid Ring of Sucrose */
      {  0.0000,  1.1622,  0.7425, -0.7221, -1.1826, 0.0000 }, /* X coordinates in A */
      {  1.2111,  0.4349, -1.0012, -1.0309,  0.3861, 0.0000 }, /* Y coordinates in A */
      { -0.0189,  0.1461, -0.2174,  0.2057, -0.1154, 0.0000 }  /* Z coordinates in A */
      /* should give:
       * CREMER POPLE theta  =   n/a
       * CREMER POPLE phi2   =   265.1   [░]
       * CREMER POPLE Q      =     0.353 [A]
       */
    },
    { /* Pyranoid Ring of Sucrose */
      {  0.0000,  1.1997,  1.2356,  0.0110, -1.2300, -1.2164 }, /* X coordinates in A */
      {  1.3839,  0.7624, -0.7040, -1.4564, -0.7208,  0.7350 }, /* Y coordinates in A */
      {  0.1976, -0.2106,  0.2393, -0.2550,  0.2420, -0.2133 }  /* Z coordinates in A */
      /* should give:
       * CREMER POPLE theta  =     5.2   [░]
       * CREMER POPLE phi2   =   183.7   [░]
       * CREMER POPLE Q      =     0.556 [A]
       */
    },
    /* Sample Molecule from mdxview */
    { /* Pyranoid Ring of Clucose */
      { -1.323,  -0.509,   1.076,   1.301,   0.395,  -1.125  }, /* X coordinates in A */
      { -0.617,  -0.675,  -0.531,   0.688,   0.693,   0.521  }, /* Y coordinates in A */
      {  0.156,   1.337,   1.004,   0.186,  -1.070,  -0.701  }  /* Z coordinates in A */
      /* should give:
       * CREMER POPLE theta  =     4.48  [░]бу
       * CREMER POPLE phi2   =   167.38  [░]
       * CREMER POPLE Q      =     0.536 [A]
       */
    }
  };
  fprintf(stdout, "***************\n*  TEST DATA  *\n***************\n");
  for ( ring = 0; ring < asize(testCoords); ring++) {
    if (testCoords[ring][0][5]==0.0 && testCoords[ring][1][5]==0.0 && testCoords[ring][1][5]==0.0000)
      ringsize=5;
    else
      ringsize=6;

    for (atom = 0; atom < 6; atom++) {
      coords_x[atom] = testCoords[ring][0][atom];
      coords_y[atom] = testCoords[ring][1][atom];
      coords_z[atom] = testCoords[ring][2][atom];
      float dist=0.0;
      if (!(ringsize==5 && atom==5)) {
        dist = sqrt(
            ((coords_x[atom]-coords_x[(atom+1)%ringsize])*(coords_x[atom]-coords_x[(atom+1)%ringsize]))
          + ((coords_y[atom]-coords_y[(atom+1)%ringsize])*(coords_y[atom]-coords_y[(atom+1)%ringsize]))
          + ((coords_z[atom]-coords_z[(atom+1)%ringsize])*(coords_z[atom]-coords_z[(atom+1)%ringsize]))
            );
        fprintf(stdout, "Atom %i:\t%8.3f%8.3f%8.3f, %8.4f\n", atom, coords_x[atom], coords_y[atom], coords_z[atom], dist);
      }
    }

    cremer_pople(ringsize, coords_x, coords_y, coords_z,  &cpQ, &cptheta, &cpphi);
    fprintf(stdout, "CREMER POPLE theta  = %8.1f \n",cptheta);
    fprintf(stdout, "CREMER POPLE phi2   = %8.1f \n",cpphi);
    fprintf(stdout, "CREMER POPLE Q      = %8.5f \n",cpQ);
  }
}

int main(int argc, char *argv[]) {

	static char Copyright[7][STRLEN] = {
		"                                                                       \n",
		"                        :-)  gmx_puckering  (-:                        \n",
		"                                                                       \n",
		"                       written 2007 by Oliver Stueker                  \n",
		"                       revised 2021 by Jicun  Li                       \n",
		"                       GNU General Public License                      \n",
		"                                                                       \n",
	};
	for(int i=0; i<7 ; i++) fprintf(stderr, Copyright[i]);

	/**********************
	* INIT Variables      *
	**********************/

	const char *desc[] = {
		"gmx_puckering computes the Cremer-Pople-Ring-Puckering Parameters of Pyranoses and Hexanoses.\n",
		"\n",
		"The index file needs to contain atom-sextuples or atom-quintuples in the order:\n",
		"* O5 C1 C2 C3 C4 C5 for Hexanoses or \n",
		"* O2 C5 C4 C3 C2    for Pyranoses\n",
		"as defined by Cremer and Pople in\n",
		"  Cremer, D.; Pople, J. A., General definition of ring puckering coordinates. ",
		"  J. Am. Chem. Soc. 1975, 97, (6), 1354-1358.\n",
		"\n",
		"If the number of atoms in the group is not divisible by the ringsize given ",
		"by -i the program will give an Error and exit. The program will also check ",
		"if the atomnames match to the scheme given above and give a warning if the ",
		"names don't match, which can be suppressed by using the -noname option.\n",
		"\n",
		"Following plots are available:\n",
		"-o    cp_Q-theta-phi.xvg    Q/Theta/Phi vs. Time\n",
		"-otp    cp_theta-phi.xvg    Theta vs. Phi\n",
		"-od   cp_dtheta-dphi.xvg    Distribution of Theta/Phi\n",
		"-or          cp_ring.xvg    IUPAC Canonical Ring Conformation\n",
		"\n",
		"\n",
		"Acknowledgments:\n",
		" * GROMACS - http://www.gromacs.org\n",
		" * g_puckering of Oliver Stueker - http://www.gromacs.org/Downloads/User_contributions/Other_software \n",
		" * mdxvu of Mark J. Forster - http://sourceforge.net/projects/mdxvu/\n",
	};

	/* Extra arguments - but note how you always get the begin/end
	* options when running the program, without mentioning them here!
	*/

	/** for TOPOLOGY **/
	t_topology top;

	/** for TRAJECTORY **/
	rvec      *xtop;
	matrix    box;
	int       frame = 0;
	int       natoms;
	real      t;
	rvec      *x;

	/** for INDEX **/
	t_atoms     *atoms;
	char        *grpnm;
	int         *index;
	int         nidx = 0, maxIndex = 0;

	/** for STATUS **/
	t_trxstatus      *status;
	gmx_output_env_t *oenv;

	int  ePBC;
	gmx_rmpbc_t gpbc;

	/** for WARNING **/
	int  maxwarn = 20;
	warninp_t wi = init_warning(TRUE, maxwarn);

	/** for ARGUMENTS **/
	int  ringsize = 6;

	bool bVerbose = FALSE;
	bool bTestMode = FALSE;
	bool bDistCheck=TRUE;
	bool bNameCheck=TRUE;
	bool bThetaPhi = 0, bDist = 0, bRing=0;

	/** for CP Parameters **/
	int         i;
	int         nRings;
	int         ring, atom;
	float       cptheta, cpphi, cpQ;
	char        sWarning[STRLEN];

	struct SFrame *first_Frame;
	struct SFrame *current_Frame;
	struct SFrame *new_Frame;
	struct SFrame *prev_Frame;
	struct SRing  *current_Ring;
	struct SRing  *new_Ring;
	struct SRing  *prev_Ring;

	t_pargs pa[] = {
		{ "-i",       FALSE, etINT,  {&ringsize},   "Size of Ring: 5 or 6" },
		{ "-v",       FALSE, etBOOL, {&bVerbose},   "Be loud and noisy" },
		{ "-dist",    FALSE, etBOOL, {&bDistCheck}, "Warn if distance between neighboring Ringatoms is larger that 0.3 nm." },
		{ "-name",    FALSE, etBOOL, {&bNameCheck}, "Warn if Atomnames don't match with Definition by Cremer & Pople." },
		{ "-maxwarn", FALSE, etINT,  {&maxwarn},    "HIDDENNumber of warnings after which input processing stops" },
		{ "-t",       FALSE, etBOOL, {&bTestMode},  "HIDDENUse Test-Modus"}
	};

	std::vector<std::string> strArr(1," ");
	t_filenm fnm[] = {
		{ efTPS, NULL,    NULL,            ffREAD , strArr},  /* this is for the topology */
		{ efTRX, "-f",    NULL,            ffREAD , strArr},  /* and this for the trajectory */
		{ efNDX, "-n",    NULL,            ffREAD , strArr},  /* and this for the index-file */
		{ efXVG, "-o",   "cp_Q-theta-phi", ffOPTWR, strArr},  /* Q/theta/phi vs. time Plot*/
		{ efXVG, "-otp", "cp_theta-phi",   ffOPTWR, strArr},  /* theta vs. phi Plot*/
		{ efXVG, "-od",  "cp_dtheta-dphi", ffOPTWR, strArr},
		{ efXVG, "-or",  "cp_ring",        ffOPTWR, strArr}
	};

	/**********************
	* Parse Arguments     *
	**********************/
	#define NFILE asize(fnm)

	/* This is the routine responsible for adding default options,
	 * calling the X/motif interface, etc. */
	if(!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT,
		NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv))
		gmx_fatal(FARGS,"Error: Invalid Option(s).");

	bThetaPhi = opt2bSet("-otp",  NFILE, fnm);
	bDist= opt2bSet("-od",  NFILE, fnm);
	bRing= opt2bSet("-or",  NFILE, fnm);

	/**********************
	* Read Topology       *
	**********************/
	read_tps_conf(ftp2fn(efTPS, NFILE, fnm), &top, &ePBC, &xtop, NULL, box, TRUE);
	sfree(xtop);

	/**********************
	* Read Index File     *
	**********************/
	atoms=&(top.atoms);
	get_index(atoms, ftp2fn_null(efNDX, NFILE, fnm), 1, &nidx, &index, &grpnm);

	/**********************
	* Checking            *
	**********************/
	/* Check if No of Elements matches with number of ring atoms */
	nRings=nidx/ringsize;
	fprintf(stdout, "Will be checking %3i Rings.\n", nRings);
	if((nidx % ringsize) != 0)
		gmx_fatal(FARGS,
		"number of index elements not multiple of %d, these can not be %s\n",
		ringsize,(ringsize==6) ? "six-membered rings" : "five-membered rings");
	else if(bVerbose) fprintf(stdout, "OK! Number of index elements matches.\n" );

	/**********************
	* Setup Data Arrays  *
	**********************/
	fprintf(stdout, "Number of Rings: %i\n", nRings);
	int **ringAtoms;
	ringAtoms=(int**)malloc(nRings*sizeof(int));
	for(i=0; i<ringsize; i++) ringAtoms[i]=(int*)malloc(ringsize*sizeof(int));
	//[nRings][ringsize];

	/**********************
	* Checking            *
	**********************/
	/* Check if all atoms in Index-File are found in Topology */
	if(bVerbose) fprintf(stdout, "Index consists of items:" );
	for(i=0; i<nidx; i++) {
		ringAtoms[i/ringsize][i%ringsize] = index[i];
		if(bVerbose) {
			if(! i%ringsize) fprintf(stdout, "\n");
			fprintf(stdout, " %6i", index[i] );
		}
		if( index[i] > maxIndex) maxIndex = index[i];
	}
	if(bVerbose) fprintf(stdout, "\n" );
	sfree(index);

	if(maxIndex<0 || maxIndex>(top.atoms.nr))
		gmx_fatal(FARGS, "Error: Atom number %d is out of range.\n", maxIndex);

	if(bVerbose) fprintf(stdout, "OK! Index doesen't exeed Topology.\n");

	if(bNameCheck) {
	/* Check if Atomnames match to
	*    hexapyranose ring (O5,C1,C2,C3,C4,C5)
	* or furanoid ring (O2', C5', C4', C3' C2')
	*/
		for(ring=0; ring<nRings; ++ring) {
			if(ringsize==6) {
				if( strcasecmp( *(top.atoms.atomname[ringAtoms[ring][0]]), "O5" )
				 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][1]]), "C1" )
				 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][2]]), "C2" )
				 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][3]]), "C3" )
				 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][4]]), "C4" )
				 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][5]]), "C5" )) {
					sprintf(sWarning,
						"Names of Atoms specified in Index Group do not match specification for hexapyranose rings.\nIs: (%s,%s,%s,%s,%s,%s), should be: (O5,C1,C2,C3,C4,C5)",
						*(top.atoms.atomname[ringAtoms[ring][0]]),
						*(top.atoms.atomname[ringAtoms[ring][1]]),
						*(top.atoms.atomname[ringAtoms[ring][2]]),
						*(top.atoms.atomname[ringAtoms[ring][3]]),
						*(top.atoms.atomname[ringAtoms[ring][4]]),
						*(top.atoms.atomname[ringAtoms[ring][5]]) );
					warning(wi, sWarning);
				}
			} else if(!(strcasecmp( *(top.atoms.atomname[ringAtoms[ring][0]]), "O2" )
					 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][1]]), "C5" )
					 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][2]]), "C4" )
					 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][3]]), "C3" )
					 || strcasecmp( *(top.atoms.atomname[ringAtoms[ring][4]]), "C2" )) ) {
				sprintf(sWarning, "Names of Atoms specified in Index Group do not match specification for furanoid rings.\nIs: (%s,%s,%s,%s,%s), should be: (O2,C5,C4,C3,C2)",
					*(top.atoms.atomname[ringAtoms[ring][0]]),
					*(top.atoms.atomname[ringAtoms[ring][1]]),
					*(top.atoms.atomname[ringAtoms[ring][2]]),
					*(top.atoms.atomname[ringAtoms[ring][3]]),
					*(top.atoms.atomname[ringAtoms[ring][4]]) );
				warning(wi, sWarning);
			}
		}
	} /*end if (bNameCheck)*/

	frame = 0;
	/* The first time we read data is a little special */
	if(bVerbose) fprintf(stdout, "Now we'll Read the first Frame:\n" );
	natoms=read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
	printf("\n");

	first_Frame=(struct SFrame *)malloc(sizeof(struct SFrame));
	if(first_Frame == NULL) {
		gmx_mem("Cannot allocate SFrame.");
		return 1;
	}

	current_Frame=first_Frame;
	prev_Frame=NULL;
	/* This is the main loop over frames */
	current_Frame->next=NULL;
	do {
		gpbc = gmx_rmpbc_init(&top.idef, ePBC, natoms);
		gmx_rmpbc(gpbc, natoms, box, x);

		current_Frame->fTime=t;
		current_Frame->nNumber=frame;
		current_Frame->rings=(struct SRing *)malloc(sizeof(struct SRing));
		if(current_Frame->rings == NULL) {
			gmx_mem("Cannot allocate SRing.");
			return 1;
		}

		current_Ring=current_Frame->rings;
		prev_Ring=NULL;
		for(ring=0; ring<nRings; ring++) { /* Loop over rings in frame */
			current_Ring->fTheta = 0.0;
			current_Ring->fPhi   = 0.0;
			current_Ring->fQ     = 0.0;
			for(atom=0; atom<ringsize; atom++) {
			/* copy coordinates of all ring atoms into the Data Structure */
				current_Ring->aXcoords[atom] = x[ringAtoms[ring][atom]][XX] ;
				current_Ring->aYcoords[atom] = x[ringAtoms[ring][atom]][YY] ;
				current_Ring->aZcoords[atom] = x[ringAtoms[ring][atom]][ZZ] ;
			}

			if(bDistCheck) {
				/* Distance Check: Warn if distance between two neighboring Ring-Atoms is > 0.3nm.  */
				for(atom=0; atom<ringsize; atom++) {
					float dX, dY, dZ, dist;
					dX= current_Ring->aXcoords[atom] - current_Ring->aXcoords[(atom+1)%ringsize];
					dY= current_Ring->aYcoords[atom] - current_Ring->aYcoords[(atom+1)%ringsize];
					dZ= current_Ring->aZcoords[atom] - current_Ring->aZcoords[(atom+1)%ringsize];
					dist = sqrt( (dX * dX) + (dY * dY) + (dZ * dZ) );

					if(dist>0.3) {
						sprintf(sWarning, "Distance between Atoms %4i (%s) and %4i (%s) is very large (%5.3G nm).",
							ringAtoms[ring][atom]+1,              *(top.atoms.atomname[ringAtoms[ring][atom]]),
							ringAtoms[ring][(atom+1)%ringsize]+1, *(top.atoms.atomname[ringAtoms[ring][(atom+1)%ringsize]]),
							dist);
						warning(wi, sWarning);
						printf("X1: %8.4f, X2: %8.4f, => dX: %f\n", current_Ring->aXcoords[atom], current_Ring->aXcoords[(atom+1)%ringsize], dX);
						printf("Y1: %8.4f, Y2: %8.4f, => dY: %f\n", current_Ring->aYcoords[atom], current_Ring->aYcoords[(atom+1)%ringsize], dY);
						printf("Z1: %8.4f, Z2: %8.4f, => dZ: %f\n", current_Ring->aZcoords[atom], current_Ring->aZcoords[(atom+1)%ringsize], dZ);
					}
				}
			}

			cremer_pople(ringsize, current_Ring->aXcoords,
						 current_Ring->aYcoords,
						 current_Ring->aZcoords,
						 &cpQ, &cptheta, &cpphi);

			current_Ring->fTheta = cptheta;
			current_Ring->fPhi   = cpphi;
			current_Ring->fQ     = cpQ;
			if(bVerbose) {
				fprintf(stdout, "CREMER POPLE theta  = %8.2f \n",cptheta);
				fprintf(stdout, "CREMER POPLE phi2   = %8.2f \n",cpphi);
				fprintf(stdout, "CREMER POPLE Q      = %8.5f \n",cpQ);
			}

			/* Append a new Ring to List */
			new_Ring = (struct SRing *)malloc(sizeof(struct SRing));
			if (new_Ring == NULL) {
				gmx_mem("Cannot allocate SRing.");
				return 1;
			}

			current_Ring->next=new_Ring;
			prev_Ring=current_Ring;
			current_Ring=new_Ring;
		}
		/* Terminate List and clean up */
		prev_Ring->next=NULL;
		free(current_Ring);
		current_Ring=NULL;

		/* Append new Frame to the List */
		new_Frame = (struct SFrame *)malloc(sizeof(struct SFrame));
		if(new_Frame == NULL) {
			gmx_mem("Cannot allocate SFrame.");
			return 1;
		}
		current_Frame->next=new_Frame;
		prev_Frame=current_Frame;
		current_Frame=new_Frame;

		frame++;
	} while(read_next_x(oenv, status, &t, x, box));

	/* Terminate List and clean up */
	prev_Frame->next=NULL;
	free(current_Frame);
	current_Frame=NULL;

	if(bTestMode) selfTest();

	/**********************
	* Output XVG FILES   *
	**********************/
	fprintf(stdout, "Ready for Output CP puckering Parameters for group: %s\n", grpnm);
	if(ringsize==5) {
		sprintf(sWarning, "Puckering Angle Theta is not defined for five-membered rings.\n");
		warning(wi, sWarning);
	}

	writeTimePlot(opt2fn("-o",NFILE,fnm), nRings, grpnm, plotQTP, first_Frame, oenv);
	if(bThetaPhi) writeTimePlot(opt2fn("-otp",NFILE,fnm), nRings, grpnm, plotTP, first_Frame, oenv);
	if(bRing)     writeTimePlot(opt2fn("-or",NFILE,fnm),  nRings, grpnm, plotConf, first_Frame, oenv);

	if(bDist)     writeDistPlot(opt2fn("-od", NFILE,fnm), grpnm, first_Frame, oenv);
}

/*
*
*  Created by Pan Chen on 10/04/15.
*  Copyright 2011 CERMAV. All rights reserved.
*

#include <valarray>
#include <iostream>
using namespace std;
int mainfc(int argc, char *argv[])
{
	int n = atoi(argv[1]);
	valarray<double> v(n*3);
	double x, y, z;
	for(int i = 0; i < n; i++){
		cin >> x >> y >> z;
		v[i*3]   = x;
		v[i*3+1] = y;
		v[i*3+2] = z;
	}
	valarray<double> center(3);
	for(int i = 0; i < n; i++) center += v[slice(i*3, 3, 1)];
	center/=n;
	for(int i = 0; i < n; i++) v[slice(i*3, 3, 1)] -=center;
	valarray<double> s(n);
	valarray<double> c(n);
	for(int i = 0; i < n; i++) {
		s[i] = sin(2*M_PI*i/n);
		c[i] = cos(2*M_PI*i/n);
		cout <<s[i]<<" "<<c[i]<<endl;
	}
	valarray<double> r1(3);
	valarray<double> r2(3);
	valarray<double> r3(3);
	for(int i = 0; i < n; i++){
		valarray<double> temp(v[slice(i*3, 3, 1)]);
		valarray<double> temp1(v[slice(i*3, 3, 1)]);
		cout << temp[0]<<" "<<temp[1]<<" "<<temp[2]<<endl;
		temp *= s[i];
		temp1 *= c[i];
		r1 += temp;
		r2 += temp1;
	}
	cout << r1[0] << " "<<r1[1]<< " "<<r1[2]<<endl;
	cout << r2[0] << " "<<r2[1]<< " "<<r2[2]<<endl;
	cout << center[0] << " "<<center[1]<< " "<<center[2]<<endl;
	for(int i = 0; i < 3; i++){
		r3[i] = r1[(i+1)%3]*r2[(i+2)%3]-r2[(i+1)%3]*r1[(i+2)%3];
	}
	double rabs = sqrt(r3[0]*r3[0]+ r3[1]*r3[1] * r3[2]*r3[2]);
	r3/=rabs;
	cout << r3[0] << " "<<r3[1]<<" "<<r3[2]<<endl;
	valarray<double> zv(n);
	for(int i = 0; i < n; i++) {
		zv[i] = v[i* 3]*r3[0] + v[i*3+1] * r3[1] + v[i*3+2] * r3[2];
		cout <<zv[i]<<endl;
	}
	valarray<double> param(n-3);
	if(n > 3){
		for(int i = 0; i < (n-4)/2.; i++){
			for(int j = 0; j < n; j++){
				param[i*2] += cos(2*M_PI*(i+2)*j/n)*zv[j];
				param[i*2+1] += cos(2*M_PI*(i+2)*j/n)*zv[j];
			}
		}
		if(!(n%2)){
            int sign = 1;
            for(int j = 0; j < n; j++,sign*=-1)
				param[n-4] += zv[j]*sign;
		}
	}
	double q2 = sqrt(param[0] * param[0] + param[1] * param[1]);
	double theta = atan(q2/param[2]);
	cout <<param[0] << " "<<param[1]<< " " << param[2] <<" "<< theta * 180/M_PI <<endl;
	cout << q2 << endl;
}
*/
