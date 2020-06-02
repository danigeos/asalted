/*
	"Two peanuts are walking on the dry Mediterranean, 
	then one is assaulted... "

	This program calculates the incision and water flow across a strait 
	communicating two basins. 
	Author: Daniel Garcia-Castellanos, 2010-2011, danielgc(AT)ictja.csic.es 
	License: This software is distributed under a Creative Commons BY-SA license. 
	You can use it and distribute it without consulting the author. You must
	cite the authorship and share your modifications alike. 
	Scientific citation: Garcia-Castellanos & Villaseñor, 2011, Nature. 

	This program is the successor of 'spillover' (Garcia-Castellanos et al., 2009, Nature). 

	Present-day	volume, area, salinity, salt:
	Mediterranean Sea: 3.77 10^6 km3 2.5e12  m2                 3.9%   1.42e17 kg (Haq et al., 2020)
	Oceans (world):    1340 10^6 km3 3.61e14 m2 (Gleick, 1993)  3.5%?  4.69e19 kg

    Present net evaporation: 0.8-1.1 m/yr, ~75,000 m3/s, similar to:
    Present water deficit in the Mediterranean of 2.4e12 m3/yr; Bryden & Stommel, 1984; Blanc, 2006

	Messinian halite mass precipitated in the Med.: 
	Blanc, 2000:       6.44e18 kg (E) 1.44e18 kg (W) 7.88e18 kg (total)
	Ryan, 2008:        3.24e18 kg (E) 1.08e18 kg (W) 4.32e18 kg (total)
	Haq et al., 2020:  1.34e18 kg (E) 0.55e18 kg (W) 1.89e18 kg (total)

	Limitations: 
	-If no uplift and basins are connected and water flowing, zsill
	decreases potentially to minus inf.  This is so even when using m=1. 
	-Second sill is not eroded (yet)

	For halite and gypsum precipitation see Bitzer, 2004. Geologica Acta, 
	Vol. 2, No 4, 2004, 321 - 331
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
//#include <malloc.h>
#include <time.h>

#define	MAXLENLINE	1024			/*Max. length for character strings, input lines, ...*/
#define TAKE_LINE_2(x, y)	{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<2) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f",       &x, &y);};      	if (lin==NULL) break;}
#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t\t2009, Daniel Garcia-Castellanos, danielgc@ictja.csic.es\n");}
#define MAX_2(x,y)	(((x)>(y))? (x) : (y))	/*Gives maximum of two values*/
#define MIN_2(x,y)	(((x)<(y))? (x) : (y))	/*Gives minimum of two values*/
#define	GYPSUMPRECIPCN	.140 			/*Salt concentration for CaSO4 precipitation, kg/l=kg/kg, Blanc, 2006, says .170*/
#define	GYPSUMCNRIVERS	.001 			/*Gypsum concentration in rivers, kg/l=kg/kg*/
#define	GYPSUMCNSEA	.0031 			/*Gypsum concentration in present oceans, kg/l=kg/kg, Bitzer, 2004*/
#define	HALITEPRECIPCN	.3715 			/*Salt concentration for NaCl precipitation, kg/l=kg/kg, E.g., Blanc, 2006*/
#define	HALITECNRIVERS	.005 			/*Halite concentration in rivers, kg/l=kg/kg*/
#define	HALITECNSEA	.035 			/*Halite concentration in present oceans, kg/l=kg/kg; NASW is .03618, Blanc, 2006*/
#define secsperyr 	(365.24  *24*3600)	/*Converts years to seconds*/

#define CAPTION "#time[yr] zsill[m] Rh[m] Slope\tV[m/s]\tQ[m3/s] \tW[m]\te[m/yr] vol0[km3] \tvol1[km3] \tvol2[km3] \tvoltr[km3] \tz0[m]\tz1[m]\tz2[m] gyppr0[10^12kg] gyppr1 gyppr2 gypcn0[kg/l] gypscn1 gypsumcn2  halpr0[10^12kg] halpr1 halpr2 halcn0[kg/l] halcn1 halcn2"
#define PI		3.1415927

int level_and_area_from_volume (float *hypso_z, float *hypso_a, int np, float vol, float *z, float *area);
int volume_and_area_from_level (float *hypso_z, float *hypso_a, int np, float z, float *vol, float *area);
int read_file_insolation(char *filename, int *n_insolation_input_points, float *insolation_mean);

int verbose_level=1;
float **insolation_var; 

main(int argc, char **argv)
{
	int	i, iarg, 
		model_eros=0, model_vel=0, 
		np0=0, np1=0, np2=0, npmax=2000, n_insolation_input_points=0, 
		nbasins=1, water_conservative=0, decreasable_width=1, turowski_width=0, linear_width=0, 
		switch_ps=1;
	float	time, 
		timeini = 0, 
		timeend = -100, 
		dt = 1*secsperyr,
		hl1 = -1000,	hl2 = -1000,
		z_sill1 = -10,	z_sill2 = -430,
		depth0=-3796, depth1 = -2500,	depth2 = -9999, 
		dist1 = 60e3, 	dist2 = 60e3, 
		Ke = 8e-6, 	expe = 1.5, 
		Kw = 1.2, 	expw = .5, 
		e0=0, e1=0, e2=0, 
		p0=0, p1=0, p2=0, 
		r0=0, r1=0, r2=0, 
		roughness = .05, 
		uplift_rate=0, uplift_incr=0, 
		erostotal=0,
		z0=0, z1=0, z2=0, 
		*hypso0_z, *hypso0_a,
		*hypso1_z, *hypso1_a,
		*hypso2_z, *hypso2_a, 
		g = 9.81, denswater = 1020, 
		gypsumds0,       gypsumds1,       gypsumds2, 
		gypsumpr0=0,     gypsumpr1=0,     gypsumpr2=0, 
		gypsumcn0=GYPSUMCNSEA, gypsumcn1=GYPSUMCNSEA, gypsumcn2=GYPSUMCNSEA,  /*kg_s/kg_w = kg/l*/
		haliteds0,       haliteds1,       haliteds2, 
		halitepr0=0,     halitepr1=0,     halitepr2=0, 
		halitecn0=HALITECNSEA, halitecn1=HALITECNSEA, halitecn2=HALITECNSEA,  /*kg_s/kg_w = kg/l*/
		mix01_ratioperyr = 0, mix12_ratioperyr = 0, silldepthmixmax = 0, silldepthmixmin = 0, /*present exchange in Gibraltar: 1.6e6 m3/s; Med vol.=3.615e15 m3; => mixrate=0.014*/
		maxshearstr1 = 0, maxshearstr2 = 0, shearcritical = 0, 
		area0, area1, area2, 
		voltot0=1.37e18, voltot1=1.5237e15, voltot2=1.247e15,  	/*world oceans, west Med, east Med*/
		voltr=0, vol0, vol1, vol2, 
		sl_amp=0, sl_per=0,
		insolation_mean=0, EPfac; 
	FILE 	*file = NULL;

	if (argc<=2) {syntax(argc, argv); exit;}
	
	
	hypso0_z = (float *) calloc(npmax, sizeof(float));
	hypso0_a = (float *) calloc(npmax, sizeof(float));
	hypso1_z = (float *) calloc(npmax, sizeof(float));
	hypso1_a = (float *) calloc(npmax, sizeof(float));
	hypso2_z = (float *) calloc(npmax, sizeof(float));
	hypso2_a = (float *) calloc(npmax, sizeof(float));
	
	/*Interpreting command line*/
	for (iarg=1; iarg<argc; iarg++) {
		if (argv[iarg][0] == '-') {
			int 	ilet;
			float 	value, value2;
			char 	prm[MAXLENLINE], filename[MAXLENLINE], *ptr;

			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);

			switch (argv[iarg][1]) {
				case 'B':
					z_sill1 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) hl1 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dist1 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					sscanf(ptr, "%s", filename);
					if (ptr != NULL) voltot1 = atof(ptr); 	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) depth1 = atof(ptr);
					else {
						if ((file = fopen(filename, "rt")) == NULL) {
							if (verbose_level>=3) fprintf(stderr, "\nWarning: Input basin1 hypsometry file '%s' not found.", filename);
						}
						else {
							if (verbose_level>=3) fprintf(stderr, "\nHypsometry file: '%s'.", filename);
							for (i=0;;i++) {
								TAKE_LINE_2(hypso1_z[i], hypso1_a[i])
								if (hypso1_a[i]==0) depth1 = hypso1_z[i];
							}
							np1=i;
						}
						break;
					}
					ptr=strtok(NULL, "/");
					if (ptr != NULL && *ptr != 'l') {
						if (*ptr == 'b') {
							/*Default box-like hypsometry*/
							hypso1_z[0]=depth1; hypso1_a[0]=0;
							hypso1_z[1]=depth1; hypso1_a[1]= voltot1/(z_sill1-depth1);
							np1 = 2;
							if (verbose_level>=3) fprintf(stderr, "\nAutomatic box-like hypsometry points for basin1: %d", np1);
						}
					}
					else {
						/*Default linear hypsometry*/
						hypso1_z[0]=depth1; hypso1_a[0]=0;
						hypso1_z[1]=z0; hypso1_a[1]= voltot1/(z_sill1-depth1);
						np1 = 2;
						if (verbose_level>=2) fprintf(stderr, "\nAutomatic linear hypsometry points for basin1: %d", np1);
					}
					z1 = depth1; 
					break;
				case 'b':
					nbasins = 2;
					z_sill2 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) hl2 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dist2 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					sscanf(ptr, "%s", filename);
					if (ptr != NULL) voltot2 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) depth2 = atof(ptr);	
					else {
						if ((file = fopen(filename, "rt")) == NULL) {
							if (verbose_level>=3) fprintf(stderr, "\nWarning: Input basin2 hypsometry file '%s' not found.", filename);
						}
						else {
							if (verbose_level>=3) fprintf(stderr, "\nHypsometry file: '%s'.", filename);
							for (i=0;;i++) {
								TAKE_LINE_2(hypso2_z[i], hypso2_a[i])
								if (hypso2_a[i]==0) depth2 = hypso2_z[i];
							}
							np2=i;
						}
						break;
					}					
					ptr=strtok(NULL, "/");
					if (ptr != NULL && *ptr != 'l') {
						if (*ptr == 'b') {
							/*Default box-like hypsometry*/
							hypso2_z[0]=depth2; hypso2_a[0]=0;
							hypso2_z[1]=depth2; hypso2_a[1]= voltot2/(z_sill2-depth2);
							np2 = 2;
							if (verbose_level>=3) fprintf(stderr, "\nAutomatic box-like hypsometry points for basin2: %d", np2);
						}
					}
					else {
						/*Default linear hypsometry*/
						hypso2_z[0]=depth2; hypso2_a[0]=0;
						hypso2_z[1]=z_sill2; hypso2_a[1]= voltot2/(z_sill1-depth2);
						np2 = 2;
						if (verbose_level>=2) fprintf(stderr, "\nAutomatic linear hypsometry points for basin2: %d", np2);
					}
					z2 = depth2; 
					break;
				case 'e':
					e0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) e1 = atof(ptr);
					else e1=e0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) e2 = atof(ptr);
					else e2=e1;
					break;
				case 'h':
					if (argc>2) syntax(argc, argv);
					timeend = dt*5;
					break;
				case 'i':
					//fprintf(stderr, "\n%s ||| %s" , prm , filename);
					strcpy(filename, "--");
					sscanf(prm, "%s", filename);
					if (strlen(prm)>2) read_file_insolation(filename, &n_insolation_input_points, &insolation_mean);
					else fprintf(stderr, "\nERROR: invalid insolation filename: %s", filename);
					//fprintf(stderr, "\n$$$$$$$ %d %.2f W/m2 %.2f W/m2", n_insolation_input_points, insolation_mean, insolation_var[0][0]);
					break;
				case 'k':
					Ke = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) expe = atof(ptr);
					ptr=strtok(NULL, "/");
					if (ptr != NULL) shearcritical = atof(ptr);
					break;
				case 'M':
					model_eros = atoi(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) model_vel = atoi(ptr);
					break;
				case 'm':
					mix01_ratioperyr = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) mix12_ratioperyr = atof(ptr);
					ptr=strtok(NULL, "/");
					if (ptr != NULL) silldepthmixmax = atof(ptr);
					ptr=strtok(NULL, "/");
					if (ptr != NULL) silldepthmixmin = atof(ptr);
					break;
				case 'P':
					switch_ps = 1;
					break;
				case 'p':
					p0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) p1 = atof(ptr);
					else p1=p0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) p2 = atof(ptr);
					else p2=p1;
					break;
				case 'R':
					roughness = value;
					break;
				case 'r':
					r0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) r1 = atof(ptr);
					else r1=r0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) r2 = atof(ptr);
					else r2=r1;
					break;
				case 'S':
					sscanf(prm, "%s", filename);
					voltot0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) 
						depth0 = atof(ptr); 
					else {
						if ((file = fopen(filename, "rt")) == NULL) {
							if (verbose_level>=3) fprintf(stderr, "\nWarning: Input basin0 hypsometry file '%s' not found. Infinite basin0 assumed", filename);
							voltot0 = 0; 
						}
						else {
							if (verbose_level>=3) fprintf(stderr, "\nHypsometry file: '%s'.", filename);
							for (i=0;;i++) {
								TAKE_LINE_2(hypso0_z[i], hypso0_a[i])
								if (hypso1_a[i]==0) depth0 = hypso0_z[i];
							}
							np0=i;
						}
						break;
					}
					ptr=strtok(NULL, "/");
					if (ptr != NULL && *ptr != 'l') {
						if (*ptr == 'b') {
							/*Default box-like hypsometry*/
							hypso0_z[0]=depth0; hypso0_a[0]=0;
							hypso0_z[1]=depth0; hypso0_a[1]= voltot0/(z0-depth0);
							np0 = 2;
							if (verbose_level>=3) fprintf(stderr, "\nAutomatic box-like hypsometry points for basin0: %d", np0);
						}
					}
					else {
						/*Default linear hypsometry*/
						hypso0_z[0]=depth0; hypso0_a[0]=0;
						hypso0_z[1]=z0; hypso0_a[1]= voltot0/(z0-depth0);
						np0 = 2;
						if (verbose_level>=2) fprintf(stderr, "\nAutomatic linear hypsometry points for basin0: %d", np0);
					}
					break;
				case 's':
					sl_amp = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) sl_per = atof(ptr);
					break;
				case 't':
					timeini = atof(strtok(prm, "/"))*secsperyr;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) timeend = atof(ptr)*secsperyr;
					else timeend=timeini+=500*secsperyr;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dt = atof(ptr)*secsperyr;
					else dt=1*secsperyr;
					break;
				case 'u':
					uplift_rate = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) uplift_incr = atof(ptr);
					else uplift_incr=0;
					break;
				case 'V':
					verbose_level = value;
					break;
				case 'W':
					water_conservative=1;
					break;
				case 'w':
					Kw = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) expw = atof(ptr);
					else {linear_width=1; break;}
					ptr=strtok(NULL, "/");
					if (ptr != NULL && *ptr == 'p') decreasable_width=0;
					if (ptr != NULL && *ptr == 'm') turowski_width=1;
					break;
				case 'z':
					z0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) z1 = atof(ptr);
					else z1=z0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) z2 = atof(ptr);
					else z2=z1;
					break;

				default:
					fprintf(stderr, "\nWarning in %s: incomprehensible parameter '%s'.", argv[0], argv[iarg]);
					break;
			}
		}
	}

	if (verbose_level>=4) for (i=0;i<np0;i++) fprintf(stderr, "\n#basin0  [%d]: z,a = %f m,  %e m2", i, hypso0_z[i], hypso0_a[i]);
	if (verbose_level>=4) for (i=0;i<np1;i++) fprintf(stderr, "\n#basin1  [%d]: z,a = %f m,  %e m2", i, hypso1_z[i], hypso1_a[i]);
	if (verbose_level>=4) for (i=0;i<np2;i++) fprintf(stderr, "\n#basin2  [%d]: z,a = %f m,  %e m2", i, hypso2_z[i], hypso2_a[i]);
	if (verbose_level>=1 && hypso0_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin0 should have area 0 in first (deepest) point: %.2f != 0 m", hypso0_a[0]);
	if (verbose_level>=1 && hypso1_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin1 should have area 0 in first (deepest) point: %.2f != 0 m", hypso1_a[0]);
	if (verbose_level>=1 && hypso2_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin2 should have area 0 in first (deepest) point: %.2f != 0 m", hypso2_a[0]);
	
	if (nbasins==1) {r2=p2=e2=0;}
	
	/*Initial volume and area of basins*/
	volume_and_area_from_level(hypso0_z, hypso0_a, np0, z0, &vol0, &area0);
	volume_and_area_from_level(hypso1_z, hypso1_a, np1, z1, &vol1, &area1);
	volume_and_area_from_level(hypso2_z, hypso2_a, np2, z2, &vol2, &area2);

	if (verbose_level>=1) 
		fprintf(stderr, "\n====== %s's initial parameters ======"
			"\nHypsometry points: basin0= %d;  basin1= %d;  basin2= %d"
			"\nmodel_eros     = %d (0 for shear stress; 1 for V in incision law)"
			"\nmodel_vel      = %d (0 for Manning's law; 1 for Critical flow; 2 for mixed)"
			"\ntimeini,end,dt = %.1f , %.1f, %.2f yr"
			"\nz0,1,2         = %.2f , %.2f , %.2f m"
			"\narea0,1,2      = %.3e , %.3e , %.3e m3"
			"\nvol0,1,2       = %.3e , %.3e , %.3e m3"
			"\nuplift_rate    = %.4f m/y"
			"\nz_sill1,2      = %.2f , %.2f m"
			"\nheadloss1,2    = %.2f , %.2f m"
			"\ndist1,dist2    = %.0f , %.0f m"
			"\nKe,expe        = %.3e m y^-1 Pa^-a, %.2f, %.1f Pa"
			"\nnbasins        = %1d"
			"\ndepth0,1,2     = %.1f , %.1f , %.1f m"
			"\nvoltot0,1,2    = %.3e , %.3e , %.3e m3"
			"\ne0,e1,e2       = %.3f , %.3f , %.3f m/y"
			"\np0,p1,p2       = %.3f , %.3f , %.3f m/y"
			"\nr0,r1,r2       = %.2f , %.2f , %.2f m3/s"
			"\nKw,expw        = %.2f , %.2f"
			"\nroughness      = %.2f"
			"\nmixing 01,12   = %.2f , %.2f ratio/yr; sill depth linera mix (max,min) = %.1f, %.1f"
			"\n", 
			argv[0], 
			np0, np1, np2, 
			model_eros, model_eros, timeini/secsperyr, timeend/secsperyr, dt/secsperyr, 
			z0, z1, z2, 
			area0, area1, area2, 
			vol0, vol1, vol2, 
			uplift_rate, 
			z_sill1, z_sill2, hl1, hl2, dist1, dist2, 
			Ke, expe, shearcritical, 
			nbasins, depth0, depth1, depth2, voltot0, voltot1, voltot2, 
			e0, e1, e2, p0, p1, p2, r0, r1, r2, Kw, expw, roughness, 
			mix01_ratioperyr, mix12_ratioperyr, silldepthmixmax, silldepthmixmin); 
	
	time=timeini; 

	/*Initial Gypsum and Halite dissolved in each water body*/
	gypsumds0 = gypsumcn0*vol0*1e3;
	gypsumds1 = gypsumcn1*vol1*1e3;
	gypsumds2 = gypsumcn2*vol2*1e3;
	haliteds0 = halitecn0*vol0*1e3;
	haliteds1 = halitecn1*vol1*1e3;
	haliteds2 = halitecn2*vol2*1e3;

	fprintf (stderr, "\nStart calculation in %s\n", argv[0]);
	fprintf (stdout, CAPTION); 


	do {
		float areasill1, areasill2, Rh1, Rh2, slope1, slope2, disch1, disch2, 
			shearstr1, shearstr2, erosrate1, erosrate2, dvoltr1, dvoltr2, 
			Dsill1, Dsill2, width1, width2, vel1, vel2; 
		float dhaliteds0, dhaliteds1, dhaliteds2, dgypsumds0, dgypsumds1, dgypsumds2, 
		z0_aux, z1_aux, z2_aux, area0_aux, area1_aux, area2_aux, mix_vol, mix01_interp, 
		vel_Manning1, vel_Manning2, vel_critical1, vel_critical2, factr=2;

		/*Cross-sectional area of the sill gate. Note z_sill1 is the MEAN depth of the sill*/
		Dsill1 = MAX_2(z0-z_sill1,0);
		Dsill2 = MAX_2(z1-z_sill2,0);
		if (time==timeini) width1 = 5*Dsill1; /*Suggested by Nature's Rev#3, 2011*/
		width2 = 5*Dsill2; 
		areasill1 = width1*Dsill1;
		areasill2 = width2*Dsill2;

		/*Hydraulic radius*/
		if (Dsill1>0) Rh1 = areasill1/(2*sqrt(4*Dsill1*Dsill1+width1*width1/4)); else Rh1=0;
		if (Dsill2>0) Rh2 = areasill2/(2*sqrt(4*Dsill2*Dsill2+width2*width2/4)); else Rh2=0;
		slope1 = (z1<z0)? (MAX_2(z1-z0,hl1))/dist1 : 0;
		slope2 = (z2<z1)? (MAX_2(z2-z1,hl2))/dist2 : 0;
		/*fprintf(stderr,"\n>>>>>>>>> %f  %f  %f  %f",slope, z1-z0, hl1, dist1);*/
		/*Manning formula or critical flow*/
		vel_Manning1 = 1/roughness*pow(Rh1,(float)2/3)*pow(-slope1,.5); 
		vel_Manning2 = 1/roughness*pow(Rh2,(float)2/3)*pow(-slope2,.5); 
		vel_critical1 = sqrt(g*Dsill1);
		vel_critical2 = sqrt(g*Dsill2);
		//fprintf(stderr, "\n%.2f %.2f", vel_Manning1, vel_critical1);
		switch (model_vel) {
		    case 0:
			/*Manning's law*/
			vel1 = vel_Manning1;
			vel2 = vel_Manning2;
			break;
		    case 1:
			/*Critical flow*/
			vel1 = vel_critical1;
			vel2 = vel_critical2;
			break;
		    case 2:
			/*Linear transition between both laws, from z1,2=z_sill1,2 / 2 to z1,2=z_sill1,2 * 2   */
			if (z1>z_sill1/factr) vel1 = vel_Manning1;
			if (z1<z_sill1*factr) vel1 = vel_critical1;
			if (z1<=z_sill1/factr && z1>=z_sill1*factr) vel1 = vel_Manning1 + (z_sill1/factr - z1)*(vel_critical1 - vel_Manning1)/(z_sill1/factr - z_sill1*factr);
			if (z2>z_sill2/factr) vel2 = vel_Manning2;
			if (z2<z_sill2*factr) vel2 = vel_critical2;
			if (z2<=z_sill2/factr && z2>=z_sill2*factr) vel2 = vel_Manning2 + (z_sill2/factr - z2)*(vel_critical2 - vel_Manning2)/(z_sill2/factr - z_sill2*factr);
			break;
		    default:
			fprintf(stderr, "\nWarning in %s: incomprehensible velocity model '%d'.", argv[0], model_vel);
			break;
		}
		disch1 = vel1 * areasill1; 
		disch2 = vel2 * areasill2; 
		shearstr1 = denswater*g*Dsill1*(-slope1); 
		//shearstr1 = 0.01 * vel1*vel1; /*e.g., Paola & Heller, 1992*/
		maxshearstr1 = MAX_2(maxshearstr1, shearstr1);
		shearstr2 = denswater*g*Dsill2*(-slope2); 
		//shearstr2 = 0.01 * vel2*vel2; /*e.g., Paola & Heller, 1992*/
		maxshearstr2 = MAX_2(maxshearstr2, shearstr2);
		//fprintf(stderr,"\nshear stress (sill1, sill2) %.2f Pa", shearstr1, shearstr2);
		if (shearstr1>shearcritical) erosrate1 = Ke*pow((shearstr1-shearcritical)*((model_eros)?vel1:1), expe); else erosrate1 =0;
//		if (shearstr2>shearcritical) erosrate2 = Ke*pow((shearstr2-shearcritical)*((model_eros)?vel2:1), expe); else erosrate2 =0;
		/*Water volume crossing sills during dt*/
		dvoltr1 = MIN_2(disch1*dt, vol0);
		dvoltr2 = MIN_2(disch2*dt, vol1); if (nbasins==1) {vel2=disch2=erosrate2=dvoltr2=0;}

		/*Limit dvoltr to avoid z2>z1 or z1>z0*/
		level_and_area_from_volume(hypso0_z, hypso0_a, np0, vol0-dvoltr1, &z0_aux, &area0_aux); 
		level_and_area_from_volume(hypso1_z, hypso1_a, np1, vol1+dvoltr1, &z1_aux, &area1_aux);
		if (z0>z_sill1 && z0_aux<z1_aux) {
			//fprintf(stderr,"\nLIMITING dvoltr1 : %.2f yr  %.2f  %.2f m\tvol=%.3e red.factor=%.3e", time/secsperyr, z1, z2, dvoltr1, .999*(z0-z1)/((z0-z1)-(z0_aux-z1_aux)));
			dvoltr1 *= .999*(z0-z1)/((z0-z1)-(z0_aux-z1_aux)); 
			disch1 = dvoltr1 / dt;
			vel1 = disch1 / areasill1;
		}
		//fprintf(stderr,"\n>>>>>>>>> %.2f yr  %.2f  %.2f     %.2f  %.2f  \t %.3e %.3e", time/secsperyr, z0, z1, z0_aux, z1_aux, dvoltr1, dvoltr1*.5*(z0-z1)/((z0-z1)-(z0_aux-z1_aux)));
		level_and_area_from_volume(hypso1_z, hypso1_a, np1, vol1-dvoltr2, &z1_aux, &area1_aux); 
		level_and_area_from_volume(hypso2_z, hypso2_a, np2, vol2+dvoltr2, &z2_aux, &area2_aux);
		if (z1>z_sill2 && z1_aux<z2_aux) {
			//fprintf(stderr,"\nLIMITING dvoltr2 : %.2f yr  %.2f  %.2f m\tvol=%.3e red.factor=%.3e", time/secsperyr, z1, z2, dvoltr2, .999*(z1-z2)/((z1-z2)-(z1_aux-z2_aux)));
			dvoltr2 *= .999*(z1-z2)/((z1-z2)-(z1_aux-z2_aux)); 
			disch2 = dvoltr2 / dt; 
			vel2 = disch2 / areasill2; 
		}
		//fprintf(stderr,"\n######### %.2f yr  %.2f  %.2f     %.2f  %.2f  \t %.3e %.3e", time/secsperyr, z1, z2, z1_aux, z2_aux, dvoltr2, dvoltr2*.5*(z1-z2)/((z1-z2)-(z1_aux-z2_aux)));

		/*TRANSFER WATER BETWEEN BASINS to obtain new volumes*/
		if (vol0) {
			vol0 -= dvoltr1; vol0=MAX_2(vol0,0); 
			gypsumds0 -= dvoltr1*gypsumcn0*1e3; gypsumds0 = MAX_2(gypsumds0,0); 
			haliteds0 -= dvoltr1*halitecn0*1e3; haliteds0 = MAX_2(haliteds0,0); 
		}
		else {
			dvoltr1 = 0;
		}
		if (z1<=z_sill2 || nbasins==1) {
			/*Keep water and salt in basin1*/
			vol1 += dvoltr1;
			gypsumds1 += dvoltr1*gypsumcn0*1e3; 
			haliteds1 += dvoltr1*halitecn0*1e3; 
		}
		else {
			/*Transfer water and salt to both basins*/
			vol1 += dvoltr1-dvoltr2;
			vol2 += dvoltr2;
			gypsumds1 += dvoltr1*gypsumcn0*1e3 - dvoltr2*gypsumcn1*1e3; 
			gypsumds2 += dvoltr2*gypsumcn1*1e3; 
			if (dvoltr2<0 || gypsumcn1<0 || gypsumds2<0) fprintf(stderr,"\nERROR in %s: time=%.2f gypsumds2=%.3e", argv[0], time/secsperyr, gypsumds2);
			haliteds1 += dvoltr1*halitecn0*1e3 - dvoltr2*halitecn1*1e3; 
			haliteds2 += dvoltr2*halitecn1*1e3; 
		}
		/*Find level and area by filling hypsometry with volume*/
		level_and_area_from_volume(hypso0_z, hypso0_a, np0, vol0, &z0, &area0);
		level_and_area_from_volume(hypso1_z, hypso1_a, np1, vol1, &z1, &area1);
		level_and_area_from_volume(hypso2_z, hypso2_a, np2, vol2, &z2, &area2);
		if (nbasins==1) {z2=z1; area2=area1; vol2=vol1;}
	


		/*Sea level change in z0*/
		if (sl_per) {
			float sl_change;
			sl_change = sl_amp*(sin(2*PI*(time-timeini)/secsperyr/sl_per) - sin(2*PI*(time-timeini-dt)/secsperyr/sl_per));
			vol0 += sl_change*area0;
		}

		voltr += dvoltr1;

		/*VARIABLE E,P*/
		if (n_insolation_input_points) {
			float insolation;
			if (time<=insolation_var[0][0] || time>=insolation_var[n_insolation_input_points-1][0]) {
				if (time<=insolation_var[0][0]) 
					{insolation = insolation_var[0][1];}
				if (time>=insolation_var[n_insolation_input_points-1][0]) 
					{insolation = insolation_var[n_insolation_input_points-1][1];}
			}
			else for (i=0; i<n_insolation_input_points-1; i++) {
				if (time>insolation_var[i][0] && time<=insolation_var[i+1][0]) {
					insolation = insolation_var[i][1]+(time-insolation_var[i][0])*(insolation_var[i+1][1]-insolation_var[i][1])/(insolation_var[i+1][0]-insolation_var[i][0]); 
					break;
				}
			}
			EPfac = insolation/insolation_mean;
			//fprintf(stdout, " # %.2f\t%.2f", insolation, EPfac);
		}
		else {
			EPfac = 1;
		}

		/*ADD RUNOFF+PRECIPITATION-EVAPORATION TO BASINS*/
		vol0 += EPfac*r0*dt + (EPfac*p0-e0)/secsperyr*area0*dt; vol0=MAX_2(vol0,0); 
		vol1 += EPfac*r1*dt + (EPfac*p1-e1)/secsperyr*area1*dt; vol1=MAX_2(vol1,0); 
		vol2 += EPfac*r2*dt + (EPfac*p2-e2)/secsperyr*area2*dt; vol2=MAX_2(vol2,0);
		if (water_conservative) {vol0 -= r1*dt + (p1-e1)/secsperyr*area1*dt + r2*dt + (p2-e2)/secsperyr*area2*dt;}

		/*ADD SALT MASS FROM RIVERS TO BUDGET OF DISSOLVED SALT*/
		gypsumds0 += GYPSUMCNRIVERS*r0*dt;
		gypsumds1 += GYPSUMCNRIVERS*r1*dt;
		gypsumds2 += GYPSUMCNRIVERS*r2*dt;
		haliteds0 += HALITECNRIVERS*r0*dt;
		haliteds1 += HALITECNRIVERS*r1*dt;
		haliteds2 += HALITECNRIVERS*r2*dt;

		/*MIX basin0 and basin1*/
		if (silldepthmixmax) {
			mix01_interp = mix01_ratioperyr * (silldepthmixmin-z_sill1)/(silldepthmixmin-silldepthmixmax); mix01_interp=MAX_2(0,mix01_interp);
		//fprintf(stderr, "\ntime= %.2f yr; mix01_interp ratio = %.2e", time/secsperyr, mix01_interp);
		}
		else {
			mix01_interp = mix01_ratioperyr;
		}
		mix_vol = mix01_interp/secsperyr * MIN_2(vol0,vol1)*dt; 
		gypsumds0 -=   mix_vol*(gypsumcn0-gypsumcn1)*1e3;
		gypsumds1 +=   mix_vol*(gypsumcn0-gypsumcn1)*1e3;
		haliteds0 -=   mix_vol*(halitecn0-halitecn1)*1e3;
		haliteds1 +=   mix_vol*(halitecn0-halitecn1)*1e3;

		/*MIX basin1 and basin2*/
		mix_vol = mix12_ratioperyr/secsperyr * MIN_2(vol1,vol2)*dt; 
		gypsumds1 -=   mix_vol*(gypsumcn1-gypsumcn2)*1e3;
		gypsumds2 +=   mix_vol*(gypsumcn1-gypsumcn2)*1e3;
		haliteds1 -=   mix_vol*(halitecn1-halitecn2)*1e3;
		haliteds2 +=   mix_vol*(halitecn1-halitecn2)*1e3;

		/*SALT CONCENTRATION*/
		gypsumcn0 = gypsumds0/vol0/1e3;
		gypsumcn1 = gypsumds1/vol1/1e3;
		gypsumcn2 = gypsumds2/vol2/1e3;
		halitecn0 = haliteds0/vol0/1e3;
		halitecn1 = haliteds1/vol1/1e3;
		halitecn2 = haliteds2/vol2/1e3;

		/*SALT PRECIPITATION*/
		if (gypsumcn0+halitecn0>GYPSUMPRECIPCN) {dgypsumds0 = MIN_2((gypsumcn0+halitecn0-GYPSUMPRECIPCN)*vol0*1e3, gypsumds0); gypsumpr0 += dgypsumds0; gypsumds0 -= dgypsumds0; }
		if (gypsumcn1+halitecn1>GYPSUMPRECIPCN) {dgypsumds1 = MIN_2((gypsumcn1+halitecn1-GYPSUMPRECIPCN)*vol1*1e3, gypsumds1); gypsumpr1 += dgypsumds1; gypsumds1 -= dgypsumds1; }
		if (gypsumcn2+halitecn2>GYPSUMPRECIPCN) {dgypsumds2 = MIN_2((gypsumcn2+halitecn2-GYPSUMPRECIPCN)*vol2*1e3, gypsumds2); gypsumpr2 += dgypsumds2; gypsumds2 -= dgypsumds2; }
		if (gypsumcn0+halitecn0>HALITEPRECIPCN) {dhaliteds0 = MIN_2((gypsumcn0+halitecn0-HALITEPRECIPCN)*vol0*1e3, haliteds0); halitepr0 += dhaliteds0; haliteds0 -= dhaliteds0; }
		if (gypsumcn1+halitecn1>HALITEPRECIPCN) {dhaliteds1 = MIN_2((gypsumcn1+halitecn1-HALITEPRECIPCN)*vol1*1e3, haliteds1); halitepr1 += dhaliteds1; haliteds1 -= dhaliteds1; }
		if (gypsumcn2+halitecn2>HALITEPRECIPCN) {dhaliteds2 = MIN_2((gypsumcn2+halitecn2-HALITEPRECIPCN)*vol2*1e3, haliteds2); halitepr2 += dhaliteds2; haliteds2 -= dhaliteds2; }

		fprintf (stdout, 
			"\n%6.2f\t%6.2f\t%6.1f\t%7.3e\t%6.2f\t%.2e\t%6.1f\t%7.5f"
			"\t%.2e\t%.2e\t%.2e\t%.2e\t%7.2f\t%6.2f\t%6.2f"
			"\t\t%6.1f\t%6.1f\t%6.1f\t%5.4f\t%5.4f\t%5.4f"
			"\t\t%6.1f\t%6.1f\t%6.1f\t%5.4f\t%5.4f\t%5.4f", 
			time/secsperyr, z_sill1, Rh1, slope1, vel1, disch1, width1, erosrate1, 
			vol0/1e9, vol1/1e9, vol2/1e9, voltr/1e9, z0, z1, z2, 
			gypsumpr0/1e12, gypsumpr1/1e12, gypsumpr2/1e12, gypsumcn0, gypsumcn1, gypsumcn2,
			halitepr0/1e12, halitepr1/1e12, halitepr2/1e12, halitecn0, halitecn1, halitecn2);

		if (linear_width) {
			width1 = MAX_2(-Kw*z_sill1, 0); 
			//fprintf(stderr, "\n?????? \t\t\t\t\tKw=%.2e  depth1=%.2e \twidth1=%e", Kw, depth1, width1);
		}
		else {
			float widthlaw;
			if (turowski_width) { /*Turowsky et al., 2007, 2009*/
				widthlaw = Kw*pow((shearcritical+uplift_rate/Ke)/denswater/g, -3./13)*pow(fabs(roughness*disch1), 6./13); 
				//fprintf(stderr, "\n\t\twidthlaw=%e, %e, %e, %e, %e, %e", widthlaw, (shearcritical+uplift_rate/Ke)/denswater/g, roughness*disch1, Kw, pow((shearcritical+uplift_rate/Ke)/denswater/g, -3./13), pow(fabs(roughness*disch1), 6./13));
				//if (expe!=1) fprintf(stderr, "\nError: Turowski's width only valid for erosion law exponent a=expe=1, but a=%.1f", expe);
			}
			else {
				widthlaw = Kw*pow(fabs(disch1), expw); 
			}
			widthlaw = MAX_2(widthlaw, Dsill1);
			if (decreasable_width) width1 = widthlaw; else width1 = MAX_2(width1, widthlaw);
		}
		z_sill1 += -erosrate1*dt/secsperyr + uplift_rate*dt/secsperyr;
//		z_sill2 += -erosrate2*dt/secsperyr;
		uplift_rate += uplift_incr*dt/1000/secsperyr;
		erostotal += erosrate1*dt/secsperyr;
		//z_sill1 = MAX_2(z_sill1, hypso0_z[0]);
		time += dt;
		if (timeend<timeini && fabs(z1-z0)<.1) timeend=time+(time-timeini)*.05;
	} while (time <= timeend || timeend<timeini);

	fprintf (stdout, "\n"CAPTION); 
	fprintf (stderr, "\nIncision total = %.1f m", erostotal);
	fprintf (stderr, "\nMax. shear stress (sill1, sill2)  %.2f , %.2f Pa", maxshearstr1, maxshearstr2);

	if (switch_ps) system("spillover.gmt.csh spillover.result");

	fprintf (stderr, "\nExiting %s\n", argv[0]);

}


int level_and_area_from_volume (float *hypso_z, float *hypso_a, int np, float vol, float *z, float *area) {
	int i;
	float vol_aux=0; 
	if (vol<0) fprintf(stderr, "\nERROR in level_and_area_from_volume: negative volume passed: %.2e m3", vol);
	for (i=1; i<np; i++) {
		float dz, dvol;
		dvol = (hypso_z[i]-hypso_z[i-1])*(hypso_a[i]+hypso_a[i-1])/2;
		//fprintf(stderr, "\n************** %.2f %.2e   %.2f  %e     %d, %d     %.3e  %.3e", hypso_z[i], hypso_a[i], *z, *area, i, np, vol, vol_aux+dvol);
		if (vol_aux+dvol >= vol) {
			double a = (hypso_a[i]-hypso_a[i-1])/(hypso_z[i]-hypso_z[i-1])/2, b=hypso_a[i-1], c=-(vol-vol_aux);
			if (a==0) {
				dz = -c/b;
			}
			else {
				if (b==0)
					dz = sqrt(-c/a);
				else
					dz = (-b+sqrt(b*b-4*a*c))/(2*a);
			}
			*z = hypso_z[i-1] + dz;
			*area = b+dz*a;
			break;
		}
		vol_aux += dvol;
	}
	/*Assumes constant area above last point in hypso_z*/
	if (i==np) {
		*z = hypso_z[np-1] + (vol-vol_aux)/hypso_a[np-1];
		*area = hypso_a[np-1];
	}
//fprintf(stderr, "\n!!!!!!!!!!!!!! %.2f %.2e   %.2f  %e", hypso_z[i], hypso_a[i], *z, *area);
}


int volume_and_area_from_level (float *hypso_z, float *hypso_a, int np, float z, float *vol, float *area) {
	int i;
	*vol=0;
	if (z<hypso_z[0]) fprintf(stderr, "\nERROR in volume_and_area_from_level: passed level below basin floor: %.2f m", z);
	for (i=1; i<np; i++) {
		if (hypso_z[i]<=z) {
			*area = hypso_a[i];
			*vol    += (hypso_z[i]-hypso_z[i-1])*(hypso_a[i]+hypso_a[i-1])/2; 
		}
		else {
			*area = hypso_a[i-1]+((hypso_a[i]-hypso_a[i-1])/(hypso_z[i]-hypso_z[i-1]))*(z-hypso_z[i-1]);
			*vol    += (z-hypso_z[i-1])*(*area+hypso_a[i-1])/2; 
			break;
		}
	}
}




int read_file_insolation(char *filename, int *n_insolation_input_points, float *insolation_mean) {
	/*
	  Reads file with insolation in W/m2
	*/

	int 	i, j, nmax_input_points=50000;
	FILE 	*file;
	float	*aux1, *aux2, insomin=1e12, insomax=-1e12;

	if ((file = fopen(filename,"r")) == NULL) {if (verbose_level>=2) fprintf(stderr, "\nInfo: Cannot read insolation file input file '%s'.", filename); return 0;}
	if (verbose_level>=2) fprintf(stderr, "\nInfo: Reading insolation at '%s'", filename);

	(*n_insolation_input_points)=0;
	aux1 = calloc(nmax_input_points, sizeof(float));
	aux2 = calloc(nmax_input_points, sizeof(float));
	
	for (;;) {
		TAKE_LINE_2(aux1[*n_insolation_input_points], aux2[*n_insolation_input_points]);
		//fprintf(stderr, "\n>>>>>> %f  %f  %d", aux1[*n_insolation_input_points], aux2[*n_insolation_input_points], *n_insolation_input_points);
		(*n_insolation_input_points)++;
		if ((*n_insolation_input_points)>=nmax_input_points-1 ) {
			fprintf(stderr, "\nERROR: Too many points (%d) in insolation file.", *n_insolation_input_points);
			break;
		}
	}
	fclose(file); 
	insolation_var = calloc((*n_insolation_input_points), sizeof(float *));
	for (i=0; i<(*n_insolation_input_points); i++) {
		insolation_var[i] = calloc(2, sizeof(float));
		insolation_var[i][0] = aux1[i]*secsperyr;
		insolation_var[i][1] = aux2[i];
		(*insolation_mean) += insolation_var[i][1];
		insomin = MIN_2(insomin, insolation_var[i][1]);
		insomax = MAX_2(insomax, insolation_var[i][1]);
	}
	(*insolation_mean) /= (*n_insolation_input_points);
	if (verbose_level>=2) 
		fprintf(stderr, "\nInfo: %d insolation points in '%s'. "
			"\nMin,Max:  %.2f, %.2f ; Mean: %.2f W/m2."
			"\nFirst time: %.2f yr", 
			*n_insolation_input_points, filename, insomin, insomax, *insolation_mean, insolation_var[0][0]/secsperyr); 
	free(aux1); free(aux2);
	return(1);
}



int syntax (int argc, char **argv) {
	fprintf(stderr, "\n"
			"\nSyntax: \t %s  -S[<voltot0>/<depth0>/<b|l>|<hypsfile0>] "
			"\n\t-B<z_sill1>/<hl1>/<dist1>[/<voltot1>/<depth1>/<b|l>|<hypsfile1>] "
			"\n\t-b<z_sill2>/<hl2>/<dist2>[/<voltot2>/<depth2>/<b|l>|<hypsfile2>] "
			"\n\t-e<e0>/<e1>/<e2> [-h] -i<insolationfile> -k<Ke>/<expe>/<critshear> "
			"\n\t-M<model_eros>[/<model_vel>] -m<mix01>/<mix12>/<depthmaxmix>/<depthminmix> "
			"\n\t-p<p0>/<p1>/<p2> -R<roughness> -r<r0>/<r1>/<r2> -s<amp>/<per> "
			"\n\t-t<timeini>/<timeend>/<dt> -u<uplift_rate>[/<accel>] [-V<level>] -W "
			"\n\t-w<Kw>/<expw>[/d] -z<z0>/<z1>/<z2> ", argv[0]);
	fprintf(stderr, "\n"
			"\nThis program calculates the water transfer from an ocean (basin0) overspilling into "
			"a receiving basin (basin1) and if needed also between basin1 and a basin2. "
			"Water circulates only in that sense 0->1->2. The incision produced at sill1 "
			"(between basin 1 and basin2) and its coupling with water flow are also calculated. "
			"\n\nNotation: "
			"\nz is elevation (negative downwards) and measured in [m]. "
			"\nbasin0, basin1, and basin 2 refer to the source basin, the upper receiving basin (first flooded by basin0) and the lower basin (to be flooded from basin1). "
			"\n"
			"\n\t-S if parameters are passed, they indicate that the source basin is not infinite (default) but it has either a linear (l) or box-like (b) hypsometry defined by <voltot0> and <depth0> or by the hypsfile0 file."
			"\n\t-B gives sill elevation between basin0 and basin1 (must be below z0), max.headloss [m] (<0), slope distance [m], and, if not given in the hypsometry file, volume [m3], max. depth of basin1 [m], and 'b' or 'l' for box or linear hypsometry. "
			"\n\t-b to account for a susidiary basin2. Here sill2 refers to the sill connecting basin1 and basin2. This sill is NOT eroded as sill1. "
			"\n"
			"\n\t-e to give the evaporation at water surface [m/y]."
			"\n\t-h for help."
			"\n\t-i to supply insolation changes through time (Milankovitch). Changes in insolation (2nd column) relative to the mean value will be used to scale P-E in time."
			"\n\t-k to specify the three constants of the incision law e=Ke*(shear-critshear)^expe. Ke is in units [m/y/(Pa)^a]. critshear is in [Pa]"
			"\n\t-M to specify the erosion <model_eros> (0 for shear stress; 1 for stream power per unit area with V explicit in incision law) and the flow velocity model (0 for Manning's eq.; 1 for critical flow; 2 for a linear transition between both according to z1 level-zsill)."
			"\n\t-m to specify the ammount of mix between basins 0 and 1, and 1 and 2 (ratio of volume per year). If <depthmaxmix> and <depthminmix> are given, then mix is linearly interpolated between <mix01> and 0 between those depths; otherwise constant. depthmaxmix<depthminmix<0. "
			"\n\t-P to call the gmt script and produce a postscript with graphics. "
			"\n\t-p to give the precipitation fallen on basins [m/y]. Takes largest area in hypsometry. "
			"\n\t-R to change the default coefficient of roughness [adimensional]."
			"\n\t-r to give the runoff collected by rivers to basins [m3/s]."
			"\n\t-s impose sea level variations on z0 of amplitude <amp> [m] and period <per> [yr]."
			"\n\t-t sets initial and final run time [yr]. <dt> is the time step. A negative timeend will find an automatic value. "
			"\n\t-u sets the sill uplift rate [m/y] and its acceleration (m/y2)). Current sea level fall is 3mm/yr. Last postglacial maximum was 12 mm/yr."
			"\n\t-W to drop the water evaporated in basins 1 and 2 at basin 0 (water conservative)."
			"\n\t-w with one argument specifies the fixed width; otherwise specifies two constants of the width law W=Kw*Q^expw [-]. Kw is between 1 and 10, expw is around 0.5. Add '/p' to impede width to decrease. Add '/m' to use Turowski et al. (2007,2009) formula (ignores exponent value; uses 6/13)."
			"\n\t-z to specify the water level in each basin [default is 0 for basin0 and basin bottom for basin 1 and basin 2]."
			"\n\nIncision law: e=Ke*tau^expe (for model=0) ; e=Ke*(tau*V)^expe (for model=1), where V is the water velocity. Note that for model=1, expe is equivalent to n (exponent of slope) in the stream power law of river incision. For model=0, expe=3/2*n. So, according to Whipple & Tucker (1999) expe=[2/3-2] for model=1 and expe=[1-3] for model=0."
			"\nReads hypsometry from files <hypsfile0> (source basin), <hypsfile1> (basin1), <hypsfile2> (basin2), in two columns (z, basin_area_above_z). z runs from basin bottom to the higghest z0 [m]. First area [m2] is thus 0. "
			"\nSource basin has an infinite hypsometry by default (ocean)."
			"\nWrites in stdout: "
			CAPTION
			"\nThese are: time [yr], sill channel depth [m], hydraulic radius [m], slope, "
			"\nwater velocity, discharge, channel width, incision rate, transferred volume, basin level");
	fprintf(stderr, "\n"
			"\nExamples: "
			"\nMediterranean desiccation, two basins:"
			"\nasalted -s0/1000 -m.014/.002/-284/-30 -u0.0049/.000018 -k1e-6/1.5/50 -W -V3 -r0/4.5e3/12e3 -p0/.6 -e0/1.1/1.3 -w6/0.5/m "
			"\t-M0/2 -Shypso0.xa.tmp -B-100/-1000/100000/hypso1.xa.tmp -b-430/-1000/100000/hypso2.xa.tmp "
			"\t-t0/120000/.2 -z0/-1/-2\n"
			"\nMediterranean flood, two basins (data from Blanc, 2002):"
			"\nasalted -B-10/-1000/100000/1.5237e15/-1509/l -b-430/-1000/100000/1.247e15/-1718/l -t0/900/1\n"
			"\nMediterranean flood, one basin (data from Blanc, 2002):"
			"\nasalted -B-10/-1000/100000/3.618e15/-1509/l -t0/900/1\n"
			"\n"
			"\nFollows a list of default parameters and the first time steps using those, as an example run:"
			"\n");
	AUTHORSHIP;
	return(0);
}


