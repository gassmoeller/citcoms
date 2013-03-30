#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include </usr/include/sys/param.h>


void perplex_data(float ref_alpha,float* Tadi,float ref_density);
void output_vtk(float*** c_morb,float*** c_harz);

    const int densityfile = 1;
    const int nodes = 65;
    int ntemp = 424;
    const int disptn_number = 2.0;
    const int gruneisen = 2.667;
    const int i_mantle_depth = 2890;
    const float mantle_depth_km = 2890.0;
    const float surf_Tadi = 1613.0;
    const float g = 9.81;
    const float mantle_depth_m = 2890.0 * 1000.0;
    const float ln_Tsubadi = 5.7; // 300 K
    const float R_km = 6371.0;
    const float R = 8.314;
    const float Delta_Ks = 0.06;
    const float Delta_Rho_CMB = 0.015;
    const float therm_diff_CMB = 1.5;
    const float refvisc = 1.0e21;
    const float ref_temperature = 2040.0;
    const float ref_cp = 1200.0;
    const float ref_density = 3340;    
    const int density_profile = 1; // 1 for Eh's profile, 2 for Juliane's
    const int mintemp = 270;
    const int deltatemp = 10;
    FILE *ifile,*ofile;
    char *bufline, bufname1[200],bufname2[200];
    size_t length = 30;
    int i,k;
    float b_alpha[2890];
    float b_H[2890];
    float b_Viscosity[22];
    float b_visc_depth[22];
    float j_rho1[257];
    float j_rho2[257];
    float j_rho_depth[257];
    float temp[2];
    float depth_km;
    int i_depth_km;
        double deltaz, alpha;
	double slope;
    float* c_alpha;
    float* Tadi;
    float* Tfinal;
    float* c_visc;
    float* c_stress_exp;
    float* c_H;
    float* c_visc_calc;
    float* c_delta_rho1;
    float* c_delta_rho2;
    float* c_rho1;
    float* c_rho2;
    float* c_therm_diff;


int main()
{
    float* c_alpha = malloc(sizeof(float) * nodes);
    float* Tadi = malloc(sizeof(float) * nodes);
    float* Tfinal = malloc(sizeof(float) * nodes);
    float* c_visc = malloc(sizeof(float) * nodes);
    float* c_stress_exp = malloc(sizeof(float) * nodes);
    float* c_H = malloc(sizeof(float) * nodes);
    float* c_visc_calc = malloc(sizeof(float) * nodes);
    float* c_delta_rho1 = malloc(sizeof(float) * nodes);
    float* c_delta_rho2 = malloc(sizeof(float) * nodes);
    float* c_rho1 = malloc(sizeof(float) * nodes);
    float* c_rho2 = malloc(sizeof(float) * nodes);
    float* c_therm_diff = malloc(sizeof(float)*nodes);



   bufline = (char*) malloc (80*sizeof(char));

sprintf(bufname1,"alpha2.d");
        ifile = fopen(bufname1,"r");

sprintf(bufname2,"refstate.dat");
        ofile = fopen(bufname2,"w");


// Einlesen
                for(k=0;k<i_mantle_depth;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%e   %e    %f \n", &b_alpha[k],&temp[0],&temp[1]);
                }

fclose(ifile);

sprintf(bufname1,"rent.d");
        ifile = fopen(bufname1,"r");

                for(k=0;k<670;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%e   %f \n", &b_H[k],&temp[0]);}

fclose(ifile);

sprintf(bufname1,"density_juliane");
        ifile = fopen(bufname1,"r");

                for(k=0;k<257;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%f %e %e %e\n", &temp[0],&j_rho_depth[k],&j_rho1[k],&j_rho2[k]);
                    }
fclose(ifile);


sprintf(bufname1,"tmelt.d");
        ifile = fopen(bufname1,"r");

                for(k=0;k<2221;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%e   %f \n", &b_H[k+670],&temp[0]);
                    b_H[k+670]*=100.0;}

fclose(ifile);

sprintf(bufname1,"rad_viscosity");
        ifile = fopen(bufname1,"r");

                for(k=0;k<22;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%e   %f \n", &b_Viscosity[k],&b_visc_depth[k]);}


   for (k=0;k<nodes;k++){
       depth_km = ((float) k)/((float)nodes-1.0)*mantle_depth_km;
       if (density_profile == 1) c_rho1[k] = 1.0 / (1.0-disptn_number*(depth_km/R_km)/gruneisen);
       if (density_profile == 2) {
         for (i=0;i<257;i++) if (mantle_depth_km-depth_km <= j_rho_depth[i]*1.007) {c_rho1[k] = j_rho1[i]; c_rho2[k] = j_rho2[i]; break;}
      } 
}

   for (k=0;k<nodes;k++){
       depth_km = ((float) k)/((float)nodes-1.0)*(mantle_depth_km-1);
       i_depth_km = (int) depth_km;
// Thermal Expansivity an Knoten festlegen
       c_alpha[k]=b_alpha[i_depth_km];
//printf("%d i_depth\n",i_depth_km);
// Free Enthalpy festlegen
       c_H[k]=b_H[i_depth_km];

// Stess Exponent festlegen
       c_stress_exp[k] = 1.0;
       if (depth_km < 750) c_stress_exp[k] = 2.0;
       if (depth_km < 650) c_stress_exp[k] = 3.0;

//adiabatische Temperatur an Knoten berechnen
       if (k==0) Tadi[k] = surf_Tadi;
       if (k>0){
        Tadi[k] = Tadi[k-1] * (1+0.5*(c_alpha[k-1]+c_alpha[k])*g*mantle_depth_m/(((float) nodes-1)*1250.0));
        for(i=1;i<11;i++) Tadi[k] = Tadi[k-1] + 0.5*(Tadi[k-1]+Tadi[k]) * 0.5*(c_alpha[k-1]+c_alpha[k])*g*mantle_depth_m/(((float) nodes-1)*1250.0);
       }
       Tfinal[k] = Tadi[k];
       //Tfinal[k] = Tadi[k] - exp(ln_Tsubadi * depth_km/mantle_depth_km);
       //if (depth < 100.0) {Tfinal[k] = Tfinal[k] - 1340 * erfc(depth/100.0);}
       //else if (depth > 2690.0) {Tfinal[k] = Tfinal[k] + 1700 * erfc((2890.0-depth)/200.0);}


//Tiefenabhängige Viskosität
       for (i=0;i<21;i++) if (1.0-depth_km/R_km > b_visc_depth[i]) c_visc[k] = b_Viscosity[i];
       c_visc_calc[k] = exp(c_H[k] / (c_stress_exp[k] * R * Tadi[k]));
       if (density_profile == 1) c_rho2[k] = 1.0 / (1.0 / (c_rho1[nodes-1] + Delta_Rho_CMB) + disptn_number * (0.45 - depth_km/R_km) / (gruneisen * (1 + Delta_Ks)));
       
       
       c_delta_rho1[k] = (c_rho2[k] - c_rho1[k]);
       c_delta_rho2[k] = c_delta_rho1[k];
       c_therm_diff[k] = 1.0 + (therm_diff_CMB - 1.0) * depth_km/mantle_depth_km; 
//       printf("delta_rho1: %f\n",c_delta_rho1[k]);
   }



//ausgeben

if (densityfile==0) for(k=nodes-1;k>=0;k--) fprintf(ofile,"%f %f %e %f %.0f %e %f %f %f %f %f\n",1.0,1.0,c_alpha[k]/c_alpha[0],1.0,Tfinal[k],c_H[k],c_visc[k]/refvisc,c_stress_exp[k],c_delta_rho1[k]/c_delta_rho1[nodes-1],c_delta_rho2[k]/c_delta_rho2[nodes-1],c_therm_diff[k]);
if (densityfile==1) for(k=nodes-1;k>=0;k--) fprintf(ofile,"%f %f %e %f %.0f %e %f %f %f\n",1.0,1.0,c_alpha[k]/c_alpha[0],1.0,Tfinal[k],c_H[k],c_visc[k]/refvisc,c_stress_exp[k],c_therm_diff[k]);

printf("ref_alpha: %e\n",c_alpha[0]);

if (densityfile == 1){


perplex_data(c_alpha[0],Tadi,ref_density);
/*
    float** per_harz = (float **) malloc(sizeof(float*) * 143);
    float** per_morb = (float **) malloc(sizeof(float*) * 143);
    for(i=0;i<143;i++){
        per_harz[i] = malloc(sizeof(float) * 210);
        per_morb[i] = malloc(sizeof(float) * 210);}
    
    float* c_harz = malloc(sizeof(float) * nodes * 210);
    float* c_morb = malloc(sizeof(float) * nodes * 210);
   
    int depth_index;

sprintf(bufname1,"morb-harz-tackley.csv");
        ifile = fopen(bufname1,"r");

                for(i=0;i<143;i++){
                    for(k=0;k<210;k++){
                    getline (&bufline, &length, ifile);
                    sscanf (bufline, "%g\n", &per_morb[i][k]);
                    }}


   for (k=0;k<nodes;k++){
       depth_km = ((float) k)/((float)nodes-1.0)*mantle_depth_km;
       i_depth_km = (int) depth_km;
       depth_index = (i_depth_km-150)/20>0 ? (i_depth_km-150)/20 : 0;
//       printf("Depth index: %d\n",depth_index);
       for (i=0;i<210;i++){
       c_morb[k*210+i] = per_morb[depth_index][i]/(c_alpha[0]*ref_temperature);
}
}



}
*/
}
}

void perplex_data(float ref_alpha,float* Tadi, float ref_density)
{
    FILE *ifile1,*ifile2,*ofile;
    char *bufline, bufname1[200],bufname2[200];

    int minpress = 2;
    int deltapress = 5517;
    int npress = 291;
    float depth_km;
    int i_depth_km;
    int depth_index;
    float pressure[nodes];
    int i,k;
    size_t length = 80;
    bufline = (char*) malloc (80*sizeof(char));

    float*** per_harz = (float ***) malloc(sizeof(float**) * npress);
    float*** per_morb = (float ***) malloc(sizeof(float**) * npress);
    for(i=0;i<npress;i++){
        per_harz[i] = (float **) malloc(sizeof(float*) * ntemp);
        per_morb[i] = (float **) malloc(sizeof(float*) * ntemp);
        for(k=0;k<ntemp;k++){
            per_harz[i][k] = (float *) malloc(sizeof(float)*5);
            per_morb[i][k] = (float *) malloc(sizeof(float)*5);
        }
    }
    
    float*** c_harz = (float ***) malloc(sizeof(float**) * nodes);
    float*** c_morb = (float ***) malloc(sizeof(float**) * nodes);
    for(i=0;i<nodes;i++){
        c_morb[i] = (float **) malloc(sizeof(float*) * ntemp);
        c_harz[i] = (float **) malloc(sizeof(float*) * ntemp);
        for(k=0;k<ntemp;k++){
            c_morb[i][k] = (float *) malloc(sizeof(float)*6);
            c_harz[i][k] = (float *) malloc(sizeof(float)*5);
        }
    }

    sprintf(bufname1,"MORB-X_2.tab");
    ifile1 = fopen(bufname1,"r");
    sprintf(bufname1,"HARZ-X_2.tab");
    ifile2 = fopen(bufname1,"r");

                for(i=0;i<13;i++){
                    getline (&bufline, &length, ifile1);
                    getline (&bufline, &length, ifile2);
                    }


                for(i=0;i<npress;i++){
                    for(k=0;k<ntemp;k++){
                    getline (&bufline, &length, ifile1);
                    sscanf (bufline, "%f %f %f %f %f\n", &per_morb[i][k][0], &per_morb[i][k][1],&per_morb[i][k][2],&per_morb[i][k][3],&per_morb[i][k][4]);
                    getline (&bufline, &length, ifile2);
                    sscanf (bufline, "%f %f %f %f %f\n", &per_harz[i][k][0], &per_harz[i][k][1],&per_harz[i][k][2],&per_harz[i][k][3],&per_harz[i][k][4]);
                    //printf("%f\n",per_harz[i][k][0]);
                    }}

   for (k=0;k<nodes;k++){
       depth_km = ((float) k)/((float)nodes-1.0)*mantle_depth_km;
       //printf("%f",Tadi[k]);
       if (k==0) pressure[k] = 3400.0;
       else {pressure[k]= pressure[k-1] + 10.0*(mantle_depth_m/(nodes-1))*c_harz[k-1][(int)((Tadi[k-1]-mintemp)/deltatemp)][0]*ref_density/1e5;}
       //i_depth_km = (int) depth_km;
       if (k == nodes-1) printf("Pressure:%f[bar]",pressure[k]);
       depth_index = MAX((int)((pressure[k]-minpress)/deltapress),0);
       depth_index = MIN(depth_index,npress-1);
       //printf("Depth index: %d\n",depth_index);
       for (i=0;i<ntemp;i++){
       if ((i>=2) && ((per_harz[depth_index][i][0]/ref_density > 2.0) || (per_harz[depth_index][i][0]/ref_density < 0.5)))
         c_harz[k][i][0] = (2 * c_harz[k][i-1][0] - c_harz[k][i-2][0]);
       else if ((k>=2) && ((per_harz[depth_index][i][0]/ref_density > 2.0) || (per_harz[depth_index][i][0]/ref_density < 0.5)))
         c_harz[k][i][0] = (2 * c_harz[k-1][i][0] - c_harz[k-2][i][0]);
       else
         c_harz[k][i][0] = per_harz[depth_index][i][0]/ref_density;

       c_harz[k][i][0] = fmax(fmin(c_harz[k][i][0],2.0),0.5);

       if ((i>=2) && ((per_harz[depth_index][i][1]/ref_alpha > 2.0) || (per_harz[depth_index][i][1]/ref_alpha < 0.05)))
         c_harz[k][i][1] = (2 * c_harz[k][i-1][1] - c_harz[k][i-2][1]);
       else if ((k>=2) && ((per_harz[depth_index][i][1]/ref_alpha > 2.0) || (per_harz[depth_index][i][1]/ref_alpha < 0.05)))
         c_harz[k][i][1] = (2 * c_harz[k-1][i][1] - c_harz[k-2][i][1]);
       else
         c_harz[k][i][1] = per_harz[depth_index][i][1]/ref_alpha;

       c_harz[k][i][1] = fmax(fmin(c_harz[k][i][1],2.0),0.05);

       if ((i>=2) && ((per_harz[depth_index][i][2]/ref_cp > 2.0) || (per_harz[depth_index][i][2]/ref_cp < 0.2)))
         c_harz[k][i][2] = (2 * c_harz[k][i-1][2] - c_harz[k][i-2][2]);
       else if ((k>=2) && ((per_harz[depth_index][i][2]/ref_cp > 2.0) || (per_harz[depth_index][i][2]/ref_cp < 0.2)))
         c_harz[k][i][2] = (2 * c_harz[k-1][i][2] - c_harz[k-2][i][2]);
       else
         c_harz[k][i][2] = per_harz[depth_index][i][2]/ref_cp;

       c_harz[k][i][2] = fmax(fmin(c_harz[k][i][2],2.0),0.2);


       c_harz[k][i][3] = fmax(fmin(per_harz[depth_index][i][3],15.0),5.9);
       c_harz[k][i][4] = fmax(fmin(per_harz[depth_index][i][4],8.0),4.5);

       if ((i>=2) && ((per_morb[depth_index][i][0]/ref_density > 2.0) || (per_morb[depth_index][i][0]/ref_density < 0.5)))
         c_morb[k][i][0] = (2 * c_morb[k][i-1][0] - c_morb[k][i-2][0]);
       else if ((k>=2) && ((per_morb[depth_index][i][0]/ref_density > 2.0) || (per_morb[depth_index][i][0]/ref_density < 0.5)))
         c_morb[k][i][0] = (2 * c_morb[k-1][i][0] - c_morb[k-2][i][0]);
       else
         c_morb[k][i][0] = per_morb[depth_index][i][0]/ref_density;

       c_morb[k][i][0] = fmax(fmin(c_morb[k][i][0],2.0),0.5);

       if ((i>=2) && ((per_morb[depth_index][i][1]/ref_alpha > 2.0) || (per_morb[depth_index][i][1]/ref_alpha < 0.05)))
         c_morb[k][i][1] = (2 * c_morb[k][i-1][1] - c_morb[k][i-2][1]);
       else if ((k>=2) && ((per_morb[depth_index][i][1]/ref_alpha > 2.0) || (per_morb[depth_index][i][1]/ref_alpha < 0.05)))
         c_morb[k][i][1] = (2 * c_morb[k-1][i][1] - c_morb[k-2][i][1]);
       else
         c_morb[k][i][1] = per_morb[depth_index][i][1]/ref_alpha;

       c_morb[k][i][1] = fmax(fmin(c_morb[k][i][1],2.0),0.05);

       if ((i>=2) && ((per_morb[depth_index][i][2]/ref_cp > 2.0) || (per_morb[depth_index][i][2]/ref_cp < 0.2)))
         c_morb[k][i][2] = (2 * c_morb[k][i-1][2] - c_morb[k][i-2][2]);
       else if ((k>=2) && ((per_morb[depth_index][i][2]/ref_cp > 2.0) || (per_morb[depth_index][i][2]/ref_cp < 0.2)))
         c_morb[k][i][2] = (2 * c_morb[k-1][i][2] - c_morb[k-2][i-2][2]);
       else
         c_morb[k][i][2] = per_morb[depth_index][i][2]/ref_cp;

       c_morb[k][i][2] = fmax(fmin(c_morb[k][i][2],2.0),0.2);

       c_morb[k][i][3] = fmax(fmin(c_morb[k][i][0]-c_harz[k][i][0],0.5),-0.5);
       c_morb[k][i][3] = c_morb[k][i][3]/(ref_alpha*ref_temperature);

       c_morb[k][i][4] = fmax(fmin(per_morb[depth_index][i][3],15.0),5.9);
       c_morb[k][i][5] = fmax(fmin(per_morb[depth_index][i][4],8.0),4.5);
}
}

// Filter
    int j,n,m,numpoints;


    float*** c_harz_f = (float ***) malloc(sizeof(float**) * nodes);
    float*** c_morb_f = (float ***) malloc(sizeof(float**) * nodes);
    for(i=0;i<nodes;i++){
        c_morb_f[i] = (float **) malloc(sizeof(float*) * ntemp);
        c_harz_f[i] = (float **) malloc(sizeof(float*) * ntemp);
        for(k=0;k<ntemp;k++){
            c_morb_f[i][k] = (float *) malloc(sizeof(float)*6);
            c_harz_f[i][k] = (float *) malloc(sizeof(float)*5);
        }
    }
    
    for (i=0;i<nodes;i++){
        for (j=0;j<ntemp;j++){
            for(k=0;k<3;k++){
                c_morb_f[i][j][k]=0.0;
                c_harz_f[i][j][k]=0.0;
                numpoints = 0;
                for (n=-1;n<=1;n++){
                    for (m=-1;m<=1;m++){
                        if ((i+n>=0) && (i+n<nodes) && (j+m>=0) && (j+m<ntemp)){
                            c_morb_f[i][j][k] += c_morb[i+n][j+m][k];
                            c_harz_f[i][j][k] += c_harz[i+n][j+m][k];
                            numpoints++;
                        }
                    }
                }
                c_morb_f[i][j][k]/=numpoints;
                c_harz_f[i][j][k]/=numpoints;
            }
            c_harz_f[i][j][3]=c_harz[i][j][3];
            c_morb_f[i][j][3]=c_morb[i][j][3];
            c_harz_f[i][j][4]=c_harz[i][j][4];
            c_morb_f[i][j][4]=c_morb[i][j][4];
            c_morb_f[i][j][5]=c_morb[i][j][5];
        }
    }
    c_morb=c_morb_f;
    c_harz=c_harz_f;

//                if ((i>0) && (i<64) && (j>0) && (j<64)){
//                                    V[age][i][j][0] = temp[0]/numpoints;
//                                                        V[age][i][j][1] = temp[1]/numpoints;
//                                                        //                } else {
//                                                        //                    V[age][i][j][0] = 0.0;
//                                                        //                    V[age][i][j][1] = 0.0;
//                                                        //                }
//                                                                    }
//                                                                            }
//


sprintf(bufname2,"perplex.dat");
        ofile = fopen(bufname2,"w");

       for(k=nodes-1;k>=0;k--) 
         for(i=0;i<ntemp;i++){
         if (c_harz[k][i][0] < 1e-2) printf("ERROR Rho too small");
         /*fprintf(ofile,"1.0 1.0 1.0\n");
         fprintf(ofile,"0.0 0.0 0.0\n");
         fprintf(ofile,"0.0 0.0 0.0\n");*/
         fprintf(ofile,"%f %f %f %f %f\n",c_harz[k][i][0],c_harz[k][i][1],c_harz[k][i][2],c_harz[k][i][3],c_harz[k][i][4]);
         fprintf(ofile,"%f %f %f %f %f\n",c_morb[k][i][0]-c_harz[k][i][0],c_morb[k][i][1]-c_harz[k][i][1],c_morb[k][i][2]-c_harz[k][i][2],c_morb[k][i][4]-c_harz[k][i][3],c_morb[k][i][5]-c_harz[k][i][4]);
         fprintf(ofile,"%f %f %f %f %f\n",c_harz[k][i][0]-c_harz[k][i][0],c_harz[k][i][1]-c_harz[k][i][1],c_harz[k][i][2]-c_harz[k][i][2],c_harz[k][i][3]-c_harz[k][i][3]);
}
fclose(ofile);

sprintf(bufname2,"densityfile.dat");
        ofile = fopen(bufname2,"w");

       for(k=nodes-1;k>=0;k--) 
         for(i=0;i<ntemp;i++){
         fprintf(ofile,"%f\n",c_morb[k][i][3]);
         fprintf(ofile,"-0.1\n");
}
fclose(ofile);

output_vtk(c_morb,c_harz);
}

/*------------------------------------------------------------------------------header----------------------------------------------------------------------*/

static void vts_file_header(FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n"
        "  <StructuredGrid WholeExtent=\"%s\">\n"
        "    <Piece Extent=\"%s\">\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
             0, nodes-1,
             0, ntemp-1,
             0, 0);

    snprintf(header, 1024, format, extent, extent);

    fputs(header, fp);

    return;
}

/*------------------------------------------------------------------------------trailer----------------------------------------------------------------------*/

static void vts_file_trailer(FILE *fp)
{
    const char trailer[] =
        "    </Piece>\n"
        "  </StructuredGrid>\n"
        "</VTKFile>\n";

    fputs(trailer, fp);

    return;
}

/*------------------------------------------------------------------------point and cell data-----------------------------------------------------------------*/

static void vtk_point_data_header(FILE *fp)
{
    fputs("      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n", fp);
    return;
}


static void vtk_point_data_trailer(FILE *fp)
{
    fputs("      </PointData>\n", fp);
    return;
}


static void vtk_cell_data_header(FILE *fp)
{
    fputs("      <CellData>\n", fp);
    return;
}


static void vtk_cell_data_trailer(FILE *fp)
{
    fputs("      </CellData>\n", fp);
    return;
}

/*------------------------------------------------------------------------------cp----------------------------------------------------------------------*/

static void vtk_output_cp(float*** c_morb,float*** c_harz,FILE *fp)
{
    int i, k;

    fputs("        <DataArray type=\"Float32\" Name=\"heat capacity\" NumberOfComponents=\"2\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            fprintf(fp, "%.4e\n", c_harz[k][i][2]);
            fprintf(fp, "%.4e\n", c_morb[k][i][2]);
        }

    fputs("        </DataArray>\n", fp);
    return;
}

/*------------------------------------------------------------------------------seismic----------------------------------------------------------------------*/

static void vtk_output_seismic(float*** c_morb,float*** c_harz,FILE *fp)
{
    int i, k;

    fputs("        <DataArray type=\"Float32\" Name=\"vp\" NumberOfComponents=\"2\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            fprintf(fp, "%.4e\n", c_harz[k][i][3]);
            fprintf(fp, "%.4e\n", c_morb[k][i][4]);
        }

    fputs("        </DataArray>\n", fp);

    fputs("        <DataArray type=\"Float32\" Name=\"vs\" NumberOfComponents=\"2\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            fprintf(fp, "%.4e\n", c_harz[k][i][4]);
            fprintf(fp, "%.4e\n", c_morb[k][i][5]);
        }

    fputs("        </DataArray>\n", fp);
    return;
}

/*------------------------------------------------------------------------------alpha----------------------------------------------------------------------*/

static void vtk_output_alpha(float*** c_morb,float*** c_harz,FILE *fp)
{
    int i, k;

    fputs("        <DataArray type=\"Float32\" Name=\"thermal expansivity\" NumberOfComponents=\"2\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            fprintf(fp, "%.4e\n", c_harz[k][i][1]);
            fprintf(fp, "%.4e\n", c_morb[k][i][1]);
        }

    fputs("        </DataArray>\n", fp);
    return;
}

/*------------------------------------------------------------------------------density----------------------------------------------------------------------*/

static void vtk_output_density(float*** c_morb,float*** c_harz,FILE *fp)
{
    int i, k;

    fputs("        <DataArray type=\"Float32\" Name=\"density\" NumberOfComponents=\"2\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            fprintf(fp, "%.4e\n", c_harz[k][i][0]);
            fprintf(fp, "%.4e\n", c_morb[k][i][0]);
        }

    fputs("        </DataArray>\n", fp);
    return;
}

/*------------------------------------------------------------------------------density----------------------------------------------------------------------*/

static void vtk_output_density_contrast(float*** c_morb,FILE *fp)
{
    int i, k;

    fputs("        <DataArray type=\"Float32\" Name=\"density contrast\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++)
            fprintf(fp, "%.4e\n", c_morb[k][i][3]);

    fputs("        </DataArray>\n", fp);
    return;
}

/*------------------------------------------------------------------------------coordinates----------------------------------------------------------------------*/

static void vtk_output_coord(FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, k;

    fputs("      <Points>\n", fp);
    fputs("        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

        for(i=0; i<ntemp; i++)
          for(k=0; k<nodes; k++){
            depth_km = ((float) k)/((float)nodes-1.0)*mantle_depth_km;
            fprintf(fp,"%.6e %.6e %.6e\n",
                    mantle_depth_km - depth_km, 
                    (float) (i * deltatemp + mintemp),
                    0.0);}

    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);

    return;
}

/*------------------------------------------------------------------------------output_open----------------------------------------------------------------------*/

FILE* output_open(char *filename, char *mode)
{
  FILE *fp1;

  /* if filename is empty, output to stderr. */
  if (*filename) {
    fp1 = fopen(filename,mode);
    if (!fp1) {
      fprintf(stderr,"Cannot open file '%s' for '%s'\n",
	      filename,mode);
    }
  }
  else
    fp1 = stderr;

  return fp1;
}

/*------------------------------------------------------------------------------vtk-output----------------------------------------------------------------------*/

void output_vtk(float*** c_morb, float*** c_harz)
{
    char output_file[255];
    FILE *fp;

    snprintf(output_file, 255, "refstate.%d.vts",0);
    fp = output_open(output_file, "w");


    /* first, write volume data to vts file */
    vts_file_header(fp);

    /* write node-based field */
    vtk_point_data_header(fp);

    vtk_output_density(c_morb,c_harz,fp);
    vtk_output_alpha(c_morb,c_harz,fp);
    vtk_output_cp(c_morb,c_harz,fp);
    vtk_output_density_contrast(c_morb,fp);
    vtk_output_seismic(c_morb,c_harz,fp);
    
    vtk_point_data_trailer(fp);

    /* write element-based field */
    vtk_cell_data_header(fp);
    vtk_cell_data_trailer(fp);

    /* write coordinate */
    vtk_output_coord(fp);

    vts_file_trailer(fp);

    /* then, write other type of data */   

    fclose(fp);


    return;
}

