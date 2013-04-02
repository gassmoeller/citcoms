/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Routine to process the output of the finite element cycles
   and to turn them into a coherent suite  files  */


#include <math.h>
#include <unistd.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "output.h"
#include "string.h"
#include "material_properties.h"

#include "zlib.h"

#define CHUNK 16384

static void write_binary_array_double(int nn, double* array, FILE * f);
static void write_binary_array(int nn, float* array, FILE * f);
static void write_int_binary_array(int nn, int* array, FILE * f);
static void write_ascii_array_float(int nn, int perLine, float *array, FILE *fp);
static void write_ascii_array_double(int nn, int perLine, double *array, FILE *fp);
static void write_ascii_array(int nn, int perLine, int ascii_precision, double *array, FILE *fp);

// TODO: Write a unit test for this function. Should be pretty easy except from including this source file.
void get_compressor_string(int is_not_binary, int string_length, char* compressor_string)
{
	if (is_not_binary == 0)
		snprintf (compressor_string,string_length," compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\"");
	else
		snprintf (compressor_string,string_length,"");
}

static void vts_file_header(struct All_variables *E, FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\"%s>\n"
        "  <StructuredGrid WholeExtent=\"%s\">\n"
        "    <Piece Extent=\"%s\">\n";

    char extent[64], compressor_string[128], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
             E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz,
             E->lmesh.exs, E->lmesh.exs + E->lmesh.elx,
             E->lmesh.eys, E->lmesh.eys + E->lmesh.ely);

    get_compressor_string (strcmp(E->output.vtk_format,"binary"),128,compressor_string);

    snprintf(header, 1024, format, compressor_string, extent, extent);

    fputs(header, fp);

    return;
}


static void vts_file_trailer(struct All_variables *E, FILE *fp)
{
    const char trailer[] =
        "    </Piece>\n"
        "  </StructuredGrid>\n"
        "</VTKFile>\n";

    fputs(trailer, fp);

    return;
}


static void vtk_point_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <PointData Scalars=\"temperature\" Vectors=\"velocity\">\n", fp);


    return;
}


static void vtk_point_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </PointData>\n", fp);
    return;
}


static void vtk_cell_data_header(struct All_variables *E, FILE *fp)
{
    fputs("      <CellData>\n", fp);
    return;
}


static void vtk_cell_data_trailer(struct All_variables *E, FILE *fp)
{
    fputs("      </CellData>\n", fp);
    return;
}


static void vtk_output_temp(struct All_variables *E, FILE *fp)
{
    int i;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\">\n", E->output.vtk_format);

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array_double(nodes,E->T[1]+1,fp);
    } else {
        write_ascii_array_double(nodes,1,E->T[1]+1,fp);
    }
    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_melttemp(struct All_variables *E, FILE *fp)
{
    int i,nz;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    double* melttemp = malloc(nodes*sizeof(double));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"melt excess temperature\" format=\"%s\">\n", E->output.vtk_format);

    for(i=0;i < nodes;i++){
        nz = i % E->lmesh.noz;
        melttemp[i] =  (*(E->T[1]+i+1) + E->control.surface_temp) * E->data.ref_temperature - E->refstate.Tm[nz*E->composition.pressure_oversampling + 1];
    }

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array_double(nodes,melttemp,fp);
    } else {
        write_ascii_array_double(nodes,1,melttemp,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(melttemp);
    return;
}

static void vtk_output_deltat(struct All_variables *E, FILE *fp)
{
    int i,nz;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    double* deltatemp = malloc(nodes*sizeof(double));

    compute_horiz_avg(E);

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"DeltaT [K]\" format=\"%s\">\n", E->output.vtk_format);

    for(i=0;i <= nodes;i++){ 
        nz = i % E->lmesh.noz + 1;
        deltatemp[i] = (*(E->T[1]+i+1)-E->Have.T[nz]) * E->data.ref_temperature;
    }
        

    if (strcmp(E->output.vtk_format,"binary") == 0) {
        write_binary_array_double(nodes,deltatemp,fp);
    } else {
        write_ascii_array_double(nodes,1,deltatemp,fp);
    }
    fputs("        </DataArray>\n", fp);
    free(deltatemp);
    return;
}

static void vtk_output_velo(struct All_variables *E, FILE *fp)
{
    int i, j;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double sint, sinf, cost, cosf;
    float *V[4];
    const int lev = E->mesh.levmax;
    double* vel = malloc(nodes*3*sizeof(double));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        V[1] = E->sphere.cap[j].V[1];
        V[2] = E->sphere.cap[j].V[2];
        V[3] = E->sphere.cap[j].V[3];

        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            vel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+0] = V[1][i]*cost*cosf - V[2][i]*sinf + V[3][i]*sint*cosf;
            vel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+1] = V[1][i]*cost*sinf + V[2][i]*cosf + V[3][i]*sint*sinf;
            vel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+2] = -V[1][i]*sint + V[3][i]*cost;
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes*3,vel,fp);
    else 
        write_ascii_array_double(nodes*3,3,vel,fp);
    fputs("        </DataArray>\n", fp);

    free(vel);
    return;
}

static void vtk_output_material(struct All_variables *E, FILE *fp)
{
    int i, m;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double* alpha = malloc(nodes*sizeof(double));
    double* rho = malloc(nodes*sizeof(double));
    double* cp = malloc(nodes*sizeof(double));
    double* radheat = malloc(nodes*sizeof(double));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"rho\" format=\"%s\">\n", E->output.vtk_format);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++) {
            rho[((m-1)*E->sphere.caps_per_proc)+i-1] = get_rho_nd(E,m,i);
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes,rho,fp);
    else 
        write_ascii_array_double(nodes,1,rho,fp);
    fputs("        </DataArray>\n", fp);


    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"alpha\" format=\"%s\">\n", E->output.vtk_format);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++) {
            alpha[((m-1)*E->sphere.caps_per_proc)+i-1] = get_alpha_nd(E,m,i);
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes,alpha,fp);
    else 
        write_ascii_array_double(nodes,1,alpha,fp);
    fputs("        </DataArray>\n", fp);


    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"cp\" format=\"%s\">\n", E->output.vtk_format);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++) {
            cp[((m-1)*E->sphere.caps_per_proc)+i-1] = get_cp_nd(E,m,i);
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes,cp,fp);
    else 
        write_ascii_array_double(nodes,1,cp,fp);
    fputs("        </DataArray>\n", fp);

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"radiogenic heating\" format=\"%s\">\n", E->output.vtk_format);

    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=1; i<=E->lmesh.nno; i++) {
            radheat[((m-1)*E->sphere.caps_per_proc)+i-1] = get_radheat_nd(E,m,i);
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes,radheat,fp);
    else 
        write_ascii_array_double(nodes,1,radheat,fp);
    fputs("        </DataArray>\n", fp);

    free(alpha);
    free(rho);
    free(cp);
    free(radheat);
    return;
}

static void vtk_output_svelo(struct All_variables *E, FILE *fp)
{
    int i, j;
    int nodes=E->sphere.caps_per_proc*E->lmesh.nno;
    double* svel = malloc(nodes*3*sizeof(double));

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"spherical velocity\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            svel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+0] = E->sphere.cap[j].V[1][i];
            svel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+1] = E->sphere.cap[j].V[2][i];
            svel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3+2] = E->sphere.cap[j].V[3][i];
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0) 
        write_binary_array_double(nodes*3,svel,fp);
    else 
        write_ascii_array_double(nodes*3,3,svel,fp);
    fputs("        </DataArray>\n", fp);

    free(svel);
    return;
}

void vtk_output_seismic(struct All_variables *E, int cycles, FILE *fp)
{
    void get_prem(double, double*, double*, double*);
    void compute_seismic_model(const struct All_variables*, double*, double*, double*);

    char output_file[255];
    int i,nz;
    FILE *fp2;
    char prem_file[255];

    double **rho, **vp, **vs;
    double *output_rho, *output_vp, *output_vs;
    double *arho, *avp, *avs;
    double premrho,premvp,premvs;
    const int nodes = E->lmesh.nno + 1;

    rho = (double **) malloc(2 * sizeof(double*));
    vp = (double **) malloc(2 * sizeof(double*));
    vs = (double **) malloc(2 * sizeof(double*));
    rho[1] = (double *) malloc(nodes * sizeof(double));
    vp[1] = (double *) malloc(nodes * sizeof(double));
    vs[1] = (double *) malloc(nodes * sizeof(double));
    output_rho = malloc(nodes * sizeof(double));
    output_vp = malloc(nodes * sizeof(double));
    output_vs = malloc(nodes * sizeof(double));
    arho = (double *)malloc( (E->lmesh.noz+1)*sizeof(double));
    avp = (double *)malloc( (E->lmesh.noz+1)*sizeof(double));
    avs = (double *)malloc( (E->lmesh.noz+1)*sizeof(double));


    if(rho==NULL || vp==NULL || vs==NULL) {
        fprintf(stderr, "Error while allocating memory\n");
        abort();
    }

    if (E->refstate.choice != 3){
        /* isotropic seismic velocity only */
        /* XXX: update for anisotropy in the future */
    /*    compute_seismic_model(E, rho, vp, vs);
        for (i=0;i<nodes;i++){
            floatrho[i] = (float) rho[i];
            floatvp[i]  = (float) vp[i];
            floatvs[i]  = (float) vs[i];
        }*/
    } else{
        for (i=1;i<=nodes;i++){
            rho[1][i]= get_rho_nd(E,1,i);
            vp[1][i] = get_vp_nd(E,1,i);
            vs[1][i] = get_vs_nd(E,1,i);
        }
        remove_horiz_ave(E,rho,arho,0);
        remove_horiz_ave(E,vp,avp,0);
        remove_horiz_ave(E,vs,avs,0);

        for (i=0;i<nodes;i++){
            nz = i % E->lmesh.noz;
            output_rho[i] = rho[1][i+1]*100.0/arho[nz+1];
            output_vp[i] = vp[1][i+1]*100.0/avp[nz+1];
            output_vs[i] = vs[1][i+1]*100.0/avs[nz+1];
        }
        if (E->parallel.me < E->parallel.nprocz){
            snprintf(prem_file, 255, "%s.prem.%d.%d.csv", E->control.data_file, E->parallel.me, cycles);
            fp2 = output_open(prem_file, "w");
            fprintf(fp2,"Radius PremVP PremVS PremRho AverageVP AverageVS AverageRho\n");
            for (nz=1;nz <= E->lmesh.noz;nz++){
                get_prem(E->sx[1][3][nz], &premvp, &premvs, &premrho);
                fprintf(fp2,"%f %f %f %f %f %f %f\n", E->sx[1][3][nz], premvp,premvs,premrho*1000,avp[nz],avs[nz],arho[nz]*E->data.density);
            }
            fclose(fp2);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0){
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Rho\" format=\"%s\">\n", E->output.vtk_format);
        write_binary_array_double(nodes-1,output_rho,fp);
        fputs("        </DataArray>\n", fp);
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Vp\" format=\"%s\">\n", E->output.vtk_format);
        write_binary_array_double(nodes-1,output_vp,fp);
        fputs("        </DataArray>\n", fp);
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Vs\" format=\"%s\">\n", E->output.vtk_format);
        write_binary_array_double(nodes-1,output_vs,fp);
        fputs("        </DataArray>\n", fp);
    }else {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Rho\" format=\"%s\">\n", E->output.vtk_format);
        write_ascii_array_double(nodes-1,1,output_rho,fp);
        fputs("        </DataArray>\n", fp);
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Vp\" format=\"%s\">\n", E->output.vtk_format);
        write_ascii_array_double(nodes-1,1,output_vp,fp);
        fputs("        </DataArray>\n", fp);
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Seismic Vs\" format=\"%s\">\n", E->output.vtk_format);
        write_ascii_array_double(nodes-1,1,output_vs,fp);
        fputs("        </DataArray>\n", fp);
    }


#if 0
    /** debug **/
    sprintf(output_file,"%s.dv.%d.%d", E->control.data_file, E->parallel.me, cycles);
    fp = output_open(output_file, "w");
    fprintf(fp, "%d %d %.5e\n", cycles, E->lmesh.nno, E->monitor.elapsed_time);
    for(i=0; i<E->lmesh.nno; i++) {
        double vpr, vsr, rhor;
        int nz = (i % E->lmesh.noz) + 1;
        get_prem(E->sx[1][3][nz], &vpr, &vsr, &rhor);

        fprintf(fp, "%.4e %.4e %.4e\n",
                rho[i]/rhor - 1.0,
                vp[i]/vpr - 1.0,
                vs[i]/vsr - 1.0);

    }
    fclose(fp);
#endif

    for (i=1;i<2;i++){
    free(rho[i]);
    free(vp[i]);
    free(vs[i]);
    }

    free(rho);
    free(vp);
    free(vs);
    free(output_rho);
    free(output_vp);
    free(output_vs);
    free(arho);
    free(avp);
    free(avs);
    return;
}


static void vtk_output_visc(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    int lev = E->mesh.levmax;

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array(nodes,&E->VI[lev][1][1],fp);
        } else {
            write_ascii_array_float(nodes,1,&E->VI[lev][1][1],fp);
        }
       
    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_dens(struct All_variables *E, FILE *fp)
{
    int lev = E->mesh.levmax;
    int i,j,k,nz,nT;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    double* density = malloc(nodes*sizeof(double));

        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nno; i++) {
                nz = ((i-1) % E->lmesh.noz) + 1;
                density[(j-1)*E->lmesh.nno+i-1] = (E->buoyancy[j][i]/(E->control.Atemp*E->refstate.gravity[nz]*get_rho_nd(E,j,i)));
	    }
        }

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"density\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array_double(nodes,density,fp);
        } else {
            write_ascii_array_double(nodes,1,density,fp);
        }
       
    fputs("        </DataArray>\n", fp);
    free(density);
    return;
}

static void vtk_output_tracer(struct All_variables *E, FILE *fp)
{
    int lev = E->mesh.levmax;
    int i,j,numtracers,flavor;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floattracer = malloc(nodes*sizeof(float));

        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nel; i++) {
                numtracers = 0;
                for (flavor=0; flavor<E->trace.nflavors; flavor++)
                    numtracers += E->trace.ntracer_flavor[j][flavor][i];
                floattracer[(j-1)*E->lmesh.nel+i-1] = (float)(numtracers);
	    }
        }

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"numtracer\" format=\"%s\">\n", E->output.vtk_format);
        if (strcmp(E->output.vtk_format, "binary") == 0) {
            write_binary_array(nodes,floattracer,fp);
        } else {
            write_ascii_array_float(nodes,1,floattracer,fp);
        }
       
    fputs("        </DataArray>\n", fp);
    free(floattracer);
    return;
}

static void vtk_output_coord(struct All_variables *E, FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, j;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    double* pos = malloc(nodes*3*sizeof(double));

    fputs("      <Points>\n", fp);
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++){
                pos[((j-1)*E->lmesh.nno+i-1)*3] = (E->x[j][1][i]);
	        pos[((j-1)*E->lmesh.nno+i-1)*3+1] = (E->x[j][2][i]);
	        pos[((j-1)*E->lmesh.nno+i-1)*3+2] = (E->x[j][3][i]);
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array_double(nodes*3,pos,fp);
    else
        write_ascii_array_double(nodes*3,3,pos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(pos);
    return;
}

static void vtk_output_stress(struct All_variables *E, FILE *fp)
{
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
   /* for stress computation */
    void allocate_STD_mem();
    void compute_nodal_stress();
    void free_STD_mem();
    float *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
    float *divv[NCS],*vorv[NCS];

    /* those are sorted like stt spp srr stp str srp  */
    allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\">\n", E->output.vtk_format);

    if (strcmp(E->output.vtk_format, "binary") == 0) {
        write_binary_array(nodes*6,&E->gstress[1][1],fp);
    } else {
        write_ascii_array_float(nodes*6,6,&E->gstress[1][1],fp);
    }

    fputs("        </DataArray>\n", fp);
    return;
}

static void vtk_output_comp_nd(struct All_variables *E, FILE *fp)
{
    int i, j, k;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    double* compo = malloc (nodes*sizeof(double));
    
    for(k=0;k<E->composition.ncomp;k++) {
        fprintf(fp, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\">\n", k+1, E->output.vtk_format);

        for(j=1; j<=E->sphere.caps_per_proc; j++) {
            for(i=1; i<=E->lmesh.nno; i++) {
                compo[(j-1)*E->lmesh.nno+i-1] = (E->composition.comp_node[j][k][i]);
	    }
        }

        if (strcmp(E->output.vtk_format, "binary") == 0)
            write_binary_array_double(nodes,compo,fp);
        else
            write_ascii_array_double(nodes,1,compo,fp);
        fputs("        </DataArray>\n", fp);
    }
    free(compo);
    return;
}


void vtk_output_surf_botm(struct All_variables *E,  FILE *fp, int cycles)
{
    int i, j, k;
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    char output_file[255];
    float* floattopo = malloc (nodes*sizeof(float));
    float* floatheating = malloc (nodes*sizeof(float));

    if((E->output.write_q_files == 0) || (cycles == 0) ||
      (cycles % E->output.write_q_files)!=0)
        heat_flux(E);
  /* else, the heat flux will have been computed already */

    if(E->control.use_cbf_topo){
        get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);
    }
    else{
        get_STD_topo(E,E->slice.tpg,E->slice.tpgb,E->slice.divg,E->slice.vort,cycles);
    }

    fprintf(fp,"        <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1;j<=E->sphere.caps_per_proc;j++){
        for(i=1;i<=E->lmesh.nsf;i++){
            for(k=1;k<=E->lmesh.noz;k++){
                floattopo[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz + k-1] = 0.0; 
            }

            if (E->output.surf && (E->parallel.me_loc[3]==E->parallel.nprocz-1)) {

                /* choose either STD topo or pseudo-free-surf topo */
                if(E->control.pseudo_free_surf)
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.freesurf[j][i];
                else
                floattopo[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.tpg[j][i];

            }

            if (E->output.botm && (E->parallel.me_loc[3]==0)){

                /* choose either STD topo or pseudo-free-surf topo */
                if(E->control.pseudo_free_surf)
                /* is this integrated at the moment? */
                /*floattopo[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz] = E->slice.freesurf[j][i];*/
                fprintf(stderr,"Bottom topography for pseudo free surface currently not possible");
                else
                floattopo[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz] = E->slice.tpgb[j][i];

            }
          
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes,floattopo,fp);
    else
        write_ascii_array_float(nodes,1,floattopo,fp);

    fputs("        </DataArray>\n", fp);

    fprintf(fp,"        <DataArray type=\"Float32\" Name=\"heatflux\" format=\"%s\">\n", E->output.vtk_format);
    for(j=1;j<=E->sphere.caps_per_proc;j++){
        for(i=1;i<=E->lmesh.nsf;i++){
            for(k=1;k<=E->lmesh.noz;k++)
                floatheating[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz + k-1] = 0.0; 
            if (E->parallel.me_loc[3]==E->parallel.nprocz-1)
              floatheating[(j-1)*E->lmesh.nno + i*E->lmesh.noz-1] = E->slice.shflux[j][i];
            if (E->parallel.me_loc[3]==0)
              floatheating[(j-1)*E->lmesh.nno + (i-1)*E->lmesh.noz] = E->slice.bhflux[j][i];
        }
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(nodes,floatheating,fp);
    else
        write_ascii_array_float(nodes,1,floatheating,fp);

    fputs("        </DataArray>\n", fp);

  free(floattopo);
  free(floatheating);

  return;
}

void write_pvtp(struct All_variables *E, int cycles)
{
    FILE *fp;
    char pvts_file[255];
    int i,j,k,l;
    snprintf(pvts_file, 255, "%s.tracer_file.%d.pvtp",
    E->control.data_file,cycles);
    fp = output_open(pvts_file, "w");

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"PPolyData\" version=\"0.1\"%s>\n"
        "  <PPolyData GhostLevel=\"#\">\n"
        "    <PPointData Scalars=\"Composition\" Vectors=\"velocity\">\n";

    char extent[64], compressor_string[128], header[1024];

    get_compressor_string (strcmp(E->output.vtk_format,"binary"),128,compressor_string);

    snprintf(header, 1024, format, compressor_string);
    fputs(header, fp);

    for (i = 0; i<E->trace.number_of_extra_quantities;i++)
        fprintf(fp, "      <DataArray type=\"Float32\" Name=\"Tracer Quantity %d\" format=\"%s\"/>\n", i, E->output.vtk_format);

    fputs("    </PPointData>\n \n"
    "    <PCellData>\n",fp);

    fputs("    </PCellData>\n \n"
    "    <PPoints>\n"
    "      <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\" />\n"
    "    </PPoints>\n", fp);

    for (i=0;i<E->parallel.nproc;i++){
        fprintf(fp, "    <Piece Source=\"%s.tracer_file.%d.%d.vtp\"/>\n",
                E->control.data_prefix, 
                i, cycles);
    }

    fputs("  </PPolyData>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}

void write_pvts(struct All_variables *E, int cycles)
{
    FILE *fp;
    char pvts_file[255];
    int i,j,k,l;
    l = E->parallel.me/(E->parallel.nprocx*E->parallel.nprocy*E->parallel.nprocz);
    snprintf(pvts_file, 255, "%s.%d.%d.pvts",
    E->control.data_file,l,cycles);
    fp = output_open(pvts_file, "w");

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        "  <PStructuredGrid WholeExtent=\"%s\" GhostLevel=\"#\">\n"
        "    <PPointData Scalars=\"temperature\" Vectors=\"velocity\">\n"
        "      <DataArray type=\"Float32\" Name=\"temperature\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"%s\"/>\n"
        "      <DataArray type=\"Float32\" Name=\"viscosity\" format=\"%s\"/>\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
        E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz*E->parallel.nprocz,
        E->lmesh.exs, E->lmesh.exs + E->lmesh.elx*E->parallel.nprocx,
        E->lmesh.eys, E->lmesh.eys + E->lmesh.ely*E->parallel.nprocy);

    snprintf(header, 1024, format, extent, E->output.vtk_format, 
             E->output.vtk_format, E->output.vtk_format);
    fputs(header, fp);

    if (E->output.stress){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"%s\"/>\n", E->output.vtk_format);
    }
    if (E->output.comp_nd && E->composition.on){
        for (i=0;i<E->composition.ncomp;i++)
            fprintf(fp,"      <DataArray type=\"Float32\" Name=\"composition%d\" format=\"%s\"/>\n",i+1, E->output.vtk_format);
    }
    if (E->output.surf || E->output.botm){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"surface\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"heatflux\" format=\"%s\"/>\n", E->output.vtk_format); 
    }
    if (E->output.density){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"density\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    if (E->output.svelo){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"spherical velocity\" NumberOfComponents=\"3\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    if (E->output.seismic){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"seismic rho\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"seismic vp\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"seismic vs\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    if (E->output.material){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"rho\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"alpha\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"cp\" format=\"%s\"/>\n", E->output.vtk_format); 
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"radiogenic heating\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    if (E->output.deltat){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"deltaT [K]\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    if (E->output.melttemp){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"melt excess temperature [K]\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    fputs("    </PPointData>\n \n"
    "    <PCellData>\n",fp);

    if (E->output.tracer && E->composition.on){
        fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Numtracer\" NumberOfComponents=\"1\" format=\"%s\"/>\n", E->output.vtk_format); 
    }

    fputs("    </PCellData>\n \n"
    "    <PPoints>\n"
    "      <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\" />\n"
    "    </PPoints>\n", fp);

    for(i=0; i < E->parallel.nprocy;i++){
        for(j=0; j < E->parallel.nprocx;j++){
            for(k=0; k < E->parallel.nprocz;k++){
                fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.proc%d.%d.vts\"/>\n",
                    (k%E->parallel.nprocz)*E->lmesh.elz, 
                    (k%E->parallel.nprocz+1)*E->lmesh.elz, 
                    (j%E->parallel.nprocx)*E->lmesh.elx, (j%E->parallel.nprocx+1)*E->lmesh.elx, 
                    (i%E->parallel.nprocy)*E->lmesh.ely, (i%E->parallel.nprocy+1)*E->lmesh.ely,
                    E->control.data_prefix, 
                    l*E->parallel.nproc/E->sphere.caps+i*E->parallel.nprocx*E->parallel.nprocz+j*E->parallel.nprocz+k, cycles);
            }
        }
    }

    fputs("  </PStructuredGrid>\n",fp);
    fputs("</VTKFile>",fp);

    fclose(fp);
}

void write_pvd(struct All_variables *E, int cycles)
{
    FILE *fp;
    char pvd_file[255];
    int i;
    snprintf(pvd_file, 255, "%s.pvd",
    E->control.data_file);

    if (cycles == 0){
        fp = output_open(pvd_file, "w");
        const char format[] =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"Collection\" version=\"0.1\"%s>\n"
            "  <Collection>\n";

        char compressor_string[128];

        get_compressor_string(strcmp(E->output.vtk_format,"binary"),128,compressor_string);

        fprintf(fp,format,compressor_string);
    }
    else
    {
        fp = output_open(pvd_file, "r+");
        fseek(fp, -27, SEEK_END);
    }

    for (i=0;i<E->sphere.caps;i++){
        fprintf(fp, "    <DataSet timestep=\"%.0f\" group=\"\" part=\"%d\" file=\"%s.%d.%d.pvts\"/>\n",E->monitor.elapsed_time*E->data.scalet*1000,i,E->control.data_prefix,i,cycles);
    }
    fflush(fp);
    const char format[] =
            "  </Collection>\n"
            "</VTKFile>\n";
    fputs(format,fp);
    fclose(fp);
}

static void vtk_tracer_extraq(struct All_variables *E, int idx_extraq, FILE *fp)
{
    int i, j;
    int tracers = 0;
    float* floatextraq = malloc (E->trace.ntracers[1]*sizeof(float)); // caps_per_proc != 1

    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Tracer Quantity %d\" format=\"%s\">\n", idx_extraq, E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->trace.ntracers[j]; i++) {
            floatextraq[tracers+i-1] = (float) (E->trace.extraq[j][idx_extraq][i]);
        }
        tracers += E->trace.ntracers[j];
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array(tracers,floatextraq,fp);
    else
        write_ascii_array_float(tracers,1,floatextraq,fp);
    fputs("        </DataArray>\n", fp);
    free(floatextraq);
    return;
}


static void vtk_tracer_coord(struct All_variables *E, FILE *ft)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. Currently just working for caps_per_proc == 1*/
    int i, j;
    int tracers = 0;
    double* pos = malloc (E->trace.ntracers[1]*3*sizeof(double)); // caps_per_proc != 1

    fputs("      <Points>\n", ft);
    fprintf(ft, "        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"%s\">\n", E->output.vtk_format);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->trace.ntracers[j]; i++){
                pos[(i-1)*3] = E->trace.basicq[j][3][i];
	        pos[(i-1)*3+1]= E->trace.basicq[j][4][i];
	        pos[(i-1)*3+2]= E->trace.basicq[j][5][i];
        }
    tracers += E->trace.ntracers[j];
    }

    if (strcmp(E->output.vtk_format, "binary") == 0)
        write_binary_array_double(tracers*3,pos,ft);
    else
        write_ascii_array_double(tracers*3,3,pos,ft);
    fputs("        </DataArray>\n", ft);
    fputs("      </Points>\n", ft);
    free(pos);
    return;
}

void write_tracer_file(struct All_variables *E, int cycles)
{
    FILE *ft;
    char vtp_file[255];
    char header[1024];
    int i;
    snprintf(vtp_file, 255, "%s.tracer_file.%d.%d.vtp",
    E->control.data_file,E->parallel.me,cycles);

    ft = output_open(vtp_file, "w");
    fprintf(ft,"<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"PolyData\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
            "  <PolyData>\n"
            "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
            "      <PointData Scalars=\"Composition\">\n", E->trace.ntracers[1],E->trace.ntracers[1]);

    for (i = 0; i<E->trace.number_of_extra_quantities;i++)
    	vtk_tracer_extraq(E,i,ft);

    vtk_point_data_trailer(E,ft);

    /* write element-based field */
    vtk_cell_data_header(E,ft);
    vtk_cell_data_trailer(E,ft);
    vtk_tracer_coord(E,ft);

    fputs("    <Verts>\n",ft);
    fputs("      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n",ft);

    int *index = malloc(sizeof(int)*(E->trace.ntracers[1]));
    for (i=0;i<E->trace.ntracers[1];i++){
        //fprintf(ft,"%d ",i);
        index[i] = i;}

    write_int_binary_array(E->trace.ntracers[1],index,ft);
    //fputs("\n",ft);
    fputs("      </DataArray>\n",ft);
    fputs("      <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n",ft);
    for (i=1;i<=E->trace.ntracers[1];i++){
        index[i-1] = i;}
      //  fprintf(ft,"%d ",i);}
    write_int_binary_array(E->trace.ntracers[1],index,ft);
    //fputs("\n",ft);
    fputs("      </DataArray>\n",ft);
    fputs("    </Verts>\n",ft);
    fputs("  </Piece>\n",ft);

    fputs("  </PolyData>\n",ft);
    fputs("</VTKFile>\n",ft);
    fclose(ft);

    if(E->parallel.me == 0) write_pvtp(E, cycles);

}

static void write_ascii_array_float(int nn, int perLine, float *array, FILE *fp)
{
	int i;
		double* double_array = malloc(nn * sizeof(double));

		for (i=0;i<nn;i++) double_array[i] = (double) array [i];

		write_ascii_array(nn,perLine,4,double_array,fp);
		free(double_array);
}

static void write_ascii_array_double(int nn, int perLine, double *array, FILE *fp)
{
	write_ascii_array(nn,perLine,6,array,fp);
}


static void write_ascii_array(int nn, int perLine, int ascii_precision, double *array, FILE *fp)
{
    int i,a;
    char format_string[255];

    switch (perLine) {
    case 1:
        a = snprintf(format_string,sizeof(format_string), "%%.%de\n",
        		ascii_precision);
        for(i=0; i<nn; i++)
            fprintf(fp, format_string, array[i]);
        break;
    case 3:
        a = snprintf(format_string,sizeof(format_string), "%%.%de %%.%de %%.%de\n",
        		ascii_precision,ascii_precision,ascii_precision);
        for(i=0; i < nn/3; i++)
            fprintf(fp,format_string,array[3*i],array[3*i+1],array[3*i+2]);
        break;
    case 6:
        a = snprintf(format_string,sizeof(format_string), "%%.%de %%.%de %%.%de %%.%de %%.%de %%.%de\n",
        		ascii_precision,ascii_precision,ascii_precision,ascii_precision,ascii_precision,ascii_precision);
        for(i=0; i < nn/6; i++)
            fprintf(fp,format_string,
                    array[6*i],array[6*i+1],array[6*i+2],
                    array[6*i+3],array[6*i+4],array[6*i+5]);
        break;
    }
    return;
}

void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
{
    /* simple float to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union FloatToUnsignedChars
        {
            float input;
            unsigned char output[4];
        } floattransform;

    for (i=0; i<nn; i++){
        floattransform.input=floatarray[i];
        chararray[4*i]=floattransform.output[0];
        chararray[4*i+1]=floattransform.output[1];
        chararray[4*i+2]=floattransform.output[2];
        chararray[4*i+3]=floattransform.output[3];
    }
    return;
}

void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
{
    /* simple int - to unsigned chararray routine via union
    nn=length(intarray) chararray has to be BIG ENOUGH! */
    int i;
    union IntToUnsignedChars
        {
            int input;
            unsigned char output[4];
        } inttransform;

    for (i=0; i<nn; i++){
        inttransform.input=intarray[i];
        chararray[4*i]=inttransform.output[0];
        chararray[4*i+1]=inttransform.output[1];
        chararray[4*i+2]=inttransform.output[2];
        chararray[4*i+3]=inttransform.output[3];
        }
}


void zlibcompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
/* function to compress "in" to "out" reducing size from nn to nn2 */
{
    /* hope compressed data will be <= uncompressed */
    *out = malloc(sizeof(unsigned char)*nn);
    int ntemp=0; 

    /* in and out of z-stream */
    unsigned char inz[CHUNK];
    unsigned char outz[CHUNK];

    /* compression level */
    int level = Z_DEFAULT_COMPRESSION;
    int ret,flush;
    int i,j,k;

    /* zlib compression stream */
    z_stream strm;   

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    /* zlib init */
    ret = deflateInit(&strm, level);
    if (ret == Z_OK){
        i=0;     // position in "in" array
        do{
            j=0; // position in "inz"
            do{
                inz[j++]=in[i++];
            } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
            strm.avail_in=j;              // set number of input chars

            flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // done?
            strm.next_in = inz;           // set input buffer

            do{
                strm.avail_out = CHUNK;   // set number of max output chars
                strm.next_out = outz;     // set output buffer

                /* zlib compress */
                ret = deflate(&strm, flush);    
                assert(ret != Z_STREAM_ERROR);  

                /* zlib changed strm.avail_out=CHUNK
                 to the number of chars we can NOT use
                 in outz */

                for (k=0;k<CHUNK-strm.avail_out;k++){
                    (*out)[ntemp+k]=outz[k];
                }

                /* increase position in "out" */
                ntemp+=(CHUNK-strm.avail_out);
            }while(strm.avail_out==0);
            assert(strm.avail_in == 0);

        }while (flush != Z_FINISH);
    }
    else{fprintf(stderr,"Error during compression init\n");}

    // now we know how short "out" should be!
    *nn2=ntemp;
    *out = realloc(*out,sizeof(unsigned char)*ntemp);

    (void)deflateEnd(&strm);

    return;
}

void base64(unsigned char * in, int nn, unsigned char* out)
{
    /*takes *in*-array and "in"-length-"nn" and fills "out"-array 
    with base64(in) "out" needs to be big enough!!!
    length(out) >= 4* |^ nn/3.0 ^| */
    char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    int len;
    int i;

    for (i=0; i < nn; i+=3){

        len = (3 < nn-i ? 3 : nn-i);
        if (len >= 3){
        /* normal base64 encoding */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[ ((in[i+1] & 0x0f) << 2) | ((in[i+2] & 0xc0) >> 6)];
            out[i/3*4+3] = cb64[ in[i+2] & 0x3f ];
        } else if (len == 2){
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) | ((in[i+1] & 0xf0) >> 4) ];
            out[i/3*4+2] = cb64[((in[i+1] & 0x0f) << 2)];
            out[i/3*4+3] = (unsigned char) '=';
        } else if (len == 1){    
        /* at the end of array fill up with '=' */
            out[i/3*4+0] = cb64[ in[i] >> 2 ];
            out[i/3*4+1] = cb64[ ((in[i] & 0x03) << 4) ];
            out[i/3*4+2] = (unsigned char) '=';
            out[i/3*4+3] = (unsigned char) '=';
        }
    }
}


void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
{
    /* writing vtk compatible zlib compressed base64 encoded data to "out" */
    int i;
    /* header of data */
    unsigned char * charhead = malloc(sizeof(unsigned char)*16);
    /* - consists of "1" (number of pieces) */
    /* - original datalength in byte */
    /* - original datalength in byte */
    /* - new datalength after z-lib compression */
    int * headInts= malloc(sizeof(int)*4);
    headInts[0]=1;
    headInts[1]=orinn;
    headInts[2]=orinn;
    headInts[3]=nn;
    // transform to unsigned char
    IntToUnsignedChar(headInts,4,charhead);

    // base64: 16byte -> 24byte
    unsigned char * b64head =  malloc(sizeof(unsigned char)*24);
    // fills b64head
    base64(charhead, 16, b64head);

    // base64 data
    int b64bodylength = 4*ceil((double) nn/3.0);
    unsigned char * b64body = malloc(sizeof(unsigned char)*b64bodylength);
    // writes base64 data to b64body
    base64(in,nn,b64body);

    // combines header and body
    for (i=0; i<24 ; i++){
        out[i]=b64head[i];
    }

    for (i=0; i<b64bodylength ; i++){
        out[24+i]=b64body[i];
    }

    if(b64body){free(b64body);}
    if(b64head){free(b64head);}
    if(headInts){free(headInts);}
    if(charhead){free(charhead);}
}

void write_vtsarray(int nn, unsigned char * array, FILE * f)
{
    /*binary output routine*/
    int i;
    fwrite(array,sizeof(unsigned char),nn,f);
    fprintf (f,"\n");
}

static void write_int_binary_array(int nn, int* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
    IntToUnsignedChar(array,nn,chararray);

    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    int base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    unsigned char * base64plusheadarray= malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);
	
    write_vtsarray(base64plusheadlength, base64plusheadarray, f);
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

static void write_binary_array_double(int nn, double* array, FILE * f)
{
	int i;
	float* float_array = malloc(nn * sizeof(float));

	for (i=0;i<nn;i++) float_array[i] = (float) array [i];

	write_binary_array(nn,float_array,f);
	free(float_array);
}


static void write_binary_array(int nn, float* array, FILE * f)
{
    /* writes vtk-data array of floats and performs zip and base64 encoding */
    int chararraylength=4*nn;	/* nn floats -> 4*nn unsigned chars */
    unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
    FloatToUnsignedChar(array,nn,chararray);

    int compressedarraylength = 0;
    unsigned char * compressedarray;
    unsigned char ** pointertocompressedarray= &compressedarray;

    /* compression routine */
    zlibcompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

    /* special header for zip compressed and bas64 encoded data
    header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3) */
    int base64plusheadlength = 24 + 4*ceil((double) compressedarraylength/3.0);
    unsigned char * base64plusheadarray= malloc(sizeof(unsigned char)* base64plusheadlength);

    /* fills base64plusheadarray with everything ready for simple writing */
    base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);
	
    write_vtsarray(base64plusheadlength, base64plusheadarray, f);
    free(chararray);
    free(base64plusheadarray);
    free(compressedarray);
}

/**********************************************************************/

void vtk_output(struct All_variables *E, int cycles)
{
    char output_file[255];
    FILE *fp;
    int procs_per_cap;

    procs_per_cap = E->parallel.nprocx*E->parallel.nprocy*E->parallel.nprocz;

    snprintf(output_file, 255, "%s.proc%d.%d.vts",
             E->control.data_file, E->parallel.me, cycles);
    fp = output_open(output_file, "w");

    /* first, write volume data to vts file */
    vts_file_header(E, fp);

    /* write node-based field */
    vtk_point_data_header(E, fp);
    vtk_output_temp(E, fp);
    
    vtk_output_velo(E, fp);

    vtk_output_visc(E, fp);

    if (E->output.stress){
    vtk_output_stress(E, fp);}

    if (E->output.comp_nd && E->composition.on) {
    vtk_output_comp_nd(E, fp);}

    if (E->output.surf){
    vtk_output_surf_botm(E, fp, cycles);}

    if (E->output.density){
    vtk_output_dens(E, fp);}

    if (E->output.svelo){
    vtk_output_svelo(E, fp);}

    if (E->output.seismic){
    vtk_output_seismic(E, cycles, fp);}

    if (E->output.material){
    vtk_output_material(E, fp);}

    if (E->output.deltat){
    vtk_output_deltat(E,fp);}

    if (E->output.melttemp){
    vtk_output_melttemp(E,fp);}
    
    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);
    /**/
    if (E->output.tracer && E->composition.on){
    vtk_output_tracer(E, fp);}

    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);

    vts_file_trailer(E, fp);

    /* then, write other type of data */   

    fclose(fp);

    /* if processor is first of cap write summary for simple reading */
    if (E->parallel.me%procs_per_cap == 0) write_pvts(E, cycles);

    /* if processor is second write pvd for real time in vtk */
    if (E->parallel.me == 1%procs_per_cap) write_pvd(E, cycles);

    if ((E->output.write_tracer_file) && (E->monitor.elapsed_time >= E->output.tracer_file_time) && (E->output.tracer_file_written != 1)){
    	write_tracer_file(E, cycles);
    	E->output.tracer_file_written = 1;
    }

    return;
}
