/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <string.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "material_properties.h"
#include "material_properties_perplex.h"
#include "parallel_related.h"
#include "output.h"

void allocate_perplex_refstate(struct All_variables *E)
{

    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;
    int i,j;
    const size_t nodes_size = (E->composition.pressure_oversampling
            * (noz - 1) + 2) * sizeof(double**);
    const size_t temperature_size = (E->composition.ntdeps + 1)
            * sizeof(double*);
    const size_t composition_size = (E->composition.ncomp + 2)
            * sizeof(double);

    /*
     * Perplex reference tables for nondimensional
     * density, thermal_expansivity, specific heat
     * and dimensional seismic velocities in dependence
     * of depth, temperature and composition
     */
    E->refstate.tab_density = (double***) malloc(nodes_size);
    E->refstate.tab_thermal_expansivity = (double ***) malloc(nodes_size);
    E->refstate.tab_heat_capacity = (double ***) malloc(nodes_size);
    E->refstate.tab_seismic_vp = (double ***) malloc(nodes_size);
    E->refstate.tab_seismic_vs = (double ***) malloc(nodes_size);

    for (i=1;i<=nodes_size-1;i++)
    {
        E->refstate.tab_density[i] = (double**) malloc(temperature_size);
        E->refstate.tab_thermal_expansivity[i] = (double **) malloc(temperature_size);
        E->refstate.tab_heat_capacity[i] = (double **) malloc(temperature_size);
        E->refstate.tab_seismic_vp[i] = (double **) malloc(temperature_size);
        E->refstate.tab_seismic_vs[i] = (double **) malloc(temperature_size);

        for (j=1;j<=temperature_size-1;j++)
        {
            E->refstate.tab_density[i][j] = (double*) malloc(composition_size);
            E->refstate.tab_thermal_expansivity[i][j] = (double *) malloc(composition_size);
            E->refstate.tab_heat_capacity[i][j] = (double *) malloc(composition_size);
            E->refstate.tab_seismic_vp[i][j] = (double *) malloc(composition_size);
            E->refstate.tab_seismic_vs[i][j] = (double *) malloc(composition_size);
        }
    }
}

void read_perplexfile(struct All_variables *E)
{
    FILE *fp = NULL;
    int i,j,k;
    char* buffer = malloc(255*sizeof(char));
    char refstate_file[255];
    double not_used1, not_used2, not_used3;

    if (E->parallel.me < E->parallel.nprocz){
        snprintf(refstate_file, 255, "%s.refstate.%d.csv", E->control.data_file, E->parallel.me);
        fp = output_open(refstate_file, "w");
        fprintf(fp,"Gravity AdiabaticTemperature InitialTemp SolidusTemp Activation_Enthalpy ViscosityPrefactor StressExponent ThermalDiffusivity\n");
        for (i=1;i <= E->lmesh.noz;i++){
            fprintf(fp,"%f %f %f %f %f %f %f\n", E->refstate.gravity[i]*E->data.grav_acc,E->refstate.Tadi[i],E->refstate.Tini[i],E->refstate.Tm[i],E->refstate.free_enthalpy[i],E->refstate.rad_viscosity[i],E->refstate.stress_exp[i],E->refstate.thermal_conductivity[i]*E->data.therm_diff);
        }
        fclose(fp);
    }

    /**** debug ****
    fprintf(stderr, "%d %f %f %f %f\n",
            i,
            E->refstate.rho[i],
            E->refstate.gravity[i],
            E->refstate.thermal_expansivity[i],
            E->refstate.heat_capacity[i]);
     end of debug */


    fp = fopen("perplex.dat", "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open perplex file: %s\n",
                "perplex.dat");
        parallel_process_termination();
    }
    /* skip these lines, which belong to other processors */
    for(i=1; i<=E->composition.pressure_oversampling*(E->lmesh.nzs-1)*E->composition.ntdeps*(E->composition.ncomp+1); i++) {
        buffer = fgets(buffer, 255, fp);
        if (buffer == NULL){
            fprintf(stderr, "Perplex file is too short or problem with processor number: %s\n",
                    "perplex.dat");
            parallel_process_termination();
        }
    }

    for(j=1; j<=E->composition.pressure_oversampling*(E->lmesh.noz-1)+1; j++) {
        for(k=1; k<=E->composition.ntdeps; k++){
            for(i=1; i<=E->composition.ncomp+1; i++){
                buffer = fgets(buffer, 255, fp);
                if(sscanf(buffer, "%lf %lf %lf %lf %lf\n",&(E->refstate.tab_density[j][k][i]),&(E->refstate.tab_thermal_expansivity[j][k][i]), &(E->refstate.tab_heat_capacity[j][k][i]), &(E->refstate.tab_seismic_vp[j][k][i]), &(E->refstate.tab_seismic_vs[j][k][i]))!=5){
                    fprintf(stderr,"Error while reading file perplex.dat\n");
                    exit(8);
                }

                /*
                 * Debug Output
                 */
                if(E->control.verbose)
                {
                    if(E->parallel.me == 0 || E->parallel.me == 1) fprintf(stderr, "me: %d noz:%d ntdeps:%d ncomp:%d rho:%f alpha:%f cp:%f\n",
                            E->parallel.me,j,k,i,
                            E->refstate.tab_density[j][k][i],
                            E->refstate.tab_thermal_expansivity[j][k][i],
                            E->refstate.tab_heat_capacity[j][k][i]);
                }
            }
        }
    }

    fclose(fp);
    return;
}

void read_error()
{
    fputs("Error while reading perplex file\n",stderr);
    parallel_process_termination();
}

/*
 * Returns the identifier of the field that is described by
 * field_name. Returns -1 if no field is recognized.
 */
int get_field_id(char* field_name)
{
    if (strstr(field_name,"T(K)") != NULL)
        return 0;
    else if (strstr(field_name,"P(bar)") != NULL)
        return 1;
    else if (strstr(field_name,"rho,kg/m3") != NULL)
        return 2;
    else if (strstr(field_name,"alpha,1/K") != NULL)
        return 3;
    else if (strstr(field_name,"cp,J/K/kg") != NULL)
        return 4;
    else if (strstr(field_name,"vp,km/s") != NULL)
        return 5;
    else if (strstr(field_name,"vs,km/s") != NULL)
        return 6;
    else if (strstr(field_name,"h,J/kg") != NULL)
        return 7;
    return -1;
}

/* Reads the order of the properties in the perplex data file
 * and saves the identifier of the ith field in input_fields_ids[i]
 * See get_field_id for description of identifiers
 */
struct IDs* read_field_order(char* bufline, int input_fields, struct IDs *field_ids)
{

    int i, field_id;
    int *input_field_ids = (int *) malloc(8*sizeof(int));
    memset(input_field_ids,-1,8*sizeof(int));

    char* field_name;

    field_name = strtok(bufline, " ");
    for (i=0;i<input_fields;i++)
    {
        field_id = get_field_id(field_name);

        if (field_id == -1)
        {
            fprintf(stderr, "Error in parsing field order of perplex file. Possibly data file"
                    "format has changed or file is corrupted?\n");
            parallel_process_termination();
        }
        input_field_ids[field_id] = i;
        field_name = strtok(NULL, " ");
    }

    for (i=3;i<6;i++)
        if (input_field_ids[i] == -1)
        {
            fprintf(stderr, "Error, necessary fields not found in perplex file. Possibly data file"
                    "format has changed or file is corrupted?\n");
            parallel_process_termination();
        }

    field_ids->rho = input_field_ids[2];
    field_ids->alpha = input_field_ids[3];
    field_ids->cp = input_field_ids[4];
    field_ids->vp = input_field_ids[5];
    field_ids->vs = input_field_ids[6];
    field_ids->h = input_field_ids[7];

    free(input_field_ids);

    return field_ids;
}

void read_perplex_header (FILE *perplex_file, struct table_properties *perplex_table)
{
    int i;
    const int buffer_size = 255;
    char* bufline = malloc(buffer_size * sizeof(char));

    // First lines, version of perplex and filename
    // Potential version check should be included here
    for (i = 0; i < 3; i++)
    {
        bufline = fgets (bufline, buffer_size, perplex_file);
    }

    // Read in independent variables Temperature, Pressure
    // in arbitrary order (Depends on Perplex input)
    for (i = 0; i < 2; i++)
    {
        bufline = fgets (bufline, buffer_size, perplex_file);
        if (strstr(bufline,"T(K)") != NULL)
        {
            perplex_table->field_ids.T = i;
        }
        else if (strstr(bufline,"P(bar)") != NULL)
        {
            perplex_table->field_ids.P = i;
        }
        else read_error();

        bufline = fgets (bufline, buffer_size, perplex_file);
        sscanf (bufline, "%lf\n", &(perplex_table->TP[i].start));

        bufline = fgets (bufline, buffer_size, perplex_file);
        sscanf (bufline, "%lf\n", &(perplex_table->TP[i].delta));

        bufline = fgets (bufline, buffer_size, perplex_file);
        sscanf (bufline, "%d\n", &(perplex_table->TP[i].ndeps));

        perplex_table->TP[i].end = perplex_table->TP[i].start + (perplex_table->TP[i].ndeps-1) * perplex_table->TP[i].delta;
    }

    // Read in number of fields
    bufline = fgets (bufline, buffer_size, perplex_file);
    sscanf (bufline, "%d\n",&(perplex_table->ninput_fields));

    // Read in fields and order
    bufline = fgets (bufline, buffer_size, perplex_file);
    perplex_table->field_ids = *(read_field_order(bufline, perplex_table->ninput_fields, &(perplex_table->field_ids)));

    free(bufline);
}

void read_perplex_body(double ****perplex_data, FILE *perplex_file, struct table_properties* perplex_table)
{
    int i,j,k;
    int buffer_size = 255;
    char* bufline = malloc(buffer_size * sizeof(char));
    char* bufval;

    for (i=0;i<perplex_table->TP[1].ndeps;i++)
        for (j=0;j<perplex_table->TP[0].ndeps;j++)
        {
            bufline = fgets (bufline, buffer_size, perplex_file);
            bufval = strtok(bufline, " ");

            for (k=0;k<perplex_table->ninput_fields;k++)
            {
                sscanf (bufval, "%le", &((*perplex_data)[i][j][k]));
                bufval = strtok(NULL, " ");
            }
        }
}

double ***allocate_perplex_data (struct table_properties *perplex_table, double *** perplex_data)
{
    int i,j;

    perplex_data = (double ***) malloc(sizeof(double**) * perplex_table->TP[1].ndeps);
    for(i=0;i<perplex_table->TP[1].ndeps;i++){
        perplex_data[i] = (double **) malloc(sizeof(double*) * perplex_table->TP[0].ndeps);
        for(j=0;j<perplex_table->TP[0].ndeps;j++){
            perplex_data[i][j] = (double *) malloc(sizeof(double) * perplex_table->ninput_fields);
        }
    }
    return perplex_data;
}

int idxPress (double refpress, double deltapress, int max)
{
    const int npress = min(max((int)(refpress/deltapress),0),max-1);
    return npress;
}

double get_depth_km (const struct All_variables *E, int k)
{
    double max_depth_km = E->data.radius_km * (E->sphere.ro - E->sphere.ri);

    if (k == 0) return max_depth_km;
    if (k == E->mesh.noz) return 0.0;

    int iz = ((k-1) / E->composition.pressure_oversampling) + 1;
    double weight =  ((double) ((k-1) % E->composition.pressure_oversampling) / (double) E->composition.pressure_oversampling);

    //fprintf(stderr,"weight: %f\n",weight);
    return  (E->sphere.ro - ((1-weight) * E->sx[1][3][iz] + weight * E->sx[1][3][iz+1])) * E->data.radius_km;
}

/*
 * Calculates the next step of a adiabatic profile. Assumes that node index k+1 already is computed
 * and computes the temperature and pressure of k according to the data present in perplex_data.
 * The dimensions and field ordering of perplex_data are described in perplex_table
 */
void calc_adiabatic_conditions(
        const struct table_properties *table,
        const double ***data,
        const int k,
        const double delta_depth_m,
        const double midpoint_gravity,
        double *Padi,
        double *Tadi)
{
    int TPidx[2];
    const struct IDs *id = &(table->field_ids);
    const struct field *TP = table->TP;

    // calc old indices
    TPidx[id->T] = idxTemp(Tadi[k+1]-TP[id->T].start,TP[id->T].delta,TP[id->T].ndeps);
    TPidx[id->P] = idxPress(Padi[k+1]-TP[id->P].start,TP[id->P].delta,TP[id->P].ndeps);

    // calc old properties
    double density = data[TPidx[0]][TPidx[1]][id->rho];
    double alpha = data[TPidx[0]][TPidx[1]][id->alpha];
    double cp = data[TPidx[0]][TPidx[1]][id->cp];

    //TODO: dh switch (see above)
   //double dHdp,dHdT;
    //dHdp = (perplex_data[pn][i][5] - perplex_data[pp][i][idH]) / (dp*TP[id.P]->delta*1e5);
    //dHdT = (perplex_data[idP][Tn][5] - perplex_data[idP][Tp][idH])/(dT*TP[id->T].start);
    //alpha = (1-data[idP][i][0]*dHdp)/((TP[id.T].start + i*TP[id->T].start)*E->data.therm_exp);
    //cp = dHdT/E->data.Cp;


        Tadi[k] = Tadi[k+1] * (1 + alpha * midpoint_gravity * delta_depth_m / cp);
        Padi[k]= Padi[k+1] + midpoint_gravity * delta_depth_m * density / 1e5;

        int i;
        double temperature;

        for(i=1;i<11;i++){
            TPidx[id->T] = idxTemp(Tadi[k+1]-TP[id->T].start,TP[id->T].delta,TP[id->T].ndeps);
            TPidx[id->P] = idxPress(Padi[k]-TP[id->P].start,TP[id->P].delta,TP[id->P].ndeps);

            density = 0.5*(density + data[TPidx[0]][TPidx[1]][id->rho]);
            alpha =  0.5*(alpha+data[TPidx[0]][TPidx[1]][id->alpha]);
            cp = 0.5*(cp+data[TPidx[0]][TPidx[1]][id->cp]);
            temperature = 0.5*(Tadi[k+1]+Tadi[k]);

            Tadi[k] = Tadi[k+1] + temperature * alpha * midpoint_gravity * delta_depth_m / cp;
            Padi[k]= Padi[k+1] + midpoint_gravity * delta_depth_m * density / 1e5;
        }

        printf("Tadi:%f pressure:%f k:%d density:%f delta_m:%f Tidx:%d Pidx:%d \n",Tadi[k],Padi[k],k, density, delta_depth_m,TPidx[id->T],TPidx[id->P]);
    }

void calculate_refstate_data (struct All_variables *E, struct table_properties *perplex_table, double ***perplex_data, const int idx_field)
{

    const struct IDs *id = &(perplex_table->field_ids);
    const struct field *TP = perplex_table->TP;


    //TODO: Input parameter to switch H calculation on/off
    //id->h = perplex_table->input_field_ids[5]

    const int limittoTm = 0;

    int TPidx[2];

    int maxnT;
    int i,k,l;
    double dHdp,dHdT;


    // Calculate adiabatic reference profile
    const int nodes = (E->mesh.noz-1) *E->composition.pressure_oversampling + 1;
    double Tadi[nodes];
    double Padi[nodes];

    Tadi[nodes] = E->control.adiabaticT0 * E->data.ref_temperature;
    Padi[nodes] = 1.0;
    double delta_depth_m,midpoint_gravity;

    for (k=nodes-1;k>0;k--){
        delta_depth_m = 1000*(get_depth_km(E,k-1)-get_depth_km(E,k));
        midpoint_gravity = 0.5 * (E->refstate.gravity[k+1]+E->refstate.gravity[k]) * E->data.grav_acc;
        calc_adiabatic_conditions(perplex_table,(const double ***) perplex_data,k,delta_depth_m,midpoint_gravity,Padi,Tadi);
    }


    // In local refstate range: fill refstate variables
    const int local_nodes = (E->lmesh.noz-1)*E->composition.pressure_oversampling + 1;
    const int start_node = E->composition.pressure_oversampling*(E->lmesh.nzs-1);
    const int end_node = start_node + local_nodes;
    for (k=start_node;k<end_node;k++)
    {
        if (limittoTm)
            maxnT = ((int) E->refstate.Tm[k] - TP[id->T].start)/TP[id->T].delta;
        else
            maxnT = TP[id->T].ndeps;

        int Tp,Tn,dT;
        int pp,pn,dp;

        TPidx[id->P] = idxPress(Padi[k],TP[id->P].delta,TP[id->P].ndeps);
        pp = max(0,TPidx[id->P]-1);
        pn = min(TP[id->P].ndeps-1,TPidx[id->P]+1);
        dp = pn-pp;

        for (l=1;l<=TP[id->T].ndeps;l++){

            i = min(l,maxnT);  // Limit index lower than T melt

            TPidx[id->T] = i;
            Tp = max(0,i-1);
            Tn = min(TP[id->T].ndeps-1,i+1);
            dT = Tn-Tp;
            fprintf(stderr,"Test\n");

            //TODO: dh switch (see above)
            //dHdp = (perplex_data[pn][i][5] - perplex_data[pp][i][idH]) / (dp*TP[id.P]->delta*1e5);
            //dHdT = (perplex_data[idP][Tn][5] - perplex_data[idP][Tp][idH])/(dT*delta_temp);
            //E->refstate[k][l][1] = (1-perplex_data[idP][i][0]*dHdp)/((start_temp + i*delta_temp)*E->data.therm_exp);
            //E->refstate[k][l][2] = dHdT/E->data.Cp;
            //E->refstate[k][l][5] = perplex_data[id[0]][id[1]][idH];

            E->refstate.tab_density[k][l][0] = perplex_data[TPidx[0]][TPidx[1]][id->rho]/E->data.density;
            E->refstate.tab_thermal_expansivity[k][l][1] = perplex_data[TPidx[0]][TPidx[1]][id->alpha]/E->data.therm_exp;
            E->refstate.tab_heat_capacity[k][l][2] = perplex_data[TPidx[0]][TPidx[1]][id->cp]/E->data.Cp;
            E->refstate.tab_seismic_vp[k][l][3] = perplex_data[TPidx[0]][TPidx[1]][id->vp];
            E->refstate.tab_seismic_vs[k][l][4] = perplex_data[TPidx[0]][TPidx[1]][id->vs];
        }
    }

// interpolate
//float borders[7][2] = {{2.0,0.8},{5.0,-5.0},{5.0,0.5},{15.0,4.5},{8.5,3.5},{5.0,-5.0},{5.0,0.5}};
//float borders[7][2] = {{2.0,0.8},{1.5,0.1},{1.5,0.5},{15.0,4.5},{8.5,3.5},{5.0,-5.0},{5.0,0.5}};
//border_field (E->refstate,borders,nodes,TP[id-T]->ndeps,5);
/*
    for (k=0;k<nodes;k++){
        for (i=0;i<numtemp;i++){
            for (l=0;l<5;l++){
                fix_value(E->refstate,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);
                fix_value(c_morb,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);
            }
        }
    }

    for (k=0;k<nodes;k++){
        for (i=0;i<numtemp;i++){
            for (l=0;l<5;l++){
                fix_vert_value(E->refstate,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);
                fix_vert_value(c_morb,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);
            }
        }
    }

    float eps = 1e-4;

    for (k=0;k<nodes;k++){
        for (i=0;i<numtemp;i++){
            for (l=0;l<5;l++){
                if ((E->refstate[k][i][l] >= borders[l][0] - eps) || (E->refstate[k][i][l] <= borders[l][1] + eps)) E->refstate[k][i][l] = interp_value(E->refstate,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);
                if ((c_morb[k][i][l] >= borders[l][0] - eps) || (c_morb[k][i][l] <= borders[l][1] + eps)) c_morb[k][i][l] = interp_value(c_morb,k,i,l,borders[l][0],borders[l][1],nodes,numtemp);

            }
        }
    }

border_field (E->refstate,borders,nodes,numtemp,5);
border_field (c_morb,borders,nodes,numtemp,5);
*/


// Filter
    int j,n,m,numpoints;


    float*** temp = (float ***) malloc(sizeof(float**) * nodes);
    for(i=0;i<nodes;i++){
        temp[i] = (float **) malloc(sizeof(float*) * TP[id->T].ndeps);
        for(k=0;k<TP[id->T].ndeps;k++){
            temp[i][k] = (float *) malloc(sizeof(float)*5);
        }
    }
/*
    for (i=0;i<nodes;i++){
        for (j=0;j<numtemp;j++){
            for(k=0;k<3;k++){
                c_morb_f[i][j][k]=0.0;
                E->refstate_f[i][j][k]=0.0;
                numpoints = 0;
                for (n=-1;n<=1;n++){
                    for (m=-5;m<=5;m++){
                        if ((i+n>=0) && (i+n<nodes) && (j+m>=0) && (j+m<numtemp)){
                            c_morb_f[i][j][k] += c_morb[i+n][j+m][k];
                            E->refstate_f[i][j][k] += E->refstate[i+n][j+m][k];
                            numpoints++;
                        }
                    }
                }
                c_morb_f[i][j][k]/=numpoints;
                E->refstate_f[i][j][k]/=numpoints;
            }
            E->refstate_f[i][j][3]=E->refstate[i][j][3];
            c_morb_f[i][j][3]=c_morb[i][j][3];
            E->refstate_f[i][j][4]=E->refstate[i][j][4];
            c_morb_f[i][j][4]=c_morb[i][j][4];
        }
    }
    for (i=0;i<nodes;i++){
        for (j=0;j<numtemp;j++){
            for(k=0;k<5;k++){
                c_morb[i][j][k] = c_morb_f[i][j][k];
                E->refstate[i][j][k] = E->refstate_f[i][j][k];
            }
        }
    }
*/
}

void read_perplex_data (struct All_variables *E, char* perplex_filename,int idx_field)
{

    struct table_properties *perplex_table = (struct table_properties*) malloc (sizeof(struct table_properties));
    double ****perplex_data = (double ****) malloc(sizeof(double***));

    fprintf (stderr, "Zero\n");

    FILE* perplex_file = fopen(perplex_filename,"r");

    // Get information about size of table, stored in perplex_table
    read_perplex_header(perplex_file,perplex_table);

    fprintf (stderr, "One\n");

    // Get perplex data stored in perplex_data
    *perplex_data = allocate_perplex_data(perplex_table, *perplex_data);

    fprintf (stderr, " OneandHalf \n");

    read_perplex_body(perplex_data,perplex_file,perplex_table);

    fclose(perplex_file);

    fprintf (stderr, "Two\n");

    // Transform to CitcomS data
    fprintf (stderr, "Three\n");


    calculate_refstate_data(E,perplex_table,*perplex_data,idx_field);

    fprintf (stderr, "Four\n");

}

const double get_property_nd_perplex(const struct All_variables *E, const double*** property, const int m, const int nn, const int temperature_accurate)
{
    int i,j;

    double prop = 0.0;
	double deltaT;

    const int nz = idxNz(nn, E->lmesh.noz);

    /* Calculating the borders in the case of pressure_oversampling == x > 1 (x times more depth nodes in material table than in model)
     * In case pressure_oversampling == 1 : nzmin == nzmax == nz
     */
    const int nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1 - E->composition.pressure_oversampling/2,1);
    const int nzmax = min(E->composition.pressure_oversampling*(nz-1) + 1 + E->composition.pressure_oversampling/2,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);

    const double refTemp = get_refTemp(E,m,nn,nz);
    const int nT = idxTemp(refTemp,E->composition.delta_temp,E->composition.ntdeps);

    const double weight = (temperature_accurate == 1) ? fmax(fmin(refTemp / E->composition.delta_temp - (nT-1),1),0) : 0.0;

    for (i=nzmin;i<=nzmax;i++){
    	prop = property[i][nT][1];
    	if (temperature_accurate == 1)
    	{
    		prop *= (1-weight);
    		prop += weight * property[i][nT+1][1];
    	}

		for(j=0;j<E->composition.ncomp;j++){
			prop +=  (1-weight) * property[i][nT][j+2]*E->composition.comp_node[m][j][nn];
			prop +=  weight * property[i][nT+1][j+2]*E->composition.comp_node[m][j][nn];
		}
    }
    prop /= (nzmax-nzmin+1);

	return prop;
}

const double get_radheat_nd_perplex(const struct All_variables *E, const int m,const int nn)
{

    const int nz = idxNz(nn,E->lmesh.noz);
    const double refTemp = get_refTemp(E,m,nn,nz);
    const int nT = idxTemp(refTemp,E->composition.delta_temp,E->composition.ntdeps);
    const double weight = fmax(fmin(refTemp / E->composition.delta_temp - (nT-1),1),0);

    double radheat = 0.0;
    double density = 0.0;
    double proportion_normal_material = 1.0;

    if (E->control.tracer_enriched){
	int j;
	for(j=0;j<E->composition.ncomp;j++){
	    proportion_normal_material -= E->composition.comp_node[m][j][nn];

	    density = E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT][j+2] + E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT][1];
	    radheat +=  (1-weight) * density * E->composition.comp_node[m][j][nn] * E->control.Q0ER[j];

	    density = E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT+1][j+2] + E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT+1][1];
	    radheat +=  weight * density * E->composition.comp_node[m][j][nn] * E->control.Q0ER[j];
	}
    }

    radheat += (1-weight) * E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT][1] * proportion_normal_material * E->control.Q0;
    radheat += weight * E->refstate.tab_density[(nz-1)*E->composition.pressure_oversampling + 1][nT+1][1] * proportion_normal_material * E->control.Q0;

    return radheat;
}

const double get_cp_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
    const int temperature_accurate = 0;
    return get_property_nd_perplex(E,(const double ***)E->refstate.tab_heat_capacity,m,nn,temperature_accurate);
}

const double get_alpha_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
    const int temperature_accurate = 0;
    return get_property_nd_perplex(E,(const double ***)E->refstate.tab_thermal_expansivity,m,nn,temperature_accurate);
}

const double get_rho_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
    const int temperature_accurate = 1;
    return get_property_nd_perplex(E,(const double ***)E->refstate.tab_density,m,nn,temperature_accurate);
}

const double get_vs_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
	const int temperature_accurate = 0;
        return get_property_nd_perplex(E,(const double ***)E->refstate.tab_seismic_vs,m,nn,temperature_accurate);
}

const double get_vp_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
	const int temperature_accurate = 0;
	return get_property_nd_perplex(E,(const double ***)E->refstate.tab_seismic_vp,m,nn,temperature_accurate);
}

