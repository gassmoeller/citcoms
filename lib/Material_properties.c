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
#include "element_definitions.h"
#include "global_defs.h"
#include "material_properties.h"
#include "parallel_related.h"

static void read_refstate(struct All_variables *E);
static void read_perplexfile(struct All_variables *E);
static void read_densityfile(struct All_variables *E);
static void read_continent_position(struct All_variables *E);
static void adams_williamson_eos(struct All_variables *E);
static void new_eos(struct All_variables *E);

int layers_r(struct All_variables *,float);

void mat_prop_allocate(struct All_variables *E)
{
    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;
    int i,j;

    /* reference profile of density */
    E->refstate.rho = (double ***) malloc((E->composition.pressure_oversampling*(noz-1)+2)*sizeof(double **));
    for (i=1;i<=E->composition.pressure_oversampling*(noz-1)+1;i++){
        E->refstate.rho[i] = (double **) malloc((E->composition.ntdeps+1)*sizeof(double *));
        for (j=1;j<=E->composition.ntdeps;j++){
            E->refstate.rho[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
        }
    }

    /* reference profile of gravity */
    E->refstate.gravity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of coefficient of thermal expansion */
    E->refstate.thermal_expansivity = (double ***) malloc((E->composition.pressure_oversampling*(noz-1)+2)*sizeof(double **));
    for (i=1;i<=E->composition.pressure_oversampling*(noz-1)+1;i++){
        E->refstate.thermal_expansivity[i] = (double **) malloc((E->composition.ntdeps+1)*sizeof(double *));
        for (j=1;j<=E->composition.ntdeps;j++){
            E->refstate.thermal_expansivity[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
        }

    }

    /* reference profile of heat capacity */
    E->refstate.heat_capacity = (double ***) malloc((E->composition.pressure_oversampling*(noz-1)+2)*sizeof(double **));
    for (i=1;i<=E->composition.pressure_oversampling*(noz-1)+1;i++){
        E->refstate.heat_capacity[i] = (double **) malloc((E->composition.ntdeps+1)*sizeof(double *));
        for (j=1;j<=E->composition.ntdeps;j++){
            E->refstate.heat_capacity[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
        }
    }

    /* reference profile of thermal conductivity */
    E->refstate.thermal_conductivity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of temperature */
    E->refstate.Tadi = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of temperature */
    E->refstate.Tini = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of gravity */
    E->refstate.Tm = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of free enthalpy */
    /* only used in viscosity option 104 */
    E->refstate.free_enthalpy = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of radial viscosity */
    /* only used in viscosity option 104 */
    E->refstate.rad_viscosity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of stress dependence of viscosity */
    /* only used in viscosity option 104 */
    E->refstate.stress_exp = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of density contrast */
    /* only used in composition option zdep_buoyancy or tdep_buoyancy*/
    E->refstate.delta_rho = (double ***) malloc((E->composition.ncomp+2)*sizeof(double **));

    if (E->composition.zdep_buoyancy == 1){
        for(i=1;i<=E->composition.ncomp;i++){
            E->refstate.delta_rho[i] = (double **) malloc((noz+1)*sizeof(double *));
            if (E->composition.tdep_buoyancy == 1)
                for(j=1;j<=noz;j++)
                    E->refstate.delta_rho[i][j] = (double *) malloc(((E->composition.end_temp-E->composition.start_temp)/E->composition.delta_temp + 2)*sizeof(double));
            else 
                for(j=1;j<=noz;j++)
                    E->refstate.delta_rho[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
        }
    }

    if (E->composition.continents){
        E->refstate.cont_position = (int **) malloc((E->lmesh.nox+1)*sizeof(int *));
        for (i=1;i<=E->lmesh.nox;i++){
            E->refstate.cont_position[i] = (int *) malloc((E->lmesh.noy+1)*sizeof(int));
        }
    }

   /* reference profile of coefficient of seismic velocity */
    E->refstate.vp = (double ***) malloc((E->composition.pressure_oversampling*(noz-1)+2)*sizeof(double **));
    E->refstate.vs = (double ***) malloc((E->composition.pressure_oversampling*(noz-1)+2)*sizeof(double **));
    for (i=1;i<=E->composition.pressure_oversampling*(noz-1)+1;i++){
        E->refstate.vp[i] = (double **) malloc((E->composition.ntdeps+1)*sizeof(double *));
        E->refstate.vs[i] = (double **) malloc((E->composition.ntdeps+1)*sizeof(double *));
        for (j=1;j<=E->composition.ntdeps;j++){
            E->refstate.vp[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
            E->refstate.vs[i][j] = (double *) malloc((E->composition.ncomp+2)*sizeof(double));
        }

    }

}


void reference_state(struct All_variables *E)
{
    int i;

    /* All refstate variables (except Tadi) must be 1 at the surface.
     * Otherwise, the scaling of eqns in the code might not be correct. */

    /* select the choice of reference state */
    switch(E->refstate.choice) {
    case 0:
        /* read from a file */
        read_refstate(E);
        break;
    case 1:
        /* Adams-Williamson EoS */
        adams_williamson_eos(E);
        fprintf(stderr,"Adams-Williamson EOS\n");
        break;
    case 2:
        /* New EoS */
        new_eos(E);
        break;
    case 3:
        /* read from a file */
        read_perplexfile(E);
        break;
    default:
        if (E->parallel.me) {
            fprintf(stderr, "Unknown option for reference state\n");
            fprintf(E->fp, "Unknown option for reference state\n");
            fflush(E->fp);
        }
        parallel_process_termination();
    }


    if(E->composition.continents && ((E->parallel.me+1) % E->parallel.nprocz == 0))
        read_continent_position(E);

    if(E->parallel.me == 0) {
      fprintf(stderr, "   nz     radius      depth    rho              layer\n");
    }
    if(E->parallel.me < E->parallel.nprocz)
        for(i=1; i<=E->lmesh.noz; i++) {
            fprintf(stderr, "%6d %11f %11f %11lf %5i\n",
                    i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                    E->refstate.rho[(i-1)*E->composition.pressure_oversampling+1][E->composition.ntdeps/2][1],layers_r(E,E->sx[1][3][i]));
        }

    return;
}


static void read_refstate(struct All_variables *E)
{
    FILE *fp;
    int i;
    char buffer[255];
    double not_used1, not_used2, not_used3;

    fp = fopen(E->refstate.filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open reference state file: %s\n",
                E->refstate.filename);
        parallel_process_termination();
    }

    /* skip these lines, which belong to other processors */
    for(i=1; i<E->lmesh.nzs; i++) {
        fgets(buffer, 255, fp);
    }

    for(i=1; i<=E->lmesh.noz; i++) {
        fgets(buffer, 255, fp);
        if(sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  &(E->refstate.rho[i]),
                  &(E->refstate.gravity[i]),
                  &(E->refstate.thermal_expansivity[i]),
                  &(E->refstate.heat_capacity[i]),
                  &(E->refstate.Tadi[i]),
                  &(E->refstate.Tm[i]),
                  &(E->refstate.free_enthalpy[i]),
                  &(E->refstate.rad_viscosity[i]),
                  &(E->refstate.stress_exp[i]),
                  &(E->refstate.thermal_conductivity[i])) != 10) {
            fprintf(stderr,"Error while reading file '%s'\n", E->refstate.filename);
            exit(8);
        }
        /**** debug ****
        fprintf(stderr, "%d %f %f %f %f\n",
                i,
                E->refstate.rho[i],
                E->refstate.gravity[i],
                E->refstate.thermal_expansivity[i],
                E->refstate.heat_capacity[i]);
        /* end of debug */
    }

    fclose(fp);
    return;
}

static void read_perplexfile(struct All_variables *E)
{
    FILE *fp;
    int i,j,k;
    char buffer[255];
    char refstate_file[255];
    double not_used1, not_used2, not_used3;

    fp = fopen(E->refstate.filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open reference state file: %s\n",
                E->refstate.filename);
        parallel_process_termination();
    }
    j = 0;
    /* skip these lines, which belong to other processors */
    for(i=1; i<=(E->lmesh.nzs-1)*E->composition.pressure_oversampling; i++) {
        fgets(buffer, 255, fp);
    }
    for(i=0; i<E->composition.pressure_oversampling*(E->lmesh.noz-1)+1; i++){
        fgets(buffer, 255, fp);
        if (i%E->composition.pressure_oversampling == 0) {
        j++;
        if(sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  &(not_used1),
                  &(E->refstate.gravity[j]),
                  &(not_used2),
                  &(not_used3),
                  &(E->refstate.Tadi[j]),
                  &(E->refstate.Tini[j]),
                  &(E->refstate.Tm[i]),
                  &(E->refstate.free_enthalpy[j]),
                  &(E->refstate.rad_viscosity[j]),
                  &(E->refstate.stress_exp[j]),
                  &(E->refstate.thermal_conductivity[j])) != 11) {
            fprintf(stderr,"Error while reading file '%s'\n", E->refstate.filename);
            exit(8);
        }
        }
    }
    fclose(fp);

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
        fgets(buffer, 255, fp);
    }
    
    for(j=1; j<=E->composition.pressure_oversampling*(E->lmesh.noz-1)+1; j++) {
      for(k=1; k<=E->composition.ntdeps; k++){
          for(i=1; i<=E->composition.ncomp+1; i++){
            fgets(buffer, 255, fp);
            if(sscanf(buffer, "%lf %lf %lf %lf %lf\n",&(E->refstate.rho[j][k][i]),&(E->refstate.thermal_expansivity[j][k][i]), &(E->refstate.heat_capacity[j][k][i]), &(E->refstate.vp[j][k][i]), &(E->refstate.vs[j][k][i]))!=5){
                fprintf(stderr,"Error while reading file perplex.dat\n");
                exit(8);
	    }
		if(E->parallel.me == 0 || E->parallel.me == 1) fprintf(stderr, "me: %d noz:%d ntdeps:%d ncomp:%d rho:%f alpha:%f cp:%f\n",
                E->parallel.me,j,k,i,
                E->refstate.rho[j][k][i],
                E->refstate.thermal_expansivity[j][k][i],
                E->refstate.heat_capacity[j][k][i]);

}}}

    fclose(fp);
    return;
}


static void read_continent_position(struct All_variables *E)
{
    FILE *fp;
    int i,j,k;
    char buffer[255];

   sprintf(buffer,"continents.%d",E->sphere.capid[1]-1);
    fp = fopen(buffer, "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open continent file: %s\n",
                buffer);
        parallel_process_termination();
    }

        for(j=0; j<E->mesh.nox; j++) {
            for(k=0; k<E->mesh.noy; k++){
                fgets(buffer, 255, fp);
                if((j>=E->parallel.me_loc[1]*E->lmesh.elx) && (j<=(E->parallel.me_loc[1]+1)*E->lmesh.elx)){ // in x-Range
                    if((k>=E->parallel.me_loc[2]*E->lmesh.ely) && (k<=(E->parallel.me_loc[2]+1)*E->lmesh.ely)){
                        if(sscanf(buffer, "%d",&(E->refstate.cont_position[(j + E->parallel.me_loc[1])%E->lmesh.nox + 1][(k + E->parallel.me_loc[2])%E->lmesh.noy + 1]))!=1){
                            fprintf(stderr,"Error while reading file '%s'\n", E->refstate.densityfilename);
            exit(8);
        }}
    }}}

    fclose(fp);


    return;
}

static void adams_williamson_eos(struct All_variables *E)
{
    int i,j,k;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.gravity[i] = 1;
	E->refstate.thermal_conductivity[i] = 1;
        E->refstate.free_enthalpy[i] = 1;
        E->refstate.rad_viscosity[i] = 1;
        E->refstate.stress_exp[i] = 1;
	E->refstate.Tadi[i] = (E->control.TBCtopval + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;
        //E->refstate.Tadi[i] = 1;
        for(k=1;k<=E->composition.ntdeps;k++){
	    E->refstate.thermal_expansivity[i][k][1] = 1;
	    E->refstate.heat_capacity[i][k][1] = 1;
	    E->refstate.rho[i][k][1] = exp(beta*z);
          for(j=2;j<=E->composition.ncomp+1;j++){
	    E->refstate.thermal_expansivity[i][k][j] = 1;
	    E->refstate.heat_capacity[i][k][j] = 1;
	    E->refstate.rho[i][k][j] = exp(beta*z)*E->composition.buoyancy_ratio[j-2];
            if (E->parallel.me == 0) fprintf(stderr,"i: %d k: %d j: %d rho: %f\n",i,k,j,E->refstate.rho[i][k][j]);
            //E->refstate.delta_rho[j][i][k] = 1.0;}}
}}
    }

    return;
}

static void new_eos(struct All_variables *E)
{
    int i,j,k;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.gravity[i] = 1;
	E->refstate.thermal_conductivity[i] = 1;
        E->refstate.free_enthalpy[i] = 1;
        E->refstate.rad_viscosity[i] = 1;
        E->refstate.stress_exp[i] = 1;
	/*E->refstate.Tadi[i] = (E->control.adiabaticT0 + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;*/
        E->refstate.Tadi[i] = 1;
        for(k=1;k<=E->composition.ntdeps;k++){
            for(j=1;j<E->composition.ncomp+1;j++){
	      E->refstate.heat_capacity[i][k][j] = 1;
	      E->refstate.thermal_expansivity[i][k][j] = 0.2 + 0.8 * (E->sx[1][3][i]- E->sphere.ri)/(E->sphere.ro - E->sphere.ri);
	      E->refstate.rho[i][k][j] = exp(beta*z);
              E->refstate.delta_rho[j][i][k] = 1.0;}}
    }

    return;
}

double get_g_el(struct All_variables *E, int m, int el)
{
    int nz,nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double g = 0;

    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        nz = ((nn-1) % E->lmesh.noz) + 1;
        g += E->refstate.gravity[nz];
    }
    g /= ends;
    return g;
}

double get_cp_el(struct All_variables *E, int m, int el)
{
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double cp = 0;
    
    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        cp += get_cp_nd(E,m,nn);
    }
    cp /= ends;
    return cp;
}

/*double get_cp_el(struct All_variables *E, int m, int el)
{
    int i,nz,nzmin,nzmax,nT,j;
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double cp = 0;
    double meantemp = 0.0;
    
    for(a=1;a<=ends;a++)
        meantemp += E->T[m][E->IEN[lev][m][el].node[a]];
    meantemp /= ends;

    nn = E->IEN[lev][m][el].node[1];
    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((meantemp + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1,1);
    nzmax = min(E->composition.pressure_oversampling*nz + 1,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);
    for (i=nzmin;i<=nzmax;i++){
    cp += E->refstate.heat_capacity[i][nT][1];
    for(j=0;j<E->composition.ncomp;j++){
        cp +=  E->refstate.heat_capacity[i][nT][j+2]*E->composition.comp_node[m][j][nn];
    }
    }
    cp /= (nzmax-nzmin+1);

    return cp;
}*/

double get_cp_nd(struct All_variables *E, int m, int nn)
{
    int i,nz,nzmin,nzmax,nT,j;
    double cp;
    assert(E->composition.delta_temp > 0);

    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    cp = 0;
    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1-E->composition.pressure_oversampling/2,1);
    nzmax = min(E->composition.pressure_oversampling*(nz-1) + 1+E->composition.pressure_oversampling/2,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);
    for (i=nzmin;i<=nzmax;i++){
    cp += E->refstate.heat_capacity[i][nT][1];
    for(j=0;j<E->composition.ncomp;j++){
        cp +=  E->refstate.heat_capacity[i][nT][j+2]*E->composition.comp_node[m][j][nn];
    }
    }
    cp /= (nzmax-nzmin+1);

    return cp;
}

double get_rho_el(struct All_variables *E, int m, int el)
{
//    int i,nz,nzmin,nzmax,nT,j;
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double rho = 0;
//    double meantemp = 0.0;

    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        rho += get_rho_nd(E,m,nn);
    }
    rho /= ends;
    return rho;
}

    
/*    for(a=1;a<=ends;a++)
        meantemp += E->T[m][E->IEN[lev][m][el].node[a]];
    meantemp /= ends;

    nn = E->IEN[lev][m][el].node[1];
    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((meantemp + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1,1);
    nzmax = min(E->composition.pressure_oversampling*nz + 1,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);
    for (i=nzmin;i<=nzmax;i++){
    rho += E->refstate.rho[i][nT][1];
    for(j=0;j<E->composition.ncomp;j++){
        rho +=  E->refstate.rho[i][nT][j+2]*E->composition.comp_node[m][j][nn];
    }
    }
    rho /= (nzmax-nzmin+1);

    return rho;
}*/

double get_rho_nd(struct All_variables *E, int m, int nn)
{
    int i,nz,nzmax,nzmin,nT,j;
    double refTemp,rho,deltaT,weight;

    nz = ((nn-1) % E->lmesh.noz) + 1;
    refTemp = ((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature));
    refTemp = fmax(fmin(refTemp,E->composition.end_temp-E->composition.start_temp),0);
    nT = (int) (refTemp / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps-1),1);

    weight = fmax(fmin(refTemp / E->composition.delta_temp - (nT-1),1),0);

    rho = 0;

    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1 - E->composition.pressure_oversampling/2,1);
    nzmax = min(E->composition.pressure_oversampling*(nz-1) + 1 + E->composition.pressure_oversampling/2,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);


    for (i=nzmin;i<=nzmax;i++){
        rho += (1-weight) * E->refstate.rho[i][nT][1];
        rho += weight * E->refstate.rho[i][nT+1][1];

        for(j=0;j<E->composition.ncomp;j++){
            rho +=  (1-weight) * E->refstate.rho[i][nT][j+2]*E->composition.comp_node[m][j][nn];
            rho +=  weight * E->refstate.rho[i][nT+1][j+2]*E->composition.comp_node[m][j][nn];
        }
    }

    rho /= (nzmax-nzmin+1);

    return rho;
}

double get_alpha_el(struct All_variables *E, int m, int el)
{
    int nn,a;
//    int i,nz,nzmax,nzmin,nT,j;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double alpha = 0;
//    double meantemp = 0;
 

    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        alpha += get_alpha_nd(E,m,nn);
    }
    alpha /= ends;
    return alpha;
}

/*    for(a=1;a<=ends;a++)
      meantemp += E->T[m][E->IEN[lev][m][el].node[a]];
    meantemp /= ends;

    nn = E->IEN[lev][m][el].node[1];
    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((meantemp + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    alpha = 0;
    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1,1);
    nzmax = min(E->composition.pressure_oversampling*nz + 1,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);
    
    for (i=nzmin;i<=nzmax;i++){
    alpha += E->refstate.thermal_expansivity[i][nT][1];

    for(j=0;j<E->composition.ncomp;j++){
        alpha +=  E->refstate.thermal_expansivity[i][nT][j+2]*E->composition.comp_node[m][j][nn];
    }
    }
    alpha /= (nzmax-nzmin+1);

    return alpha;
}*/

double get_alpha_nd(struct All_variables *E, int m, int nn)
{
    int i,nz,nzmax,nzmin,nT,j;
    double alpha = 0.0;

    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    nzmin = max(E->composition.pressure_oversampling*(nz-1) + 1 - E->composition.pressure_oversampling/2,1);
    nzmax = min(E->composition.pressure_oversampling*(nz-1) + 1 + E->composition.pressure_oversampling/2,(E->lmesh.noz-1)*E->composition.pressure_oversampling+1);
    
    for (i=nzmin;i<=nzmax;i++){
    alpha += E->refstate.thermal_expansivity[i][nT][1];

    for(j=0;j<E->composition.ncomp;j++){
        alpha +=  E->refstate.thermal_expansivity[i][nT][j+2]*E->composition.comp_node[m][j][nn];
    }
    }
    alpha /= (nzmax-nzmin+1);

    return alpha;
}

double get_vs_el(struct All_variables *E, int m, int el)
{
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double vs = 0;
    
    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        vs += get_vs_nd(E,m,nn);
    }
    vs /= ends;
    return vs;
}

double get_vs_nd(struct All_variables *E, int m, int nn)
{
    int nz,nT,j;
    double vs;

    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    vs = E->refstate.vs[(nz-1)*E->composition.pressure_oversampling + 1][nT][1];

    for(j=0;j<E->composition.ncomp;j++){
        vs +=  E->refstate.vs[(nz-1)*E->composition.pressure_oversampling + 1][nT][j+2]*E->composition.comp_node[m][j][nn];
    }

    return vs;
}

double get_vp_el(struct All_variables *E, int m, int el)
{
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double vp = 0;
    
    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        vp += get_vp_nd(E,m,nn);
    }
    vp /= ends;
    return vp;
}

double get_vp_nd(struct All_variables *E, int m, int nn)
{
    int nz,nT,j;
    double vp;

    nz = ((nn-1) % E->lmesh.noz) + 1;
    nT = ((int)((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature)) / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps),1);

    vp = E->refstate.vp[(nz-1)*E->composition.pressure_oversampling + 1][nT][1];

    for(j=0;j<E->composition.ncomp;j++){
        vp +=  E->refstate.vp[(nz-1)*E->composition.pressure_oversampling + 1][nT][j+2]*E->composition.comp_node[m][j][nn];
    }

    return vp;
}

double get_radheat_el(struct All_variables *E, int m, int el)
{
    int nn,a;
    const int ends=enodes[E->mesh.nsd];
    const int lev=E->mesh.levmax;
    double radheat=0;
    
    for(a=1;a<=ends;a++){
        nn = E->IEN[lev][m][el].node[a];
        radheat += get_radheat_nd(E,m,nn);
    }
    radheat /= ends;
    return radheat;
}

double get_radheat_nd(struct All_variables *E, int m, int nn)
{
    int nz,nT,j;
    double refTemp,radheat,deltaT,weight;

    assert(E->composition.delta_temp > 0);

    nz = ((nn-1) % E->lmesh.noz) + 1;
    refTemp = ((E->T[m][nn] + E->control.surface_temp) * E->data.ref_temperature - E->composition.start_temp + (1-fmax(0,E->control.disptn_number)/fmax(1e-7,E->control.disptn_number))*(E->refstate.Tadi[nz]-E->control.adiabaticT0*E->data.ref_temperature));
    refTemp = fmax(fmin(refTemp,E->composition.end_temp-E->composition.start_temp),0);
    nT = (int) (refTemp / E->composition.delta_temp + 1);
    nT = max(min(nT,E->composition.ntdeps-1),1);

    weight = fmax(fmin(refTemp / E->composition.delta_temp - (nT-1),1),0);

    radheat = (1-weight) * E->refstate.rho[(nz-1)*E->composition.pressure_oversampling + 1][nT][1] * E->control.Q0;
    radheat += weight * E->refstate.rho[(nz-1)*E->composition.pressure_oversampling + 1][nT+1][1] * E->control.Q0;

    if (E->control.tracer_enriched){
        for(j=0;j<E->composition.ncomp;j++){
            radheat +=  (1-weight) * E->refstate.rho[(nz-1)*E->composition.pressure_oversampling + 1][nT][j+2]*E->composition.comp_node[m][j][nn] * E->control.Q0ER[j];
            radheat +=  weight * E->refstate.rho[(nz-1)*E->composition.pressure_oversampling + 1][nT+1][j+2]*E->composition.comp_node[m][j][nn] * E->control.Q0ER[j];
        }
    }

    return radheat;
}

