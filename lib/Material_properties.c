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

#include "global_defs.h"
#include "material_properties.h"
#include "parallel_related.h"

static void read_refstate(struct All_variables *E);
static void adams_williamson_eos(struct All_variables *E);
static void new_eos(struct All_variables *E);

int layers_r(struct All_variables *,float);

void mat_prop_allocate(struct All_variables *E)
{
    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;
    int i;

    /* reference profile of density */
    E->refstate.rho = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of gravity */
    E->refstate.gravity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of coefficient of thermal expansion */
    E->refstate.thermal_expansivity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of heat capacity */
    E->refstate.heat_capacity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of thermal conductivity */
    E->refstate.thermal_conductivity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of temperature */
    E->refstate.Tadi = (double *) malloc((noz+1)*sizeof(double));

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
    /* only used in tracer option */
    E->refstate.delta_rho = (double **) malloc((E->composition.ncomp+1)*sizeof(double*));
    for(i=1;i<=E->composition.ncomp;i++){
    E->refstate.delta_rho[i] = (double *) malloc((noz+1)*sizeof(double));}

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
        break;
    case 2:
        /* New EoS */
        new_eos(E);
        break;
    default:
        if (E->parallel.me) {
            fprintf(stderr, "Unknown option for reference state\n");
            fprintf(E->fp, "Unknown option for reference state\n");
            fflush(E->fp);
        }
        parallel_process_termination();
    }

    if(E->parallel.me == 0) {
      fprintf(stderr, "   nz     radius      depth    rho              layer\n");
    }
    if(E->parallel.me < E->parallel.nprocz)
        for(i=1; i<=E->lmesh.noz; i++) {
            fprintf(stderr, "%6d %11f %11f %11e %5i\n",
                    i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                    E->refstate.rho[i],layers_r(E,E->sx[1][3][i]));
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
        if(sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                  &(E->refstate.rho[i]),
                  &(E->refstate.gravity[i]),
                  &(E->refstate.thermal_expansivity[i]),
                  &(E->refstate.heat_capacity[i]),
                  &(E->refstate.Tadi[i]),
                  &(E->refstate.free_enthalpy[i]),
                  &(E->refstate.rad_viscosity[i]),
                  &(E->refstate.stress_exp[i]),
                  &(E->refstate.delta_rho[1][i]),
                  &(E->refstate.delta_rho[2][i]),
                  &(E->refstate.thermal_conductivity[i])) != 11) {
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


static void adams_williamson_eos(struct All_variables *E)
{
    int i,j;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.rho[i] = exp(beta*z);
	E->refstate.gravity[i] = 1;
	E->refstate.thermal_expansivity[i] = 1;
	E->refstate.heat_capacity[i] = 1;
	E->refstate.thermal_conductivity[i] = 1;
        E->refstate.free_enthalpy[i] = 1;
        E->refstate.rad_viscosity[i] = 1;
        E->refstate.stress_exp[i] = 1;
	E->refstate.Tadi[i] = (E->control.TBCtopval + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;
        //E->refstate.Tadi[i] = 1;
        for(j=1;j<E->composition.ncomp+1;j++){
            E->refstate.delta_rho[j][i] = 1.0;}

    }

    return;
}

static void new_eos(struct All_variables *E)
{
    int i;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.rho[i] = exp(beta*z);
	E->refstate.gravity[i] = 1;
	E->refstate.thermal_expansivity[i] = 0.2 + 0.8 * (E->sx[1][3][i]- E->sphere.ri)/(E->sphere.ro - E->sphere.ri);
	E->refstate.heat_capacity[i] = 1;
	E->refstate.thermal_conductivity[i] = 1;
        E->refstate.free_enthalpy[i] = 1;
        E->refstate.rad_viscosity[i] = 1;
        E->refstate.stress_exp[i] = 1;
	/*E->refstate.Tadi[i] = (E->control.adiabaticT0 + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;*/
        E->refstate.Tadi[i] = 1;
        for(j=1;j<E->composition.ncomp+1;j++){
            E->refstate.delta_rho[j][i] = 1.0;}
    }

    return;
}

