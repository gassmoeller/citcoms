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
#include "material_properties_perplex.h"
#include "parallel_related.h"


static void read_refstate(struct All_variables *E);
static void read_continent_position(struct All_variables *E);
static void adams_williamson_eos(struct All_variables *E);
static void new_eos(struct All_variables *E);
const double get_alpha_nd_refstate(const struct All_variables *E, const int m, const int nn);
const double get_rho_nd_refstate(const struct All_variables *E, const int m, const int nn);
const double get_cp_nd_refstate(const struct All_variables *E, const int m, const int nn);
const double get_radheat_nd_refstate(const struct All_variables *E, const int m, const int nn);

int layers_r(struct All_variables *,float);


void allocate_refstate(struct All_variables *E)
{

    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;

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

    /* reference profile of adiabatic temperature */
    E->refstate.Tadi = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of initial temperature */
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

    if (E->composition.continents){
        E->refstate.cont_position = (int **) malloc((E->lmesh.nox+1)*sizeof(int *));
        int i;
        for (i=1;i<=E->lmesh.nox;i++){
            E->refstate.cont_position[i] = (int *) malloc((E->lmesh.noy+1)*sizeof(int));
        }
    }
}

void mat_prop_allocate(struct All_variables *E)
{
	allocate_refstate(E);
	if (E->refstate.choice == 3)
		allocate_perplex_refstate(E);
}

void set_refstate(struct All_variables *E)
{
    E->get_alpha_nd = get_alpha_nd_refstate;
    E->get_cp_nd = get_cp_nd_refstate;
    E->get_rho_nd = get_rho_nd_refstate;
    E->get_vp_nd = NULL;
    E->get_vs_nd = NULL;
    E->get_radheat_nd = get_radheat_nd_refstate;
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
    case 3:
        /* read from a file */
    	read_refstate(E);

        if (E->control.tracer == 0)
            E->perplex.nfields = 1;
        else
            E->perplex.nfields = max(1,E->trace.nflavors);

        read_perplexfile(E);
        break;
    case 4:
        read_refstate(E);

        if (E->control.tracer == 0)
            E->perplex.nfields = 1;
        else
            E->perplex.nfields = max(1,E->trace.nflavors);

        read_perplex_data(E);
        break;
    default:
        if (E->parallel.me) {
            fprintf(stderr, "Unknown option for reference state\n");
            fprintf(E->fp, "Unknown option for reference state\n");
            fflush(E->fp);
        }
        parallel_process_termination();
        break;
    }

    if ((E->refstate.choice == 3) || (E->refstate.choice == 4))
        set_perplex(E);
    else
        set_refstate(E);


    if(E->composition.continents && ((E->parallel.me+1) % E->parallel.nprocz == 0))
        read_continent_position(E);

    if(E->parallel.me == 0) {
        fprintf(stderr, "   nz     radius      depth    rho              layer\n");
    }
    if(E->parallel.me < E->parallel.nprocz)
        for(i=1; i<=E->lmesh.noz; i++) {
            fprintf(stderr, "%6d %11f %11f %11lf %5i\n",
                    i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                    E->refstate.rho[idxNz(i,E->lmesh.noz)],layers_r(E,E->sx[1][3][i]));
        }

    return;
}


static void read_refstate(struct All_variables *E)
{
    FILE *fp;
    int i,j;
    char buffer[255];
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
                    &(E->refstate.rho[j]),
                    &(E->refstate.gravity[j]),
                    &(E->refstate.thermal_expansivity[j]),
                    &(E->refstate.heat_capacity[j]),
                    &(E->refstate.Tadi[j]),
                    &(E->refstate.Tini[j]),
                    &(E->refstate.Tm[j]),
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
                        fprintf(stderr,"Error while reading continent file\n");
                        exit(8);
                    }
                }
            }
        }
    }

    fclose(fp);


    return;
}

static void adams_williamson_eos(struct All_variables *E)
{
    int i;
    double r, z, beta;

    if(E->parallel.me == 0) {
        fprintf(stderr,"Adams-Williamson EOS\n");
    }

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
        E->refstate.thermal_expansivity[i] = 1;
        E->refstate.heat_capacity[i] = 1;
        E->refstate.rho[i] = exp(beta*z);
    }

    return;
}

static void new_eos(struct All_variables *E)
{
    int i,j,k;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++)
    {
        r = E->sx[1][3][i];
        z = 1 - r;
        E->refstate.gravity[i] = 1;
        E->refstate.thermal_conductivity[i] = 1;
        E->refstate.free_enthalpy[i] = 1;
        E->refstate.rad_viscosity[i] = 1;
        E->refstate.stress_exp[i] = 1;
        /*E->refstate.Tadi[i] = (E->control.adiabaticT0 + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;*/
        E->refstate.Tadi[i] = 1;
        E->refstate.heat_capacity[i] = 1;
        E->refstate.thermal_expansivity[i] = 0.2 + 0.8 * (E->sx[1][3][i]- E->sphere.ri)/(E->sphere.ro - E->sphere.ri);
        E->refstate.rho[i] = exp(beta*z);
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

const double get_dimensionalT(const double dimensionlessT, const double dimensionlessT_surf, const double refT)
{
	return (dimensionlessT + dimensionlessT_surf) * refT;
}

/*
 * Returns the dimensional adiabatic correction Temperature. This is the temperature rise
 * that is missing  in the Temperature equation when using the Boussinesq approximation.
 */
const double get_adiabatic_correction(const struct All_variables *E, const int nz)
{
    if (E->control.disptn_number <= F_EPS)
        return E->refstate.Tadi[nz] - E->control.adiabaticT0 * E->data.ref_temperature;
    else
        return 0.0;
}

const double get_refTemp(const struct All_variables *E, const int m, const int nn, const int nz)
{
	const double compressible_factor = (E->control.disptn_number <= F_EPS) ? 0.0 : 1.0;
	const double compressible_correction = E->refstate.Tadi[nz]
			- E->control.adiabaticT0 * E->data.ref_temperature;

    double refTemp = get_dimensionalT(E->T[m][nn],E->control.surface_temp,E->data.ref_temperature)
    		- E->perplex.start_temp
    		+ get_adiabatic_correction(E,nz);

    refTemp = fmax(fmin(refTemp,E->perplex.end_temp-E->perplex.start_temp),0);
    return refTemp;
}

const int idxTemp(const double refTemp, const float delta_temp, const int ntempsteps)
{

	if (delta_temp == 0 || ntempsteps == 0)
		return -1;

	const int nT = (int) (refTemp / delta_temp + 1);
	const int bounded_nT = max(min(nT,ntempsteps-1),1);

	return bounded_nT;
}

const int idxNz (const int nn, const int noz) {
	if (noz != 0)
		return ((nn-1) % noz) + 1;
	else
	{
		fprintf(stderr,"Error in passing noz to idxNz(). noz must be larger than 0.");
		fflush(stderr);
		return -1;
	}
}

const double get_property_nd_refstate(const struct All_variables *E, const double* property, const int m, const int nn)
{
	return property[idxNz(nn,E->lmesh.noz)];
}

const double get_property_el(const struct All_variables *E, const double (*get_property_nd)(const struct All_variables *, const int, const int), const int m, const int el)
{
	int a,nn;

	const int ends=enodes[E->mesh.nsd];
	const int lev=E->mesh.levmax;
	double prop = 0;

	for(a=1;a<=ends;a++){
		nn = E->IEN[lev][m][el].node[a];
		prop += get_property_nd(E,m,nn);
	}

	prop /= ends;
	return prop;
}

const double get_cp_el(const struct All_variables *E, const int m, const int el)
{
    return get_property_el(E, E->get_cp_nd,m,el);
}

const double get_alpha_el(const struct All_variables *E, const int m, const int el)
{
    return get_property_el(E, E->get_alpha_nd,m,el);
}

const double get_rho_el(const struct All_variables *E, const int m, const int el)
{
    return get_property_el(E, E->get_rho_nd,m,el);
}

const double get_radheat_el(const struct All_variables *E, const int m, const int el)
{
    return get_property_el(E, E->get_radheat_nd,m,el);
}

const double get_cp_nd_refstate(const struct All_variables *E, const int m, const int nn)
{
    return get_property_nd_refstate(E, E->refstate.heat_capacity,m,nn);
}

const double get_alpha_nd_refstate(const struct All_variables *E, const int m, const int nn)
{
    return get_property_nd_refstate(E, E->refstate.thermal_expansivity,m,nn);
}

const double get_rho_nd_refstate(const struct All_variables *E, const int m, const int nn)
{
    return get_property_nd_refstate(E, E->refstate.rho,m,nn);
}

const double get_radheat_nd_refstate(const struct All_variables *E, const int m, const int nn)
{
    return E->control.Q0 * get_rho_nd_refstate(E,m,nn);
}
