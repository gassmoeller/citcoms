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
    FILE *fp;
    int i,j,k;
    char buffer[255];
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
        fgets(buffer, 255, fp);
    }

    for(j=1; j<=E->composition.pressure_oversampling*(E->lmesh.noz-1)+1; j++) {
        for(k=1; k<=E->composition.ntdeps; k++){
            for(i=1; i<=E->composition.ncomp+1; i++){
                fgets(buffer, 255, fp);
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
    return get_property_nd_perplex(E,E->refstate.tab_heat_capacity,m,nn,temperature_accurate);
}

const double get_alpha_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
    const int temperature_accurate = 0;
    return get_property_nd_perplex(E,E->refstate.tab_thermal_expansivity,m,nn,temperature_accurate);
}

const double get_rho_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
    const int temperature_accurate = 1;
    return get_property_nd_perplex(E,E->refstate.tab_density,m,nn,temperature_accurate);
}

const double get_vs_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
	const int temperature_accurate = 0;
        return get_property_nd_perplex(E,E->refstate.tab_seismic_vs,m,nn,temperature_accurate);
}

const double get_vp_nd_perplex(const struct All_variables *E, const int m, const int nn)
{
	const int temperature_accurate = 0;
	return get_property_nd_perplex(E,E->refstate.tab_seismic_vp,m,nn,temperature_accurate);
}

