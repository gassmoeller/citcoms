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

#if !defined(CitcomS_material_properties_perplex_h)
#define CitcomS_material_properties_perplex_h

#ifdef __cplusplus
extern "C" {
#endif

void allocate_perplex_refstate(struct All_variables *E);

void read_perplexfile(struct All_variables *E);

const double get_cp_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_alpha_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_rho_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_vp_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_vs_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_radheat_nd_perplex(const struct All_variables *E, const int m, const int nn);
const double get_adabatic_density_correction(const struct All_variables *E, const int m,const int nn);

struct IDs
    {
int T;
int P;
int rho;
int alpha;
int cp;
int vp;
int vs;
int h;
    };

struct field
    {
double start;
double end;
double delta;
int ndeps;
    };

struct table_properties
    {
int ninput_fields;

struct IDs field_ids;
struct field TP[2];
    };


#ifdef __cplusplus
}
#endif

#endif
