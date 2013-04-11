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


void mat_prop_allocate(struct All_variables *);
void reference_state(struct All_variables *);

const double get_cp_el(const struct All_variables *, const int m, const int el);
const double get_alpha_el(const struct All_variables *, const int m, const int el);
const double get_rho_el(const struct All_variables *, const int m, const int el);
double get_g_el(struct All_variables *, int m, int el);
double get_vp_el(struct All_variables *, int m, int el);
double get_vs_el(struct All_variables *, int m, int el);
const double get_radheat_el(const struct All_variables *, const int m, const int el);

const double get_adiabatic_correction(const struct All_variables *E, const int nz);
const double get_refTemp(const struct All_variables *, const int m, const int nn, const int nz);
const int idxTemp(const double refTemp, const float delta_temp, const int ntempsteps);
const int idxNz (const int nn, const int noz);



