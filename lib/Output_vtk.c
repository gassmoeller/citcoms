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


#include <stdlib.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "output.h"

static void vts_file_header(struct All_variables *E, FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n"
        "  <StructuredGrid WholeExtent=\"%s\">\n"
        "    <Piece Extent=\"%s\">\n";

    char extent[64], header[1024];

    snprintf(extent, 64, "%d %d %d %d %d %d",
             E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz,
             E->lmesh.exs, E->lmesh.exs + E->lmesh.elx,
             E->lmesh.eys, E->lmesh.eys + E->lmesh.ely);

//           E->lmesh.exs, E->lmesh.exs + E->lmesh.elx,
//           E->lmesh.eys, E->lmesh.eys + E->lmesh.ely,
//           E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz);

    snprintf(header, 1024, format, extent, extent);

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
    int i, j;

    fputs("        <DataArray type=\"Float32\" Name=\"temperature\" format=\"ascii\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            fprintf(fp, "%.6e\n", E->T[j][i]);
        }
    }

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_velo(struct All_variables *E, FILE *fp)
{
    int i, j;
    double sint, sinf, cost, cosf;
    float *V[4];
    const int lev = E->mesh.levmax;

    fputs("        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        V[1] = E->sphere.cap[j].V[1];
        V[2] = E->sphere.cap[j].V[2];
        V[3] = E->sphere.cap[j].V[3];

        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            fprintf(fp, "%.6e %.6e %.6e\n",
                    V[1][i]*cost*cosf - V[2][i]*sinf + V[3][i]*sint*cosf,
                    V[1][i]*cost*sinf + V[2][i]*cosf + V[3][i]*sint*sinf,
                    -V[1][i]*sint + V[3][i]*cost);
        }
    }

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_visc(struct All_variables *E, FILE *fp)
{
    int i, j;
    int lev = E->mesh.levmax;

    fputs("        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"ascii\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++)
            fprintf(fp, "%.4e\n", E->VI[lev][j][i]);
    }

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_coord(struct All_variables *E, FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, j;

    fputs("      <Points>\n", fp);
    fputs("        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++)
            fprintf(fp,"%.6e %.6e %.6e\n",
                    E->x[j][1][i],
                    E->x[j][2][i],
                    E->x[j][3][i]);
    }

    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);

    return;
}

static void vtk_output_stress(struct All_variables *E, FILE *fp)
{
  int m, node;
 /* for stress computation */
  void allocate_STD_mem();
  void compute_nodal_stress();
  void free_STD_mem();
  float *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
  float *divv[NCS],*vorv[NCS];
  /*  */
    allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);

fputs("        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"ascii\">\n", fp);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    /* those are sorted like stt spp srr stp str srp  */
    for (node=1;node<=E->lmesh.nno;node++) {   
      fprintf(fp, "%.4e %.4e %.4e %.4e %.4e %.4e\n",
              E->gstress[m][(node-1)*6+1],
              E->gstress[m][(node-1)*6+2],
              E->gstress[m][(node-1)*6+3],
              E->gstress[m][(node-1)*6+4],
              E->gstress[m][(node-1)*6+5],
              E->gstress[m][(node-1)*6+6]);
    }
  }

fputs("        </DataArray>\n", fp);
return;
}

void output_vtk_stress(struct All_variables *E, int cycles)
{
  int m, node;
  char output_file[255];
  FILE *fp1;
 /* for stress computation */
  void allocate_STD_mem();
  void compute_nodal_stress();
  void free_STD_mem();
  float *SXX[NCS],*SYY[NCS],*SXY[NCS],*SXZ[NCS],*SZY[NCS],*SZZ[NCS];
  float *divv[NCS],*vorv[NCS];
  /*  */
    allocate_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    compute_nodal_stress(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
    free_STD_mem(E, SXX, SYY, SZZ, SXY, SXZ, SZY, divv, vorv);
  sprintf(output_file,"%s.stress.%d.%d", E->control.data_file,
          E->parallel.me, cycles);
  fp1 = output_open(output_file, "w");

  fprintf(fp1,"%d %d %.5e\n",cycles,E->lmesh.nno,E->monitor.elapsed_time);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    fprintf(fp1,"%3d %7d\n",m,E->lmesh.nno);
    /* those are sorted like stt spp srr stp str srp  */
    for (node=1;node<=E->lmesh.nno;node++)
      fprintf(fp1, "%.4e %.4e %.4e %.4e %.4e %.4e\n",
              E->gstress[m][(node-1)*6+1],
              E->gstress[m][(node-1)*6+2],
              E->gstress[m][(node-1)*6+3],
              E->gstress[m][(node-1)*6+4],
              E->gstress[m][(node-1)*6+5],
              E->gstress[m][(node-1)*6+6]);
  }
  fclose(fp1);
}



static void vtk_output_comp_nd(struct All_variables *E, FILE *fp)
{
    int i, j, k;
    char name[255];

    for(k=0;k<E->composition.ncomp;k++) {

    snprintf(name, 255, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"ascii\">\n", (k+1));

//    name1 = "composition";
//    itoa(i,name2,10);
//    strcat(name1,name2);

//    name2 = "        <DataArray type=\"Float32\" Name=\"";
//    strcat(name2, name1);
//    strcat(name2, "\" format=\"ascii\">\n");

    fputs(name, fp); 

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            fprintf(fp,"%.6e\n",E->composition.comp_node[j][k][i]);
            }
        }

    fputs("        </DataArray>\n", fp);

    }

    return;
}


void vtk_output_surf(struct All_variables *E,  FILE *fp, int cycles)
{
  int i, j, k, s;
  char output_file[255];
  float *topo;

  if((E->output.write_q_files == 0) || (cycles == 0) ||
     (cycles % E->output.write_q_files)!=0)
      heat_flux(E);
  /* else, the heat flux will have been computed already */

  if(E->control.use_cbf_topo){
    get_CBF_topo(E,E->slice.tpg,E->slice.tpgb);

  }else{
    get_STD_topo(E,E->slice.tpg,E->slice.tpgb,E->slice.divg,E->slice.vort,cycles);
  }

  if (E->output.surf && (E->parallel.me_loc[3]==E->parallel.nprocz-1)) {

    fputs("        <DataArray type=\"Float32\" Name=\"surface\" format=\"ascii\">\n", fp);

    for(j=1;j<=E->sphere.caps_per_proc;j++)  {
        /* choose either STD topo or pseudo-free-surf topo */
        if(E->control.pseudo_free_surf)
            topo = E->slice.freesurf[j];
        else
            topo = E->slice.tpg[j];


        for(i=1;i<=E->lmesh.nsf;i++)   {
          for(k=1;k<=E->lmesh.noz-1;k++)   
            fprintf(fp, "%.4e\n", 0.0);
          fprintf(fp, "%.4e\n", topo[i]);
        }

    }

    fputs("        </DataArray>\n", fp);

  }

  return;
}


/**********************************************************************/

void vtk_output(struct All_variables *E, int cycles)
{
    char output_file[255];
    FILE *fp;

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

    vtk_output_surf(E, fp, cycles); 
    
    vtk_point_data_trailer(E, fp);

    /* write element-based field */
    vtk_cell_data_header(E, fp);
    /**/
    vtk_cell_data_trailer(E, fp);

    /* write coordinate */
    vtk_output_coord(E, fp);

    vts_file_trailer(E, fp);

    /* then, write other type of data */   


    fclose(fp);


    return;
}
