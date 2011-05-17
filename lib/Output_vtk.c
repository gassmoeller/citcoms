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

#include "zlib.h"
#include <assert.h>

#define CHUNK 16384

static void write_array(int nn, float * array, FILE * f);
void DoubleToFloatArray(double* doublearray, int nn, float* floatarray);

static void vts_file_header(struct All_variables *E, FILE *fp)
{

    const char format[] =
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"StructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
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
    int lev = E->mesh.levmax;
    int array_size = E->sphere.caps_per_proc*E->lmesh.nno;
    float* floattemp = malloc(sizeof(E->T[1][1])*array_size);
    fputs("        <DataArray type=\"Float32\" Name=\"temperature\" format=\"binary\">\n", fp);

/*    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            fprintf(fp, "%.6e\n", E->T[j][i]);
        }
    }
*/
    DoubleToFloatArray(&E->T[1][1], array_size, floattemp);
    write_array(array_size,floattemp,fp);

    fputs("        </DataArray>\n", fp);
    free(floattemp);
    return;
}


static void vtk_output_velo(struct All_variables *E, FILE *fp)
{
    int i, j;
    double sint, sinf, cost, cosf;
    float *V[4];
    const int lev = E->mesh.levmax;
    float* tempvel = malloc(sizeof(float)*E->sphere.caps_per_proc*E->lmesh.nno*3);
    fputs("        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"binary\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        V[1] = E->sphere.cap[j].V[1];
        V[2] = E->sphere.cap[j].V[2];
        V[3] = E->sphere.cap[j].V[3];

        for(i=1; i<=E->lmesh.nno; i++) {
            sint = E->SinCos[lev][j][0][i];
            sinf = E->SinCos[lev][j][1][i];
            cost = E->SinCos[lev][j][2][i];
            cosf = E->SinCos[lev][j][3][i];

            tempvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3 +0] = (float)(V[1][i]*cost*cosf - V[2][i]*sinf + V[3][i]*sint*cosf);
            tempvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3 +1] = (float)(V[1][i]*cost*sinf + V[2][i]*cosf + V[3][i]*sint*sinf);
            tempvel[(((j-1)*E->sphere.caps_per_proc)+i-1)*3 +2] = (float)(-V[1][i]*sint + V[3][i]*cost);
        }
    }
    write_array(E->sphere.caps_per_proc*E->lmesh.nno*3,tempvel,fp);
    fputs("        </DataArray>\n", fp);
    free(tempvel);
    return;
}


static void vtk_output_visc(struct All_variables *E, FILE *fp)
{
    int i, j;
    int lev = E->mesh.levmax;

    fputs("        <DataArray type=\"Float32\" Name=\"viscosity\" format=\"binary\">\n", fp);

/*    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++)
            fprintf(fp, "%.4e\n", E->VI[lev][j][i]);
    }
*/
    
    write_array((int)(E->sphere.caps_per_proc*E->lmesh.nno),&E->VI[lev][1][1],fp);

    fputs("        </DataArray>\n", fp);
    return;
}


static void vtk_output_coord(struct All_variables *E, FILE *fp)
{
    /* Output Cartesian coordinates as most VTK visualization softwares
       assume it. */
    int i, j;
    float* pos = malloc(sizeof(float)*E->sphere.caps_per_proc*E->lmesh.nno*3);

    fputs("      <Points>\n", fp);
    fputs("        <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\">\n", fp);

    for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++){
                   pos[((j-1)*E->sphere.caps_per_proc+i-1)*3] = (float)(E->x[j][1][i]);
		   pos[((j-1)*E->sphere.caps_per_proc+i-1)*3+1]=(float)(E->x[j][2][i]);
		   pos[((j-1)*E->sphere.caps_per_proc+i-1)*3+2]=(float)(E->x[j][3][i]);
    	}
    }
    write_array(E->sphere.caps_per_proc*E->lmesh.nno*3,pos,fp);
    fputs("        </DataArray>\n", fp);
    fputs("      </Points>\n", fp);
    free(pos);
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

fputs("        <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"binary\">\n", fp);

/*  for(m=1;m<=E->sphere.caps_per_proc;m++) {*/
    /* those are sorted like stt spp srr stp str srp  */
    /*for (node=1;node<=E->lmesh.nno;node++) {   
      fprintf(fp, "%.4e %.4e %.4e %.4e %.4e %.4e\n",
              E->gstress[m][(node-1)*6+1],
              E->gstress[m][(node-1)*6+2],
              E->gstress[m][(node-1)*6+3],
              E->gstress[m][(node-1)*6+4],
              E->gstress[m][(node-1)*6+5],
              E->gstress[m][(node-1)*6+6]);
    }
  }*/
    write_array(E->sphere.caps_per_proc*E->lmesh.nno*6,&E->gstress[1][1],fp);

fputs("        </DataArray>\n", fp);
return;
}

static void vtk_output_comp_nd(struct All_variables *E, FILE *fp)
{
    int i, j, k;
    char name[255];
    int nodes = E->sphere.caps_per_proc*E->lmesh.nno;
    float* compo[nodes];
    for(k=0;k<E->composition.ncomp;k++) {

    snprintf(name, 255, "        <DataArray type=\"Float32\" Name=\"composition%d\" format=\"binary\">\n", (k+1));

    fputs(name, fp); 
    
    DoubleToFloatArray(E->composition.comp_node, nodes, compo);
    /*for(j=1; j<=E->sphere.caps_per_proc; j++) {
        for(i=1; i<=E->lmesh.nno; i++) {
            fprintf(fp,"%.6e\n",E->composition.comp_node[j][k][i]);
            }
        }*/
    write_array(nodes,compo,fp);

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

void FloatToUnsignedChar(float * floatarray, int nn, unsigned char * chararray)
// simple float - to unsigned chararray routine via union
// nn=length(intarray)
// chararray has to be BIG ENOUGH!
{
	int i;
        union FloatToUnsignedChars      // will handle transformation float -> 4 unsigned char
        {
                float input;
                unsigned char output[ 4 ];
        } floattransform;

        for (i=0; i<nn; i++)
        {
                floattransform.input=floatarray[i];
                chararray[4*i]=floattransform.output[0];
                chararray[4*i+1]=floattransform.output[1];
                chararray[4*i+2]=floattransform.output[2];
                chararray[4*i+3]=floattransform.output[3];
        }
}

void DoubleToFloatArray(double* doublearray, int nn, float* floatarray)
/*converts array of double values to floatvalues */
{
	int i;
	
	for (i=0; i<nn; i++)
	{
		floatarray[i]=(float)(doublearray[i]);
	}
}

void IntToUnsignedChar(int * intarray, int nn, unsigned char * chararray)
// simple int - to unsigned chararray routine via union
// nn=length(intarray)
// chararray has to be BIG ENOUGH!
{
        int i;
	union IntToUnsignedChars        // will handle transformation int -> 4 unsigned char
        {
                int input;
                unsigned char output[ 4 ];
        } inttransform;

        for (i=0; i<nn; i++)
        {
                inttransform.input=intarray[i];
                chararray[4*i]=inttransform.output[0];
                chararray[4*i+1]=inttransform.output[1];
                chararray[4*i+2]=inttransform.output[2];
                chararray[4*i+3]=inttransform.output[3];
        }
}


void usercompress(unsigned char* in, int nn, unsigned char** out, int *nn2)
// compression routine
{
        // creates temporarily output array
        // hope compressed data will be <= uncompressed
        // otherwise...  ??

        unsigned char* outtemp = malloc(sizeof(unsigned char)*nn);
        int ntemp=0; // position in "outtemp"

        // in and out of z-stream
        unsigned char inz[CHUNK];
        unsigned char outz[CHUNK];

        // compression level
        // could be better, don't know the options
        int level = Z_DEFAULT_COMPRESSION;
        int ret,flush;
	int i;
        // the compression stream
        z_stream strm;

        // nothing to see here
        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;

        // maschine is running?
        ret = deflateInit(&strm, level);
        if (ret == Z_OK)
        {
                // ok lets go

                i=0;       // will store position in "in" array
                int j=0;       // will store position in "inz" array
		int k;
                do
                {
                        j=0; // new filling of "inz"
                        do
                        {
                                inz[j++]=in[i++];
                        } while((j<CHUNK) && (i<nn)); // stopps if "inz"-buffer is full or "in" array empty
                        strm.avail_in=j;        // zlib can use j chars in "inz"

                        flush = (i==nn) ? Z_FINISH : Z_NO_FLUSH; // will tell zlib if we are done
                        strm.next_in = inz; // data array, zlib should use

                        do
                        {
                                strm.avail_out = CHUNK; // zlib can write "CHUNK" chars
                                strm.next_out = outz;   // in "outz"

                                // zlib compression
				ret = deflate(&strm, flush);    /* no bad return value */
                                assert(ret != Z_STREAM_ERROR);  /* state not clobbered */

                                // zlib changed strm.avail_out=CHUNK
                                // to the number of chars we can NOT use
                                // in outz

                                // saving good part of "outz" to "outtemp"
                                for (k=0;k<CHUNK-strm.avail_out;k++)
                                        {outtemp[ntemp+k]=outz[k];}
                                // have to increase position in "outtemp"
                                ntemp+=(CHUNK-strm.avail_out);
                        } while(strm.avail_out==0);
                        assert(strm.avail_in == 0);

                } while (flush != Z_FINISH);
                // DONE
        }
        //else{cout << "ARRG! Problem in Compression!"<<endl;}
	else{}
        // now we know how short "out" should be!
        *nn2=ntemp;
        // and can create an array with the right size
        *out = malloc(sizeof(unsigned char)*ntemp);

        // copying (memcpy could be better, anyway)
        for (i=0;i<ntemp; i++){(*out)[i]=outtemp[i];}

        (void)deflateEnd(&strm);

        if(outtemp){free(outtemp);};
	
//      old routine for changing nothing
//      *nn2=nn;
//      *out = new unsigned char[nn];
//      for (int i(0);i<nn;i++){(*out)[i]=in[i];} 
}

void base64(unsigned char * in, int nn, unsigned char* out)
// takes *in*-array and "in"-length-"nn"
// and fills "out"-array with base64(in) 
// "out" needs to be big enough!!!
// length(out) >= 4* |^ nn/3.0 ^| 
{
        char cb64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        unsigned char * inblock = malloc(sizeof(unsigned char)*3);        // input array for base64
        //unsigned char * outblock = new unsigned char [4];     // output array for base64      
        int len=3;
	int i;
        // idee:
        // repeat       3 chars -> 4 chars      with len=3
        // if not enough chars left, flood with 0'
        // and set len to (1,2)

        for (i=0; i < nn; i+=3)
        {
                inblock[0]=in[i];
                if ( i+1 < nn)
                {
                        inblock[1]=in[i+1];
                        if ( i+2 < nn)
                        {
                                inblock[2]=in[i+2];
                                len=3;
                        }
                        else
                        {
                                len=2;
                                inblock[2]=0;
                        }
                }
                else
                {
                        len=1;
                        inblock[1]=0; inblock[2]=0;
                }
                //normal base64 endcoding
                out[i/3*4] = cb64[ inblock[0] >> 2 ];
                out[i/3*4+1] = cb64[ ((inblock[0] & 0x03) << 4) | ((inblock[1] & 0xf0) >> 4) ];
                out[i/3*4+2] = (unsigned char) (len > 1 ? cb64[ ((inblock[1] & 0x0f) << 2) | ((inblock[2] & 0xc0) >> 6) ] : '=');
                out[i/3*4+3] = (unsigned char) (len > 2 ? cb64[ inblock[2] & 0x3f ] : '=');
        }
	free(inblock);
}


void base64plushead(unsigned char * in, int nn, int orinn, unsigned char* out)
// writing vtk compatible zlib compressed base64 encoded data to "out"
{
	int i;
        // cout << "Vor C: " <<  orinn << " nach C: " << nn << endl;

        // header of data
        unsigned char * charhead = malloc(sizeof(unsigned char)*16);
        // - consists of "1" (number of pieces)
        // - original datalength in byte
        // - original datalength in byte
        // - new datalength after z-lib compression
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
        int b64bodylength = 4*ceil((float) nn/3.0);
        unsigned char * b64body = malloc(sizeof(unsigned char)*b64bodylength);
        // writes base64 data to b64body
        base64(in,nn,b64body);

        // combines header and body
        for (i=0; i<24 ; i++)
        {
                out[i]=b64head[i];
        }

        for (i=0; i<b64bodylength ; i++)
        {
                out[24+i]=b64body[i];
        }

        if(b64body){free(b64body);}
	if(b64head){free(b64head);}
	if(headInts){free(headInts);}
	if(charhead){free(charhead);}
}

void write_vtsarray(int nn, unsigned char * array, FILE * f)
// binary output routine, just writing no thinking
{
	int i;
        for(i=0; i < nn; i++)
        {
                fprintf(f,"%c",(char) array[i]);
        }
        fprintf (f,"\n");
}


static void write_array(int nn, float* array, FILE * f)
/* writes a typical vtk-data array of floats in "f" and performs zip and base64 encoding if (bin) */
{
	int chararraylength=4*nn;	/*nn floats -> 4*nn unsigned chars*/
	unsigned char * chararray = malloc (chararraylength * sizeof(unsigned char));
	FloatToUnsignedChar(array,nn,chararray);

	int compressedarraylength = 0;
	unsigned char * compressedarray;
	unsigned char ** pointertocompressedarray= &compressedarray;

	/*compression routine: creates compressedarray and gives the new length of data in compressedarraylength*/
	usercompress(chararray,chararraylength,pointertocompressedarray,&compressedarraylength);

	/*special header for zip compressed and bas64 encoded data
	header needs 4 int32 = 16 byte -> 24 byte due to base64 (4*16/3)*/
	int base64plusheadlength = 24 + 4*ceil((float) compressedarraylength/3.0);
	unsigned char * base64plusheadarray= malloc(sizeof(unsigned char)* base64plusheadlength);

	/* fills base64plusheadarray with everything ready for simple writing*/
	base64plushead(compressedarray,compressedarraylength, chararraylength, base64plusheadarray);
	
	write_vtsarray(base64plusheadlength, base64plusheadarray, f);
	free(chararray);
	free(base64plusheadarray);
}

void write_pvts (struct All_variables *E, int cycles)
{
	FILE *fp;
	char pvts_file[255];
	int i,j,k;
	snprintf(pvts_file, 255, "%s.%d.pvts",
	E->control.data_file, cycles);
	fp = output_open(pvts_file, "w");

	const char format[] =
        	"<?xml version=\"1.0\"?>\n"
        	"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" compressor=\"vtkZLibDataCompressor\" byte_order=\"LittleEndian\">\n"
        	"  <PStructuredGrid WholeExtent=\"%s\" GhostLevel=\"#\">\n"
		"    <PPointData Scalars=\"temperature\" Vectors=\"velocity\">\n"
		"      <DataArray type=\"Float32\" Name=\"temperature\" format=\"binary\"/>\n"
		"      <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"binary\"/>\n"
		"      <DataArray type=\"Float32\" Name=\"viscosity\" format=\"binary\"/>\n";
    	
	char extent[64], header[1024];

    	snprintf(extent, 64, "%d %d %d %d %d %d",
             	E->lmesh.ezs, E->lmesh.ezs + E->lmesh.elz*E->parallel.nprocz,
            	E->lmesh.exs, E->lmesh.exs + E->lmesh.elx*E->parallel.nprocx,
             	E->lmesh.eys, E->lmesh.eys + E->lmesh.ely*E->parallel.nprocy);

    	snprintf(header, 1024, format, extent);

    	fputs(header, fp);
	
	if (E->output.stress){
		fputs("      <DataArray type=\"Float32\" Name=\"stress\" NumberOfComponents=\"6\" format=\"binary\"/>\n", fp);}
	if (E->output.comp_nd && E->composition.on){
		fputs("      <DataArray type=\"Float32\" Name=\"composition1\" format=\"binary\"/>\n", fp);}

	fputs(	"    </PPointData>\n \n"
		"    <PCellData>\n"
		"    </PCellData>\n \n"
		"    <PPoints>\n"
		"      <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"binary\" />\n"
		"    </PPoints>\n", fp);

	for(i=0; i < E->parallel.nprocz;i++){
		for(j=0; j < E->parallel.nprocx;j++){
			for(k=0; k < E->parallel.nprocy;k++){

		fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s.proc%d.%d.vts\"/>\n",
			(k%E->parallel.nprocz)*E->lmesh.elz, (k%E->parallel.nprocz + 1)*E->lmesh.elz, 
			(j%E->parallel.nprocx)*E->lmesh.elx, (j%E->parallel.nprocx+1)*E->lmesh.elx, 
			(i%E->parallel.nprocy)*E->lmesh.ely, (i%E->parallel.nprocy+1)*E->lmesh.ely,
			E->control.data_prefix, 
			i*E->parallel.nprocx*E->parallel.nprocy+j*E->parallel.nprocy+k, cycles);

			}
		}
	}

	fputs(	"  </PStructuredGrid>\n",fp);
	fputs(	"</VTKFile>",fp);

	fclose(fp);
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

    if (E->parallel.me == 0) {write_pvts(E, cycles);}


    return;
}
