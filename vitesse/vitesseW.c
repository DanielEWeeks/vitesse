

#ifdef INFO_HEADER
/*
* Program:       VITESSE 
* Function:      Likelihood Calculation on Pedigrees
* Version:       1.0
* Date:          1995.11.29
* Author:        Jeff O'Connell 
* Copyright:     (c) 1995 Jeff O'Connell 
* Collaborators: Daniel E. Weeks
* Language:      C (ANSI standard)
* Distribution:  via anonymous FTP from watson.hgen.pitt.edu/pub/vitesse.tar.Z
* Registration:  via e-mail to jeff@sherlock.hgen.pitt.edu
* Documentation: README file 
*/
#endif
#include "v_prog.h"

/*
#define INCL 
*/


/*********************************************************/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/************************************************************/
/* struct no longer has these fields  
void free_v_loci(ncU_kBgzMo * A[],int LT_URkV)
{
	int    k;

	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->child_patall[k][LT_URkV]);
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->child_matall[k][LT_URkV]);
	}
	  fprintf(OUTFILE, " <-> ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->dad_recomb[k][LT_URkV]);
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->mom_recomb[k][LT_URkV]);
	}
	fprintf(OUTFILE, "  #  ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->child_position[k][LT_URkV]);
	}
	  
	fprintf(OUTFILE, "\n");
}
*/



void generate_founder(GLIST *A[],int LT_URkV)
{
	int    k;

	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->aoNzgDs);
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->gmBiiBb_JmwFc_Qgi);
	}
	  fprintf(OUTFILE, " <-> ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->noBhhFh);
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->mnQgi);
	}
	  fprintf(OUTFILE, "  #  ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->aoFov_uiFj);
	}
	fprintf(OUTFILE, "\n");
}




void generate_genotypes(GLIST2 *A[],int LT_URkV)
{
	int    k;

	for (k = 0; k < hnPabHlgF; k++)
	{
	  /*
	  fprintf(OUTFILE, "%d ", A[k]->aoNzgDs);
	  */
	  fprintf(OUTFILE, "  ");
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  /*
	  fprintf(OUTFILE, "%d ", A[k]->gmBiiBb_JmwFc_Qgi);
	  */
	  fprintf(OUTFILE, "  ");
	}
	  fprintf(OUTFILE, " <-> ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->noBhhFh);
	}
	  fprintf(OUTFILE, "|");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->mnQgi);
	}
	  fprintf(OUTFILE, "  #  ");
	for (k = 0; k < hnPabHlgF; k++)
	{
	  fprintf(OUTFILE, "%d ", A[k]->aoFov_uiFj);
	}
	fprintf(OUTFILE, "\n");
}



/*********************************************************/
void delete_allele_up(ncU_kBgzMo * A[]) 
{
}
