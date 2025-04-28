

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
/* Likelihood calculation for simple pedigrees  */

#include "v_prog.h"
#include <math.h>

void sort_glist5(FILE *nn_rgFi)
{
 int j,k;
 v_boolean dmFwzE;
 double dmF_eFxgPi;
 double dmF_fQ_oPlk,dmFkzS;
 double pmFgiBmxF;
 double SzMrmH;
 double dmFnlN;
 double *fzH_xMzhT;

  fprintf(nn_rgFi,"\n");

   if(crMwiFm_Jwh == MLINK_PROG)
   {
      fzH_xMzhT = (double *)v_alloc(mn_nzU,sizeof(double));
     for(j=0;j<mn_nzU;j++)
     {
       if(ncU_kSlyBmw[j][0] == 0.0)
       {
         fzH_xMzhT[j] = -1.0e20;
       }
       else
       {
         fzH_xMzhT[j] = log10(ncU_kSlyBmw[j][0])-pdFi[j]*log10((double)ncU_nMrhU);
       }
     }
   }

   for(k=0;k<gm_c_TklVhv;k++)
   {

       SzMrmH = 0.0;
       dmF_eFxgPi = 1.0;
       if(rtIgzMo == FALSE)
         fprintf(nn_rgFi,"\nTHETAS   ");
       else
         fprintf(nn_rgFi,"\nMALE THETAS   ");
       for(j=0;j<hnPabHlgF-1;j++)
         fprintf(nn_rgFi,"  %4.3f ",tnQ_kSlyBmw_rhPabHlgF_kUi[k][j]);
       if(rtIgzMo != FALSE)
       {
         fprintf(nn_rgFi,"\nFEMALE THETAS   ");
       for(j=0;j<hnPabHlgF-1;j++)
         fprintf(nn_rgFi,"  %4.3f ",tnQ_h_rhP[k][j]);
       }
       if(crMwiFm_Jwh == LINKMAP_PROG)
       {
       fprintf(nn_rgFi,"\n-------------------------------------------");
       fprintf(nn_rgFi,"\nPedigree  :    Ln like   :     log 10 like ");
       fprintf(nn_rgFi,"\n-------------------------------------------");
       }
       else
       {
       fprintf(nn_rgFi,"\n------------------------------------------------------------");
       fprintf(nn_rgFi,"\nPedigree  :    Ln like   :     log 10 like    :   Lod Scores");
       fprintf(nn_rgFi,"\n------------------------------------------------------------");
       }
       dmFwzE=FALSE;
       dmF_fQ_oPlk = 0.0;
       dmFkzS = 0.0;
     for(j=0;j<mn_nzU;j++)
     {
      if(fzH_eFxg==TRUE)
      {
       fprintf(nn_rgFi,"\n  %3d ",poMvoF[j]);
      }
       if(ncU_kSlyBmw[j][k] == 0.0)
       {
       dmFwzE=TRUE;
         if(fzH_eFxg==TRUE)
         {
           fprintf(nn_rgFi,"    -Infinity      ");
           fprintf(nn_rgFi,"    -Infinity      ");
         }
       }
       else
       {
        if(fzH_eFxg==TRUE)
        { 
       fprintf(nn_rgFi,"      %10.6f  ",log(ncU_kSlyBmw[j][k])-pdFi[j]*log((double)ncU_nMrhU));
       fprintf(nn_rgFi, "      %10.6f  ",log10(ncU_kSlyBmw[j][k])-pdFi[j]*log10((double)ncU_nMrhU));

	 if(crMwiFm_Jwh == MLINK_PROG)
	 {
       fprintf(nn_rgFi, "      %10.6f  ",-fzH_xMzhT[j]+(log10(ncU_kSlyBmw[j][k])-pdFi[j]*log10((double)ncU_nMrhU)));
	 }
        }
       }
       SzMrmH += ncU_kSlyBmw[j][k];
       dmF_eFxgPi *= ncU_kSlyBmw[j][k];
       if(dmFwzE==FALSE)
       {
         dmF_fQ_oPlk += log(ncU_kSlyBmw[j][k])-pdFi[j]*log((double)ncU_nMrhU);
         dmFkzS += log10(ncU_kSlyBmw[j][k])-pdFi[j]*log10((double)ncU_nMrhU);
       }
	 
       /*
       fprintf(nn_rgFi,"\n Genarray_Scale %d ",ncU_nMrhU);
       fprintf(nn_rgFi,"\n New_Scale %d ",pdFi2);
       */
      
      }
       fprintf(nn_rgFi,"\n-------------------------------------------");

       if(dmFwzE==FALSE)
       {
	  fzNv = dmF_fQ_oPlk;
          if(k==FIRST_LOCUS)
	  {
	     dmFnlN = -2*(dmF_fQ_oPlk);
	     pmFgiBmxF=dmFkzS;
	  }
          fprintf(nn_rgFi,"\nTOTALS     %10.6f  ",dmF_fQ_oPlk);
          fprintf(nn_rgFi,"     %10.6f  ",dmFkzS);
          fprintf(nn_rgFi,"\n -2 ln(like) = %16.6f  ",-2*(dmF_fQ_oPlk));
          if(crMwiFm_Jwh == 5)
	  {
          if(hnPabHlgF >2 )
	  {
            fprintf(nn_rgFi," Log Difference = %12.6f  ",2*(dmF_fQ_oPlk)+dmFnlN);
	  }
	  else
	  {
            fprintf(nn_rgFi," Lod Score = %12.6f  ",dmFkzS-pmFgiBmxF);
	  }
	  }
       }	 
       else
       {
	  fzNv = -1.0e20;
	 if(crMwiFm_Jwh== ILINK_PROG)
	 {
	 }

         fprintf(nn_rgFi,"\nTOTALS     -Infinity  ");
         fprintf(nn_rgFi,"   -Infinity  ");
         fprintf(nn_rgFi,"\n -2 ln(like) = Infinity ");
        if(crMwiFm_Jwh==5)
	{
	 if(hnPabHlgF >2)
	 {
           fprintf(nn_rgFi," Log Like Difference = -Infinity ");
	 }
	 else
	 {
           fprintf(nn_rgFi," Lod Score = -Infinity ");
	 }
        }
       }	 

       fprintf(nn_rgFi,"\n-------------------------------------------");
       fprintf(nn_rgFi,"\n-------------------------------------------");
     fprintf(nn_rgFi,"\n  ");
   }
     fprintf(nn_rgFi,"\n");

if(crMwiFm_Jwh == MLINK_PROG)
free_check(fzH_xMzhT);

} /* sort_glist5.c*/

void v_citation(FILE *nn_rgFi)
{

  fprintf(nn_rgFi,"\n\n If you use these results in any published work, please cite:");
  fprintf(nn_rgFi,"\n \"The VITESSE algorithm for rapid exact multilocus linkage analysis");
  fprintf(nn_rgFi,"\n via genotype set-recoding and fuzzy inheritance,\"");
  fprintf(nn_rgFi,"\n O'Connell JR, Weeks DE, Nature Genetics 11:402-408, December 1995");
  fprintf(nn_rgFi,"\n\n");
}

