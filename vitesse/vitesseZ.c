

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


/*
#define ALL_OUT
#define SYMMETRIC
#define FINAL_OUT
#define GRAPH
*/


void set_per_pentrance2(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,int *ce_zwK,PERSON *ao1)
{
   int k,plCzmE_yBhv;
   int  pi_xvOg;
   double hkMl_GivR;
   v_boolean dmF_oPlk,ow_zoMvoF;
   double nn_ezMfvT;
   int siU_rOwrDvh_kgS;    

do_lik = TRUE;

if(do_lik == TRUE)
{
     for(k=0;k<gm_c_TklVhv;k++)
        mkGfm[k]=0.0L;

     for(k=0;k<hnPabHlgF;k++)
     {
       if(ciS_kCrgT[k] == TRUE)
       {
          siU_rOwrDvh_kgS = k;
          k = hnPabHlgF;
       }
     }

    
     /*
     plCzmE_yBhv=QuicksortFounderD(ceBi_EvgFin,ce_zwK);
     */
     FmBo_nfgQfg=QuicksortFounder(ceBi_EvgFin);

/* 
     if( ciS_kPhrUrlO > hnPabHlgF-1 || ciS_kPhrUrlO == -1)
*/
     if( fnJorFh == 0)
     {
	 
	if(free_genotypes(dnNb,ceBi_EvgFin,ao1) == TRUE)
	{
     /*
     plCzmE_yBhv=QuicksortFounderD(ceBi_EvgFin,ce_zwK);
     */
     plCzmE_yBhv=QuicksortDouble(ceBi_EvgFin,ce_zwK);
          checkresult(ceBi_EvgFin[plCzmE_yBhv]);
          ceBi_EvgFin[plCzmE_yBhv]=processparent(dnNb,ceBi_EvgFin[plCzmE_yBhv],&ce_zwK[plCzmE_yBhv]);
	  /*
          ceBi_EvgFin[plCzmE_yBhv]=processparent(ceBi_EvgFin[plCzmE_yBhv],ce_zwK);
	  */
        }
	 
        delete_first_bit_array(dnNb,ceBi_EvgFin,ao1); 
     }
     else
     {
	/* ciSvmU_wVzo_kiPyzOw lmHgs present */
        ow_zoMvoF = FALSE;
      free_genotypes(dnNb,ceBi_EvgFin,ao1);
      /*
      if(TRUE) 
      */
      if(qeBirBmxF == NULL && nd_liEvi == NULL)
      {
         delete_first_bit_array(dnNb,ceBi_EvgFin,ao1); 
      }
      else
      {
       if(qeBirBmxF != NULL)
       {
	 ow_zoMvoF = TRUE;
       /*
         fprintf(OUTFILE, "\n Symmetric Part  Present");
             fflush(OUTFILE);
       */
       }
        dmF_oPlk =FALSE;

if(screen_out == TRUE)
{
        fprintf(stdout,"\n Peeling Non_Symmetric_Pair");
}

	/*
        nd_liEvi = NULL;
	*/
        if(nd_liEvi != NULL)
        {
           ceBi_EvgFin[siU_rOwrDvh_kgS]=nd_liEvi;
       
       if(qeBirBmxF != NULL)
         nuclear_fam_elim(dnNb,qeBirBmxF);
      
          delete_first_bit_array(dnNb,ceBi_EvgFin,ao1); 
          dmF_oPlk =TRUE;
       }
       else
       {
         ceBi_EvgFin[siU_rOwrDvh_kgS]=qeBirBmxF;
if(screen_out)
{
         fprintf(stdout,"\n Peeling Symmetric Part");
}

          FmBo_nfgQfg=QuicksortFounder(ceBi_EvgFin);

     /*
     plCzmE_yBhv=QuicksortFounderD(ceBi_EvgFin,ce_zwK);
     */
     plCzmE_yBhv=QuicksortDouble(ceBi_EvgFin,ce_zwK);
          checkresult(ceBi_EvgFin[plCzmE_yBhv]);
	  /*
          ceBi_EvgFin[plCzmE_yBhv]=processparent(ceBi_EvgFin[plCzmE_yBhv],ce_zwK);
	  */
          ceBi_EvgFin[plCzmE_yBhv]=processparent(dnNb,ceBi_EvgFin[plCzmE_yBhv],&ce_zwK[plCzmE_yBhv]);
	  
	 
          delete_first_bit_array(dnNb,ceBi_EvgFin,ao1); 
          ow_zoMvoF = FALSE;
       }


      if(ow_zoMvoF == TRUE)
      {
        if(dmF_oPlk == TRUE)
        {
         for(k=0;k<hnPabHlgF;k++)
         {
           ceBi_EvgFin[k]=display_chidren(dnNb,k,ao1);
           if(ceBi_EvgFin[k] == NULL)
           {  
             fprintf(OUTFILE,"\n Top[%d] = NULL \n",k);
             exit(1);
           }
	 }
        }   
       free_genotypes(dnNb,ceBi_EvgFin,ao1);
       ceBi_EvgFin[siU_rOwrDvh_kgS]=qeBirBmxF;
       
     /*
     plCzmE_yBhv=QuicksortFounderD(ceBi_EvgFin,ce_zwK);
     */
     plCzmE_yBhv=QuicksortDouble(ceBi_EvgFin,ce_zwK);
          FmBo_nfgQfg=QuicksortFounder(ceBi_EvgFin);

       if(nd_liEvi != NULL)
         nuclear_fam_elim(dnNb,nd_liEvi);
       	
       checkresult(ceBi_EvgFin[plCzmE_yBhv]);
       /*
       ceBi_EvgFin[plCzmE_yBhv]=processparent(ceBi_EvgFin[plCzmE_yBhv],ce_zwK);
       */
       ceBi_EvgFin[plCzmE_yBhv]=processparent(dnNb,ceBi_EvgFin[plCzmE_yBhv],&ce_zwK[plCzmE_yBhv]);
       delete_first_bit_array(dnNb,ceBi_EvgFin,ao1); 
     }
    }
   } 
      qeBirBmxF = NULL;
      nd_liEvi = NULL;

     /*
     isozygote_combined(dnNb,ceBi_EvgFin);
     ceBi_EvgFin=pick_pairs(dnNb,ceBi_EvgFin);
     */
   /* 
   for(k=0;k<hnPabHlgF;k++)
   {
    freeGLIST(ceBi_EvgFin[k]);
    printloci(ceBi_EvgFin[k]);
    ceBi_EvgFin[k]=processparent(dnNb,ceBi_EvgFin[k],ce_zwK);
    freeGLIST(ceBi_EvgFin[k]);
   }
   */ 
}

}

/*##################  */


void set_per_pentrance(NUC_FAM *dnNb,ncU_kBgzMo **T,PERSON *ao1)
{
     int k;

     for(k=0;k<gm_c_TklVhv;k++)
        mkGfm[k]=0.0L;

        FmBo_nfgQfg=QuicksortFounder(T);
	 
        delete_first_bit_array(dnNb,T,ao1); 

}


void set_pos_index(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,int *ce_zwK,PERSON *ao1)
{
   int k,plCzmE_yBhv;
   int  pi_xvOg;
   double hkMl_GivR;
   v_boolean dmF_oPlk,ow_zoMvoF;
   double nn_ezMfvT;
    

     for(k=0;k<gm_c_TklVhv;k++)
        mkGfm[k]=0.0L;

    
     /*
     plCzmE_yBhv=QuicksortFounderD(ceBi_EvgFin,ce_zwK);
     */
     plCzmE_yBhv=QuicksortDouble(ceBi_EvgFin,ce_zwK);
     FmBo_nfgQfg=QuicksortFounder(ceBi_EvgFin);


}


void set_pos_index2(NUC_FAM *dnNb,ncU_kBgzMo **ceBi_EvgFin,PERSON *ao1)
{
   int k,j;
   long lt_rmGrmJgb;
   double *nnQvw;
    

    /*
    do_couple_list_inter2(dnNb,ceBi_EvgFin);
    */
    lt_rmGrmJgb=ao1->nn_fhFh;
   
     for(k=0;k<gm_c_TklVhv;k++)
        mkGfm[k]=0.0L;

   for(j=0;j<gm_c_TklVhv;j++)
   {
     nnQvw=ao1->tnQy[j];
     for(k=0;k<lt_rmGrmJgb;k++)
     {
       mkGfm[j]+=nnQvw[k];
     }
    }
     for(k=0;k<gm_c_TklVhv;k++)
        mkGfm[k]*=2.0L;



}


