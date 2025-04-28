

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
/*  this has the extra statements for the
   debugging  */
/* 
 This version uses chars instead of longs for 
 the sort.
*/
  
#include "v_prog.h"

/*
#define CHILD_LEN
#define ORDERLIST 
*/


/**********************************************/

/* This handles the code for the founders  */


void  mark_genotypes(NUC_FAM  *pi_hvY,ncU_kBgzMo *T[],MLIST **P )
{


    int fnJob_hrAv,j,k,slVhv_klT;
    int  **gmPgbQv;

    ncU_kBgzMo  **junk3;
    char  **Il_xoBhh_evDglS;
    char *siU_lVgkVg_JmwFc,*co_m;
    ncU_kBgzMo  *F_ptr;
    int nn_xlMfnOh;
    MLIST  *nd_nzMv_Svx;

   fnJob_hrAv = pi_hvY->fnJob_hrAv;

   gmPgbQv=(int**)v_alloc(fnJob_hrAv ,sizeof(int*));

   for(j=0;j<fnJob_hrAv;j++)
   {

   gmPgbQv[j]=(int*)v_alloc(hnPabHlgF,sizeof(int));

     gmPgbQv[j][hnPabHlgF-1]=1;
     for(k=hnPabHlgF-2;k>=0;k--)
     {
       /*
       gmPgbQv[j][k]=gmPgbQv[j][k+1]*(T[k+1]->hg_nfMgrQorFi[j]->Il_GiBmh_RgFi);
       */
       gmPgbQv[j][k]=gmPgbQv[j][k+1]*(T[k+1]->tg_ulVmwFi_Dmg[j]);
     }
   }


     junk3 = (ncU_kBgzMo**)v_alloc(FmBo_nfgQfg+OFFSET, sizeof(ncU_kBgzMo*));

     Il_xoBhh_evDglS = (char**)v_alloc(FmBo_nfgQfg+OFFSET, sizeof(char*));

     for(k=0;k<=FmBo_nfgQfg;k++)
     {

       Il_xoBhh_evDglS[k]=(char*)v_alloc(fnJob_hrAv,sizeof(char));
     }

 /* added to handle like calculations  */

     /*
     P = (MLIST **)v_alloc(hnPabHlgF, sizeof(MLIST*));
     */

     /*
     fprintf(OUTFILE, "\n Sorting Glists  ");
     */
     for (k = 0; k < hnPabHlgF; k++)
     {

      /*
      fprintf(OUTFILE, "\n\n Locus  %d ",k);
      */
      F_ptr=T[k];
      slVhv_klT=0;
      while(F_ptr != NULL)
      {
	junk3[slVhv_klT]=F_ptr;
	/*
	convergence(F_ptr,k);
	*/
	for(j = 0; j < fnJob_hrAv; j++)
	{
	    Partition3(F_ptr->iwFc[j],F_ptr->hg_nfMgrQorFi[j],k);
	    PartitionFound(F_ptr->iwFc[j],F_ptr->hg_nfMgrQorFi[j],gmPgbQv[j][k]);
	    Il_xoBhh_evDglS[slVhv_klT][j]=F_ptr->lt10_UlgBo_Tfn[j];
	}
	F_ptr=F_ptr->link;
        slVhv_klT++;
      }

     /* NCH */
     Add_Alist(Il_xoBhh_evDglS,junk3,0,slVhv_klT-1,fnJob_hrAv);



      P[k] = (MLIST *)v_alloc(1,sizeof(MLIST));
     
        nd_nzMv_Svx=P[k];
        P[k]->founder_pair=junk3[0];
	nn_xlMfnOh=1;
	nd_nzMv_Svx->plCzmE_yBhv=nn_xlMfnOh;
        co_m=Il_xoBhh_evDglS[0];
	
	for(j=1;j<slVhv_klT;j++)
	{
          siU_lVgkVg_JmwFc=Il_xoBhh_evDglS[j];
          if(!assign_rel(siU_lVgkVg_JmwFc,co_m,fnJob_hrAv))
          {
	      nd_nzMv_Svx->link = (MLIST *)v_alloc(1,sizeof(MLIST));

                nd_nzMv_Svx = nd_nzMv_Svx->link;
	        nd_nzMv_Svx->plCzmE_yBhv=++nn_xlMfnOh;
		nd_nzMv_Svx->founder_pair=junk3[j];
                junk3[j-1]->link=NULL;
          }
	  else
	     junk3[j-1]->link=junk3[j];

           co_m=siU_lVgkVg_JmwFc;

        }
	 junk3[slVhv_klT-1]->link=NULL;
      

   }


   /* free gm_hlSg_JmwJxvT  */

   for(j=0;j<fnJob_hrAv;j++)
   {
       free_check(gmPgbQv[j]);
   }

   free_check(gmPgbQv);
    
   free_check(junk3);


     for(k=0;k<=FmBo_nfgQfg;k++)
     {
      free_check(Il_xoBhh_evDglS[k]);
     }
      free_check(Il_xoBhh_evDglS);


/*   end of code for new combinations  */

}
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

