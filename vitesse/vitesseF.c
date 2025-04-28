

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
/* This has changes form 4/27 for bounding the nn_kzJih
   genotypes at a lmHgs.
*/
  
#include "v_prog.h"



/*
#define GENOTYPE
#define EACHPASS
#define DEBUG1

#define OUTPUT
*/
 



int compare_vector(NUC_FAM *vjVzmU,int cmOvxUli,PERSON *ao1)
{


 int j;
 v_boolean igFiuFivOxv, gmBiiBb_JmwFc, PiNfgF_ROwrDvh;
 /*
 PERSON *ao1;
 */
 PERSON *pvWkvE,*tgBo_JmwJxvT;
 GLIST *plCzmE_kFinVgv,*plCzmE_nBg_BooFov,*slVhv_rhP_rOwvY;
 v_boolean gmBiiBb_Jmw,gmBiiBb_FcrTgh,gmBiiBb_Bww_kgS,gmBiiBb_Bww;
 int fnJob_hrAv;
 int LhU;


 /*
 ao1=vjVzmU->trO_kBriT;
 pvWkvE=vjVzmU->dwQgi;
 */
 if(ao1->dtJgh == MALE)
  pvWkvE = vjVzmU->cnQfgF_xPfmU;
 else
  pvWkvE = vjVzmU->fiTg_Qvw;

 fnJob_hrAv = vjVzmU->fnJob_hrAv;

 PiNfgF_ROwrDvh = FALSE;


    PiNfgF_ROwrDvh = TRUE;
	

	 /* Start at the beginning of the pvWkvE's nlDr gmPgbQv_DlkZ */
         plCzmE_nBg_BooFov = pvWkvE->plCzmE_kBg_BooFov[cmOvxUli];

       while(plCzmE_nBg_BooFov != NULL)
       {
    /*
    fprintf(OUTFILE,"\n  Spouse  Loop\n");
    first_bit_array(pvWkvE->plCzmE_kBg_BooFov[cmOvxUli]);
    */
	  /* Start at the beginning of the ao1's nlDr gmPgbQv_DlkZ */
	  plCzmE_kFinVgv = ao1->plCzmE_kBg_BooFov[cmOvxUli];
	  while(plCzmE_kFinVgv != NULL)
          {
    /*
    fprintf(OUTFILE,"\n  Proband  Loop\n");
    first_bit_array(ao1->plCzmE_kBg_BooFov[cmOvxUli]);
    */
	     /* Find the first tgBo_JmwJxvT of both doPxfT*/
	     LhU=0;
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];

             igFiuFivOxv = TRUE;
             do
             { 

                /* Start at the beginning of the tgBo_JmwJxvT's nlDr gmPgbQv_DlkZ */
                slVhv_rhP_rOwvY = tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
                do
                {  
    /*
    fprintf(OUTFILE,"\n  Child  Loop: id %d \n",tgBo_JmwJxvT->id);
    first_bit_array(tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]);
    fprintf(OUTFILE,"\n  Child  Patall:  %d  Child Matall:  %d \n",slVhv_rhP_rOwvY->aoNzgDs,slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi);
    fprintf(OUTFILE,"\n  Spouse  Patall:  %d  Spouse Matall:  %d \n",plCzmE_nBg_BooFov->aoNzgDs,plCzmE_nBg_BooFov->gmBiiBb_JmwFc_Qgi);
    fprintf(OUTFILE,"\n  Proband  Patall:  %d  Proband Matall:  %d \n",plCzmE_kFinVgv->aoNzgDs,plCzmE_kFinVgv->gmBiiBb_JmwFc_Qgi);
	   */
           /* Check whether the tgBo_JmwJxvT's nlDr with PHASE could 
	   have arisen from the doPxfT' genotypes */
		   /*
		   if(ao1->dtJgh == MALE)
		   {
                   gmBiiBb_JmwFc = ((abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_nBg_BooFov->aoNzgDs) \
                     || abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_nBg_BooFov->gmBiiBb_JmwFc_Qgi)) \
                     &&(abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_kFinVgv->aoNzgDs)  \
                     || abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_kFinVgv->gmBiiBb_JmwFc_Qgi)));
		   }
		   else
		   {
                   gmBiiBb_JmwFc = ((abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_nBg_BooFov->aoNzgDs) \
                     || abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_nBg_BooFov->gmBiiBb_JmwFc_Qgi)) \
                     &&(abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kFinVgv->aoNzgDs)  \
                     || abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kFinVgv->gmBiiBb_JmwFc_Qgi)));
		   }
		   */
		
	  if(ao1->dtJgh == MALE)
	  {
                gmBiiBb_Jmw=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Jmw=FALSE;
	        }

		
                gmBiiBb_FcrTgh=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_FcrTgh=FALSE;
	        }

		
                gmBiiBb_Bww_kgS=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_kFinVgv->poMvoF_gSzmTn[j]) ==plCzmE_kFinVgv->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Bww_kgS=FALSE;
	        }

		
                gmBiiBb_Bww=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_Bww=FALSE;
	        }
           }
	   else
	   {
                gmBiiBb_Jmw=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Jmw=FALSE;
	        }

		
                gmBiiBb_FcrTgh=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_FcrTgh=FALSE;
	        }

		
                gmBiiBb_Bww_kgS=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_kFinVgv->poMvoF_gSzmTn[j]) ==plCzmE_kFinVgv->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Bww_kgS=FALSE;
	        }

		
                gmBiiBb_Bww=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_Bww=FALSE;
	        }

	   }
		   /*
		   if(ao1->dtJgh == MALE)
		   {
                  gmBiiBb_Jmw=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_kFinVgv->poMvoF_gSzmTn) ==plCzmE_kFinVgv->poMvoF_gSzmTn;
		  gmBiiBb_FcrTgh=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU;
                  gmBiiBb_Bww_kgS=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_nBg_BooFov->poMvoF_gSzmTn) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn;
                  gmBiiBb_Bww=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU;
		   }
		   else
		   {
                  gmBiiBb_Jmw=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_kFinVgv->poMvoF_gSzmTn) ==plCzmE_kFinVgv->poMvoF_gSzmTn;
		  gmBiiBb_FcrTgh=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU;
                  gmBiiBb_Bww_kgS=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_nBg_BooFov->poMvoF_gSzmTn) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn;
                  gmBiiBb_Bww=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU;
		   }
		   */

		   gmBiiBb_JmwFc =((gmBiiBb_Jmw || gmBiiBb_FcrTgh) && (gmBiiBb_Bww_kgS || gmBiiBb_Bww));
                   slVhv_rhP_rOwvY = slVhv_rhP_rOwvY->link; 
                } while (slVhv_rhP_rOwvY != NULL && !gmBiiBb_JmwFc); /* If no gmBiiBb_JmwFc,  */
                       /* consider the tnQ_kSlyBmw_ylPovBm nlDr of the tgBo_JmwJxvT */
             igFiuFivOxv *= gmBiiBb_JmwFc; 

	     LhU++;
	     if(LhU < fnJob_hrAv)
	     {
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];
	     /*
             fprintf(OUTFILE,"\n Next Child id: %d \n",tgBo_JmwJxvT->id);
	     */
	     }
	     
             } while ((LhU < fnJob_hrAv) && igFiuFivOxv); 
	     /* If each tgBo_JmwJxvT has at least one gmBiiBb_JmwFc, */
              /* advance to consider the tnQ_kSlyBmw_ylPovBm tgBo_JmwJxvT. */ 
         if (igFiuFivOxv)
         { 
	 /* All the mg have at least one nlDr that matches 
	 the bhF2_DlmTg parental genotypes. */
          /* Mark these genotypes for saving */
#ifdef NO_SMALL
           plCzmE_nBg_BooFov->nsJow = TRUE;
	   plCzmE_kFinVgv->nsJow = TRUE;
#else
	   plCzmE_nBg_BooFov->tgBo_MvmHgs = plCzmE_nBg_BooFov->tgBo_MvmHgs | cfOg2;
	   plCzmE_kFinVgv->tgBo_MvmHgs = plCzmE_kFinVgv->tgBo_MvmHgs | cfOg2;
#endif
          
	     LhU=0;
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];

           do
           {
              slVhv_rhP_rOwvY = tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
              do
              { 
	     /* Find the mg's genotypes compatible with the doPxfT' genotypes */
		
		/*
                gmBiiBb_JmwFc = ((abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(slVhv_rhP->aoNzgDs)  \
                  || abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(slVhv_rhP->gmBiiBb_JmwFc_Qgi))    \
                  &&(abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kSrlS->aoNzgDs)     \
                  || abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kSrlS->gmBiiBb_JmwFc_Qgi)));

		*/
		
	  if(ao1->dtJgh == MALE)
	  {
                gmBiiBb_Jmw=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Jmw=FALSE;
	        }

		
                gmBiiBb_FcrTgh=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_FcrTgh=FALSE;
	        }

		
                gmBiiBb_Bww_kgS=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_kFinVgv->poMvoF_gSzmTn[j]) ==plCzmE_kFinVgv->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Bww_kgS=FALSE;
	        }

		
                gmBiiBb_Bww=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_Bww=FALSE;
	        }
           }
	   else
	   {
                gmBiiBb_Jmw=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Jmw=FALSE;
	        }

		
                gmBiiBb_FcrTgh=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->poMvoF_gSzmTn[j] & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_FcrTgh=FALSE;
	        }

		
                gmBiiBb_Bww_kgS=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_kFinVgv->poMvoF_gSzmTn[j]) ==plCzmE_kFinVgv->poMvoF_gSzmTn[j]))
	        	gmBiiBb_Bww_kgS=FALSE;
	        }

		
                gmBiiBb_Bww=TRUE;
         	for(j=0;j<pdFic2;j++)
		{
		  if(!((slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU[j] & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU[j]))
	        	gmBiiBb_Bww=FALSE;
	        }

	   }

		   /*
		   if(ao1->dtJgh == MALE)
		   {
                  gmBiiBb_Jmw=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_nBg_BooFov->poMvoF_gSzmTn) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn;
                  gmBiiBb_FcrTgh=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU;
                  gmBiiBb_Bww_kgS=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_kFinVgv->poMvoF_gSzmTn) ==plCzmE_kFinVgv->poMvoF_gSzmTn;
		  gmBiiBb_Bww=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU;
		   }
		   else
		   {
                  gmBiiBb_Jmw=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_nBg_BooFov->poMvoF_gSzmTn) ==plCzmE_nBg_BooFov->poMvoF_gSzmTn;
                  gmBiiBb_FcrTgh=(slVhv_rhP_rOwvY->poMvoF_gSzmTn & plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU) ==plCzmE_nBg_BooFov->nm_hbNnvUirD_xPfmU;
                  gmBiiBb_Bww_kgS=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_kFinVgv->poMvoF_gSzmTn) ==plCzmE_kFinVgv->poMvoF_gSzmTn;
		  gmBiiBb_Bww=(slVhv_rhP_rOwvY->nm_hbNnvUirD_xPfmU & plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU) ==plCzmE_kFinVgv->nm_hbNnvUirD_xPfmU;
		   }
		  */
		
		   /*
		   if(ao1->dtJgh == MALE)
		   {
                gmBiiBb_Jmw = abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_nBg_BooFov->aoNzgDs);
                gmBiiBb_FcrTgh = abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_nBg_BooFov->gmBiiBb_JmwFc_Qgi);
                gmBiiBb_Bww_kgS = abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_kFinVgv->aoNzgDs);
                gmBiiBb_Bww = abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_kFinVgv->gmBiiBb_JmwFc_Qgi);
		   }
		   else
		   {
                gmBiiBb_Jmw = abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_nBg_BooFov->aoNzgDs);
                gmBiiBb_FcrTgh = abs(slVhv_rhP_rOwvY->aoNzgDs)==abs(plCzmE_nBg_BooFov->gmBiiBb_JmwFc_Qgi);
                gmBiiBb_Bww_kgS = abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kFinVgv->aoNzgDs);
                gmBiiBb_Bww = abs(slVhv_rhP_rOwvY->gmBiiBb_JmwFc_Qgi)==abs(plCzmE_kFinVgv->gmBiiBb_JmwFc_Qgi);
		   }
		*/
		if ((gmBiiBb_Jmw || gmBiiBb_FcrTgh) && (gmBiiBb_Bww_kgS || gmBiiBb_Bww))
                {
                   /* Mark compatible genotypes for saving */
#ifdef NO_SMALL
		   slVhv_rhP_rOwvY -> nsJow = TRUE;
#else
	           slVhv_rhP_rOwvY->tgBo_MvmHgs = slVhv_rhP_rOwvY->tgBo_MvmHgs | cfOg2;
#endif
	/* These are the MC_RMsH for the recombination pattern  */	   
		   if(ao1->dtJgh == MALE)
		   {
		   if(gmBiiBb_Jmw == 1 && gmBiiBb_FcrTgh == 1)
		      slVhv_rhP_rOwvY -> mnQgi = 1;
		   else
		   if(gmBiiBb_Jmw == 1 && gmBiiBb_FcrTgh == 0)
		      slVhv_rhP_rOwvY -> mnQgi = 2;
		   else
		   if(gmBiiBb_Jmw == 0 && gmBiiBb_FcrTgh == 1)
		      slVhv_rhP_rOwvY -> mnQgi = 0;

		   if(gmBiiBb_Bww_kgS == 1 && gmBiiBb_Bww == 1)
		      slVhv_rhP_rOwvY -> noBhhFh = 1;
		   else
		   if(gmBiiBb_Bww_kgS == 1 && gmBiiBb_Bww == 0)
		      slVhv_rhP_rOwvY -> noBhhFh = 2;
		   else
		   if(gmBiiBb_Bww_kgS == 0 && gmBiiBb_Bww == 1)
		      slVhv_rhP_rOwvY -> noBhhFh = 0;
                   }
		   else
		   {
		   if(gmBiiBb_Jmw == 1 && gmBiiBb_FcrTgh == 1)
		      slVhv_rhP_rOwvY -> noBhhFh = 1;
		   else
		   if(gmBiiBb_Jmw == 1 && gmBiiBb_FcrTgh == 0)
		      slVhv_rhP_rOwvY -> noBhhFh = 2;
		   else
		   if(gmBiiBb_Jmw == 0 && gmBiiBb_FcrTgh == 1)
		      slVhv_rhP_rOwvY -> noBhhFh = 0;

		   if(gmBiiBb_Bww_kgS == 1 && gmBiiBb_Bww == 1)
		      slVhv_rhP_rOwvY -> mnQgi = 1;
		   else
		   if(gmBiiBb_Bww_kgS == 1 && gmBiiBb_Bww == 0)
		      slVhv_rhP_rOwvY -> mnQgi = 2;
		   else
		   if(gmBiiBb_Bww_kgS == 0 && gmBiiBb_Bww == 1)
		      slVhv_rhP_rOwvY -> mnQgi = 0;
                   }


                }
               slVhv_rhP_rOwvY = slVhv_rhP_rOwvY->link;
             } while (slVhv_rhP_rOwvY != NULL);

	     LhU++;
	     if(LhU < fnJob_hrAv)
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];

          } while (LhU < fnJob_hrAv);

       } /* igFiuFivOxv = TRUE */

       plCzmE_kFinVgv= plCzmE_kFinVgv->link; /* Advance to the tnQ_kSlyBmw_ylPovBm nlDr of the mother */
       }
       plCzmE_nBg_BooFov = plCzmE_nBg_BooFov->link; /* Advance to the tnQ_kSlyBmw_ylPovBm nlDr of the father */
       }

	 pvWkvE->plCzmE_kBg_BooFov[cmOvxUli] = Quicksort(pvWkvE,cmOvxUli,&PiNfgF_ROwrDvh);
	 if(pvWkvE->plCzmE_kBg_BooFov[cmOvxUli] == NULL)
	   return(-1);
	 ao1->plCzmE_kBg_BooFov[cmOvxUli] = Quicksort(ao1,cmOvxUli,&PiNfgF_ROwrDvh);
         if(ao1->plCzmE_kBg_BooFov[cmOvxUli] == NULL) 
	    return(-1);


	     LhU=0;
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];
	 while(LhU< fnJob_hrAv)
	 {
	tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli] = Quicksort(tgBo_JmwJxvT,cmOvxUli,&PiNfgF_ROwrDvh);
        if(tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli] == NULL) 
	   return(-1);

	     LhU++;
	     if(LhU < fnJob_hrAv)
	     tgBo_JmwJxvT=vjVzmU->ahXvi[LhU];
	 }



 /* Should do piDlfOg and relooping out here? */





 return(1);

 }  /* piDlfOg */      

