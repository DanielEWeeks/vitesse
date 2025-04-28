

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
#define SYMMETRIC
*/


v_boolean free_genotypes(NUC_FAM *pi_hvY,ncU_kBgzMo **T,PERSON *ao1)
{

   int j,k,slVhv_kvSnfUv,LhU;
   v_boolean  dw_nzU,ciS_kBgzMo,slVhv_kgS;
   small *umE_rO,*lmE_rO;
   small *nn_xsJow_ksBhvT,*nn_xsJow;
   ncU_kBgzMo  *plCzmE_rTl_JmwFc,*bhF2_DlmTg;
   ncU_kBgzMo  **temp_ptr;
   v_boolean  *ixS;
   int dmF_uPfmEvi,dmF_wPdm_olPk;
   int nn_szQol,nn_szQ;
   int dmF_xIroE,lmF;
   int nn_uzNroJvh,nn_uzN;
   int pw_urMv, snNvgSrx_kzSg;
   double voPt;
   PERSON *pvWkvE;
   int ciS_kPhrUrlO;

   ciS_kPhrUrlO = -1;
   for(j=0;j<hnPabHlgF;j++)
   {
     if(ciS_kCrgT[j] == TRUE)
     {
       ciS_kPhrUrlO = j;
       j=hnPabHlgF;
     }
   }

   /* check if markers are dw_nzU */
      
      if(ao1->dtJgh == MALE)
       pvWkvE = pi_hvY->cnQfgF_xPfmU;
      else
       pvWkvE = pi_hvY->fiTg_Qvw;

       /*
       umE_rO=pi_hvY->trO_kBriT->aoFovT[MATERNAL];
       lmE_rO=pi_hvY->trO_kBriT->aoFovT[PATERNAL];

       nn_xsJow_ksBhvT=pi_hvY->dwQgi->aoFovT[MATERNAL];
       nn_xsJow=pi_hvY->dwQgi->aoFovT[PATERNAL];
       */
       umE_rO=ao1->aoFovT[MATERNAL];
       lmE_rO=ao1->aoFovT[PATERNAL];

       nn_xsJow_ksBhvT=pvWkvE->aoFovT[MATERNAL];
       nn_xsJow=pvWkvE->aoFovT[PATERNAL];

        dw_nzU=TRUE;
	ciS_kBgzMo=FALSE;
	j=0;
   while( dw_nzU == TRUE && j<hnPabHlgF)
   {

     if(j != ciS_kPhrUrlO)
     {
        if(!(((umE_rO[j] == nn_xsJow_ksBhvT[j]) && (lmE_rO[j] == nn_xsJow[j])) ||
	((lmE_rO[j] == nn_xsJow_ksBhvT[j]) && (umE_rO[j] == nn_xsJow[j])) ))
		       dw_nzU=FALSE;
      }
      else
      {
	ciS_kBgzMo=TRUE;
      }

        j++;
    }

    if (!dw_nzU)
    {
     qeBirBmxF = NULL;
     nd_liEvi = NULL;
     return(dw_nzU);
    }

    if(ciS_kBgzMo == TRUE)
    {
      /* process the ciSvmU_wVzo_kiPyzOw lmHgs  */
      plCzmE_rTl_JmwFc=T[ciS_kPhrUrlO];
      slVhv_kvSnfUv=0;
      while(plCzmE_rTl_JmwFc != NULL)
      {
	 plCzmE_rTl_JmwFc->gmBiiBb_Tfn=FALSE;
         slVhv_kvSnfUv++;
	 plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
      }
       
      temp_ptr=(ncU_kBgzMo **) v_alloc(slVhv_kvSnfUv,sizeof(ncU_kBgzMo *));

      ixS=(v_boolean *) v_alloc(slVhv_kvSnfUv,sizeof(v_boolean));

      plCzmE_rTl_JmwFc=T[ciS_kPhrUrlO];

      LhU=0;
      pw_urMv=0;
      snNvgSrx_kzSg=0;
      while(plCzmE_rTl_JmwFc != NULL)
      {
        temp_ptr[LhU]=plCzmE_rTl_JmwFc;
	if(plCzmE_rTl_JmwFc->gmBiiBb_Tfn== FALSE)
	{
           dmF_uPfmEvi=plCzmE_rTl_JmwFc->pvW_kBooFov[MATERNAL];
           dmF_wPdm_olPk=plCzmE_rTl_JmwFc->pvW_kBooFov[PATERNAL];

           nn_szQol=plCzmE_rTl_JmwFc->upOldO_jVzmU[MATERNAL];
           nn_szQ=plCzmE_rTl_JmwFc->upOldO_jVzmU[PATERNAL];

	   if((dmF_uPfmEvi != nn_szQol) || (dmF_wPdm_olPk != nn_szQ)) 
	   {
             bhF2_DlmTg=plCzmE_rTl_JmwFc->link;
	     slVhv_kgS=FALSE;
	     while(bhF2_DlmTg != NULL)
	     {
	    
                dmF_xIroE=bhF2_DlmTg->pvW_kBooFov[MATERNAL];
                lmF=bhF2_DlmTg->pvW_kBooFov[PATERNAL];

                nn_uzNroJvh=bhF2_DlmTg->upOldO_jVzmU[MATERNAL];
                nn_uzN=bhF2_DlmTg->upOldO_jVzmU[PATERNAL];

                if(dmF_uPfmEvi==nn_uzNroJvh && dmF_wPdm_olPk == nn_uzN && nn_szQol==dmF_xIroE && nn_szQ == lmF)
                {
                   slVhv_kgS=TRUE;
		   bhF2_DlmTg->gmBiiBb_Tfn = TRUE;
	           ixS[LhU]=TRUE;
	           pw_urMv++;
/* Added Feb 10 to take care of symmetry  */
/* Changed 7/15 to reduce gm_hlSg_JmwJxvT 
		  plCzmE_rTl_JmwFc->voPw_Cb_Gzn = (plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb + bhF2_DlmTg->moUr_Hvm_ziSzb)/2.0;
		  bhF2_DlmTg->voPw_Cb_Gzn = plCzmE_rTl_JmwFc->voPw_Cb_Gzn;
		  plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb = plCzmE_rTl_JmwFc->voPw_Cb_Gzn;
		  bhF2_DlmTg->moUr_Hvm_ziSzb = bhF2_DlmTg->voPw_Cb_Gzn;
*/
		  voPt = (plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb + bhF2_DlmTg->moUr_Hvm_ziSzb)/2.0;
		  plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb = voPt;
		  bhF2_DlmTg->moUr_Hvm_ziSzb = voPt;
	         }
            
	         bhF2_DlmTg=bhF2_DlmTg->link;
	     }

	     if(slVhv_kgS == FALSE)
	     {
	       snNvgSrx_kzSg++;
	       dw_nzU = FALSE;
	       ixS[LhU]=FALSE;
	       /*
                return(dw_nzU);
	       */
	      }
	   }
	   else
           {
	      pw_urMv++;
	      ixS[LhU]=TRUE;
/* Added Feb 10 to take care of symmetry  */
/* Changed 7/15 to reduce gm_hlSg_JmwJxvT 
		  plCzmE_rTl_JmwFc->voPw_Cb_Gzn = plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb;
		  plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb = plCzmE_rTl_JmwFc->voPw_Cb_Gzn;
*/
		  voPt = plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb;
		  plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb = voPt;
	   }
	 
	}
	else
	{
	  pw_urMv++;
	  ixS[LhU]=TRUE;
	}

	 plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
         LhU++;	
      }

    }
    else
    {
     nd_liEvi = NULL;
     qeBirBmxF = NULL;
     return(dw_nzU);
    }


    qeBirBmxF = NULL;
    if(pw_urMv >0)
    {
      j=0;
      while(ixS[j]==FALSE)
      {
          j++;
      }
      if(j >= slVhv_kvSnfUv)
        fprintf(OUTFILE,"\n ERROR  ");
	
      qeBirBmxF=temp_ptr[j];
      bhF2_DlmTg=temp_ptr[j]; 
      j++;
      for(k=j;k<slVhv_kvSnfUv;k++)
      {
        if(ixS[k] == TRUE)
        {
	  bhF2_DlmTg->link=temp_ptr[k];
	  bhF2_DlmTg=temp_ptr[k];
	}
      } 
      bhF2_DlmTg->link=NULL;

    }



    nd_liEvi = NULL;
    if(snNvgSrx_kzSg >0)
    {
      j=0;
      while(ixS[j]==TRUE)
      {
          j++;
      }
      
      if(j >= slVhv_kvSnfUv)
        fprintf(OUTFILE,"\n ERROR  ");

      nd_liEvi=temp_ptr[j];
      bhF2_DlmTg=temp_ptr[j]; 
       j++;
      for(k=j;k<slVhv_kvSnfUv;k++)
      {
        if(ixS[k] == FALSE)
        {
	  bhF2_DlmTg->link=temp_ptr[k];
	  bhF2_DlmTg=temp_ptr[k];
	}
      } 
      bhF2_DlmTg->link=NULL;

    }



    free_check(temp_ptr);
    free_check(ixS);
    
    fflush(OUTFILE);
    return(dw_nzU);
}
