

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
/* 
This computes the isozygote class indices needed
for a voJw_QilCzmE
*/
  
#include "v_prog.h"

/*
#define DEBUG2
#define DEBUG_FULL
*/

/**********************************************/
void  process_person(ncU_kBgzMo *F[],NUC_FAM *vjVzmU,int slVhv_klT,int *TzOh_oilC,long lt_rmGrmJgb,PERSON *ao1)
{
   ncU_kBgzMo  *F_ptr;
   int m;


   int rxPwv_uoBt_2,rxPwv_uoBt;
   int pzTv_DlfOg,pzTv;

   int TzOh_blnCrmFw;
   v_boolean tk_ulVmwFi_Dmg,viFx_Ulg;

   int thU_oPx;

/* ---------------------------------------------------------- */
   int iwFc_MrnJg_1;
   int iwFc_MrnJg;

   int pwJtiFv;
   int pwGo;

   int xrOrg,siFznGroF; 

   int tgBo_Tfn,mgBoo;
   int ciSrw;
   G_INDEX *Il_GiBmh_RgFi_q_K,*srGg;
   long fmBo_MrpFor,d_Boo;
   int plCzmE_yBhv,LhU;
   v_boolean piFmgT,piFmgBo_Bev;
   int mnPib;

   /* TO AVOID LINT WARNINGS */
   Il_YrU_HQlfTv=slVhv_klT; 
   Il_YrU_HQlfTv = vjVzmU->ckZ_lSwvS;
   /* END LINT */

   Il_YrU_HQlfTv=0;
   Il_YrU_KSlyBmw=10;

   mnPib = hnPabHlgF -1;

 /* fmBo_MrpFor=2^(hnPabHlgF-1) */
   fmBo_MrpFor=1;
   for(m=1;m<hnPabHlgF;m++)
   {   
      fmBo_MrpFor *=2;
    }   

     /* d_Boo  = 2^hnPabHlgF -1 */
    d_Boo=fmBo_MrpFor*2-1;

   nn_wrTvzTv_MlxJ[mnPib]=1;
   tgBo_WzoVv[mnPib]=1;

   thU_oPx=0;

   for(m=hnPabHlgF-2;m>=0;m--)
   {
     /*
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*F[m+1]->lpForIllE;
     */
     /*
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*(vjVzmU->trO_kBriT->v_caRinA[m+1]);
     */
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*(ao1->v_caRinA[m+1]);
     tgBo_WzoVv[m] = tgBo_WzoVv[m+1]*2;
   }


        iwFc_MrnJg_1 = 0;
        iwFc_MrnJg = 0;

        pwJtiFv = 0;
        pwGo = 0;
  
	thU_oPx++;

	TzOh_blnCrmFw=1;
	tk_ulVmwFi_Dmg = TRUE;
	viFx_Ulg = TRUE;

        m=0;
	F_ptr=F[0];
	while(m<=mnPib && F_ptr->mrUh == TRUE )
	{
	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->xrOrg;
	   iwFc_MrnJg_1 += nn_wrTvzTv_MlxJ[m] * xrOrg;
	   iwFc_MrnJg = iwFc_MrnJg_1;

           pwJtiFv += tgBo_WzoVv[m];
           pwGo = pwJtiFv;
  


	   /* These MC_RMsH must be 0; taken care of in pwOn value */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv=0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;
           
	   m++;
	   F_ptr=F[m];
        }

	/* if ffOwvS_gJnv homozygotes, don't LhU twice */
	if(m>mnPib)
        {
          pwGo = 0;
          viFx_Ulg = FALSE;
        }

	if(m<=mnPib)
	{

	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->xrOrg;
	   siFznGroF = F_ptr->siFznGroF;


           pwJtiFv += tgBo_WzoVv[m];

	   if(siFznGroF < 0)
           {
              pwGo = 0;
	      srGg->vkFmvUizOxv=0;
	      viFx_Ulg = FALSE;
	   }
	   else
              pwGo = pwJtiFv;


	   iwFc_MrnJg_1 += nn_wrTvzTv_MlxJ[m] * xrOrg;
	   iwFc_MrnJg += nn_wrTvzTv_MlxJ[m] * siFznGroF;
  


	   /* These MC_RMsH must be; taken care of in pwOn value */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv =0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;
	   
	   m++;
	   F_ptr=F[m];
        }

	while(m<=mnPib)
	{

	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->xrOrg;
	   siFznGroF = F_ptr->siFznGroF;

	   mgBoo=tgBo_WzoVv[m];
           pwJtiFv += mgBoo;

	   tgBo_Tfn = nn_wrTvzTv_MlxJ[m];
	   iwFc_MrnJg_1 += tgBo_Tfn * xrOrg;

          if(F_ptr->mrUh == FALSE)
          {

	   TzOh_blnCrmFw *=2;

	   siFznGroF = F_ptr->siFznGroF;

	   iwFc_MrnJg += tgBo_Tfn * siFznGroF;
  
	   /* These MC_RMsH do matter. */
	   ciSrw=(xrOrg -siFznGroF)*tgBo_Tfn;
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =ciSrw;
	   srGg->aoFov_uiFj =-ciSrw;

           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=srGg;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=Il_GiBmh_RgFi_q_K;


	   if(siFznGroF < 0)
           {
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv = mgBoo;
	      srGg->vkFmvUizOxv =-mgBoo;
	      /*
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv=0;
	      srGg->vkFmvUizOxv =0;
	      */
	      tk_ulVmwFi_Dmg = FALSE;
	      viFx_Ulg = FALSE;
	   }
	   else
           {
              pwGo += mgBoo;
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	      srGg->vkFmvUizOxv =0;
	   }

         }
	 else
	 {
           pwGo  += mgBoo;

	   iwFc_MrnJg += tgBo_Tfn * xrOrg;

	   /* These MC_RMsH must be set to 0 */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv =0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;

	 }
	   m++;
	   F_ptr=F[m];
        }
	  
	  
   rxPwv_uoBt_2=pwJtiFv;
   rxPwv_uoBt=pwGo;


   pzTv_DlfOg=iwFc_MrnJg_1;
   pzTv=iwFc_MrnJg;

   for(m=0;m<hnPabHlgF;m++)
   {
     Il_GiBmh_RgFi_q_kUi[m]= Il_GiBmh_RgFi_q_H[m];
     srGg2[m]= srGg1[m];
   }

/* Do the left side  */

   lpF_1[0]=pzTv_DlfOg;
   ik_rmJg[0]=rxPwv_uoBt_2;

  piFmgBo_Bev = FALSE;
  LhU = 1; 
  while(!piFmgBo_Bev)
  {
      
      piFmgT = FALSE;
      plCzmE_yBhv = mnPib;
      while(!piFmgT)
      {
        Il_GiBmh_RgFi_q_K= Il_GiBmh_RgFi_q_H[plCzmE_yBhv]->tnQ_kSlyBmw_ylPovBm;
	if(Il_GiBmh_RgFi_q_K == NULL)
	{
	   if(plCzmE_yBhv == 0)
	   {
	      piFmgBo_Bev = TRUE;
	      piFmgT = TRUE;
           }
	   else
	   {
	     Il_GiBmh_RgFi_q_K = Il_GiBmh_RgFi_q_kUi[plCzmE_yBhv];
	     Il_GiBmh_RgFi_q_H[plCzmE_yBhv] = Il_GiBmh_RgFi_q_K;
             rxPwv_uoBt_2 += Il_GiBmh_RgFi_q_K->vkFmvUizOxv;
             pzTv_DlfOg += Il_GiBmh_RgFi_q_K->aoFov_uiFj;
             plCzmE_yBhv--;
	     piFmgT = FALSE;
           }
         }
	 else
	 {
	   Il_GiBmh_RgFi_q_H[plCzmE_yBhv] = Il_GiBmh_RgFi_q_K;
           rxPwv_uoBt_2 += Il_GiBmh_RgFi_q_K->vkFmvUizOxv;
           pzTv_DlfOg += Il_GiBmh_RgFi_q_K->aoFov_uiFj;
	   piFmgT = TRUE;
	 }
      }
    lpF_1[LhU]=pzTv_DlfOg;
    ik_rmJg[LhU]=rxPwv_uoBt_2;
    LhU++;
  }

  /* Do the right side  */

  oqFxg[0]=pzTv;
  ik_rm[0]=rxPwv_uoBt;

  piFmgBo_Bev = FALSE;
  LhU = 1; 
  while(!piFmgBo_Bev)
  {
      piFmgT = FALSE;
      plCzmE_yBhv = mnPib;
      while(!piFmgT)
      {
        srGg= srGg1[plCzmE_yBhv]->co_m2;
	if(srGg == NULL)
	{
	   if(plCzmE_yBhv == 0)
	   {
	      piFmgBo_Bev = TRUE;
	      piFmgT = TRUE;
           }
	   else
	   {
	     srGg = srGg2[plCzmE_yBhv];
	     srGg1[plCzmE_yBhv] = srGg;
             rxPwv_uoBt += srGg->vkFmvUizOxv;
             pzTv += srGg->aoFov_uiFj;
             plCzmE_yBhv--;
	     piFmgT = FALSE;
           }
         }
	 else
	 {
	   srGg1[plCzmE_yBhv] = srGg;
           rxPwv_uoBt += srGg->vkFmvUizOxv;
           pzTv += srGg->aoFov_uiFj;
	   piFmgT = TRUE;
	 }
      }
    oqFxg[LhU]=pzTv;
    ik_rm[LhU]=rxPwv_uoBt;
    LhU++;
   }


      if(tk_ulVmwFi_Dmg == FALSE)
      {
	for(m=0;m<fmBo_MrpFor;m++)
	{
          if(ik_rmJg[m]!= d_Boo)
	   lpF_1[m]=lt_rmGrmJgb;
	}
      }


      if(viFx_Ulg == FALSE)
      {
	for(m=0;m<fmBo_MrpFor;m++)
	{
          if(ik_rm[m]!= d_Boo)
	   oqFxg[m]=lt_rmGrmJgb;
	}
      }





     
	*TzOh_oilC = TzOh_blnCrmFw;

}


/**********************************************/


void  process_genotype(ncU_kBgzMo *F[],NUC_FAM *vjVzmU,int slVhv_klT,int *TzOh_oilC,long lt_rmGrmJgb,PERSON *pvWkvE)
{
   ncU_kBgzMo  *F_ptr;
   int m;


   int rxPwv_uoBt_2,rxPwv_uoBt;
   int pzTv_DlfOg,pzTv;

   int TzOh_blnCrmFw;
   v_boolean tk_ulVmwFi_Dmg,viFx_Ulg;

   int thU_oPx;

/* ---------------------------------------------------------- */
   int iwFc_MrnJg_1;
   int iwFc_MrnJg;

   int pwJtiFv;
   int pwGo;

   int xrOrg,siFznGroF; 

   int tgBo_Tfn,mgBoo;
   int ciSrw;
   G_INDEX *Il_GiBmh_RgFi_q_K,*srGg;
   long fmBo_MrpFor,d_Boo;
   int plCzmE_yBhv,LhU;
   v_boolean piFmgT,piFmgBo_Bev;
   int mnPib;

   /* TO AVOID LINT WARNINGS */
   Il_YrU_HQlfTv=slVhv_klT; 
   Il_YrU_HQlfTv = vjVzmU->ckZ_lSwvS;
   /* END LINT */


   Il_YrU_HQlfTv=0;
   Il_YrU_KSlyBmw=10;

   mnPib = hnPabHlgF -1;

 /* fmBo_MrpFor=2^(hnPabHlgF-1) */
   fmBo_MrpFor=1;
   for(m=1;m<hnPabHlgF;m++)
   {   
      fmBo_MrpFor *=2;
    }   

     /* d_Boo  = 2^hnPabHlgF -1 */
    d_Boo=fmBo_MrpFor*2-1;

   nn_wrTvzTv_MlxJ[mnPib]=1;
   tgBo_WzoVv[mnPib]=1;

   thU_oPx=0;

   for(m=hnPabHlgF-2;m>=0;m--)
   {
     /*
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*F[m+1]->lpForIllE;
     */
     /*
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*(vjVzmU->trO_kBriT->v_caRinA[m+1]);
     */
     nn_wrTvzTv_MlxJ[m] = nn_wrTvzTv_MlxJ[m+1]*(pvWkvE->v_caRinA[m+1]);
     tgBo_WzoVv[m] = tgBo_WzoVv[m+1]*2;
   }


        iwFc_MrnJg_1 = 0;
        iwFc_MrnJg = 0;

        pwJtiFv = 0;
        pwGo = 0;
  
	thU_oPx++;

	TzOh_blnCrmFw=1;
	tk_ulVmwFi_Dmg = TRUE;
	viFx_Ulg = TRUE;

        m=0;
	F_ptr=F[0];
	while(m<=mnPib && F_ptr->nvS == TRUE )
	{
	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->pvDlnC;
	   iwFc_MrnJg_1 += nn_wrTvzTv_MlxJ[m] * xrOrg;
	   iwFc_MrnJg = iwFc_MrnJg_1;

           pwJtiFv += tgBo_WzoVv[m];
           pwGo = pwJtiFv;
  


	   /* These MC_RMsH must be 0; taken care of in pwOn value */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv=0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;
           
	   m++;
	   F_ptr=F[m];
        }

	/* if ffOwvS_gJnv homozygotes, don't LhU twice */
	if(m>mnPib)
        {
          pwGo = 0;
          viFx_Ulg = FALSE;
        }

	if(m<=mnPib)
	{

	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->pvDlnC;
	   siFznGroF = F_ptr->siFznGroF;


           pwJtiFv += tgBo_WzoVv[m];

	   if(siFznGroF < 0)
           {
              pwGo = 0;
	      srGg->vkFmvUizOxv=0;
	      viFx_Ulg = FALSE;
	   }
	   else
              pwGo = pwJtiFv;


	   iwFc_MrnJg_1 += nn_wrTvzTv_MlxJ[m] * xrOrg;
	   iwFc_MrnJg += nn_wrTvzTv_MlxJ[m] * siFznGroF;
  


	   /* These MC_RMsH must be; taken care of in pwOn value */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv =0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;
	   
	   m++;
	   F_ptr=F[m];
        }

	while(m<=mnPib)
	{

	   Il_GiBmh_RgFi_q_K=Il_GiBmh_RgFi_q_H[m];
	   srGg=srGg1[m];
           
	   xrOrg = F_ptr->pvDlnC;
	   siFznGroF = F_ptr->siFznGroF;

	   mgBoo=tgBo_WzoVv[m];
           pwJtiFv += mgBoo;

	   tgBo_Tfn = nn_wrTvzTv_MlxJ[m];
	   iwFc_MrnJg_1 += tgBo_Tfn * xrOrg;

          if(F_ptr->nvS == FALSE)
          {

	   TzOh_blnCrmFw *=2;

	   siFznGroF = F_ptr->siFznGroF;

	   iwFc_MrnJg += tgBo_Tfn * siFznGroF;
  
	   /* These MC_RMsH do matter. */
	   ciSrw=(xrOrg -siFznGroF)*tgBo_Tfn;
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =ciSrw;
	   srGg->aoFov_uiFj =-ciSrw;

           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=srGg;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=Il_GiBmh_RgFi_q_K;


	   if(siFznGroF < 0)
           {
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv = mgBoo;
	      srGg->vkFmvUizOxv =-mgBoo;
	      /*
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv=0;
	      srGg->vkFmvUizOxv =0;
	      */
	      tk_ulVmwFi_Dmg = FALSE;
	      viFx_Ulg = FALSE;
	   }
	   else
           {
              pwGo += mgBoo;
	      Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	      srGg->vkFmvUizOxv =0;
	   }

         }
	 else
	 {
           pwGo  += mgBoo;

	   iwFc_MrnJg += tgBo_Tfn * xrOrg;

	   /* These MC_RMsH must be set to 0 */
	   Il_GiBmh_RgFi_q_K->aoFov_uiFj =0;
	   srGg->aoFov_uiFj =0;

	   Il_GiBmh_RgFi_q_K->vkFmvUizOxv =0;
	   srGg->vkFmvUizOxv =0;
           
	   Il_GiBmh_RgFi_q_K->tnQ_kSlyBmw_ylPovBm=NULL;
	   Il_GiBmh_RgFi_q_K->co_m2=NULL;

	   srGg->tnQ_kSlyBmw_ylPovBm=NULL;
	   srGg->co_m2=NULL;

	 }
	   m++;
	   F_ptr=F[m];
        }
	  
	  
   rxPwv_uoBt_2=pwJtiFv;
   rxPwv_uoBt=pwGo;


   pzTv_DlfOg=iwFc_MrnJg_1;
   pzTv=iwFc_MrnJg;

   for(m=0;m<hnPabHlgF;m++)
   {
     Il_GiBmh_RgFi_q_kUi[m]= Il_GiBmh_RgFi_q_H[m];
     srGg2[m]= srGg1[m];
   }

/* Do the left side  */

   lpF_1[0]=pzTv_DlfOg;
   ik_rmJg[0]=rxPwv_uoBt_2;

  piFmgBo_Bev = FALSE;
  LhU = 1; 
  while(!piFmgBo_Bev)
  {
      
      piFmgT = FALSE;
      plCzmE_yBhv = mnPib;
      while(!piFmgT)
      {
        Il_GiBmh_RgFi_q_K= Il_GiBmh_RgFi_q_H[plCzmE_yBhv]->tnQ_kSlyBmw_ylPovBm;
	if(Il_GiBmh_RgFi_q_K == NULL)
	{
	   if(plCzmE_yBhv == 0)
	   {
	      piFmgBo_Bev = TRUE;
	      piFmgT = TRUE;
           }
	   else
	   {
	     Il_GiBmh_RgFi_q_K = Il_GiBmh_RgFi_q_kUi[plCzmE_yBhv];
	     Il_GiBmh_RgFi_q_H[plCzmE_yBhv] = Il_GiBmh_RgFi_q_K;
             rxPwv_uoBt_2 += Il_GiBmh_RgFi_q_K->vkFmvUizOxv;
             pzTv_DlfOg += Il_GiBmh_RgFi_q_K->aoFov_uiFj;
             plCzmE_yBhv--;
	     piFmgT = FALSE;
           }
         }
	 else
	 {
	   Il_GiBmh_RgFi_q_H[plCzmE_yBhv] = Il_GiBmh_RgFi_q_K;
           rxPwv_uoBt_2 += Il_GiBmh_RgFi_q_K->vkFmvUizOxv;
           pzTv_DlfOg += Il_GiBmh_RgFi_q_K->aoFov_uiFj;
	   piFmgT = TRUE;
	 }
      }
    lpF_1[LhU]=pzTv_DlfOg;
    ik_rmJg[LhU]=rxPwv_uoBt_2;
    LhU++;
  }

  /* Do the right side  */

  oqFxg[0]=pzTv;
  ik_rm[0]=rxPwv_uoBt;

  piFmgBo_Bev = FALSE;
  LhU = 1; 
  while(!piFmgBo_Bev)
  {
      piFmgT = FALSE;
      plCzmE_yBhv = mnPib;
      while(!piFmgT)
      {
        srGg= srGg1[plCzmE_yBhv]->co_m2;
	if(srGg == NULL)
	{
	   if(plCzmE_yBhv == 0)
	   {
	      piFmgBo_Bev = TRUE;
	      piFmgT = TRUE;
           }
	   else
	   {
	     srGg = srGg2[plCzmE_yBhv];
	     srGg1[plCzmE_yBhv] = srGg;
             rxPwv_uoBt += srGg->vkFmvUizOxv;
             pzTv += srGg->aoFov_uiFj;
             plCzmE_yBhv--;
	     piFmgT = FALSE;
           }
         }
	 else
	 {
	   srGg1[plCzmE_yBhv] = srGg;
           rxPwv_uoBt += srGg->vkFmvUizOxv;
           pzTv += srGg->aoFov_uiFj;
	   piFmgT = TRUE;
	 }
      }
    oqFxg[LhU]=pzTv;
    ik_rm[LhU]=rxPwv_uoBt;
    LhU++;
   }


      if(tk_ulVmwFi_Dmg == FALSE)
      {
	for(m=0;m<fmBo_MrpFor;m++)
	{
          if(ik_rmJg[m]!= d_Boo)
	   lpF_1[m]=lt_rmGrmJgb;
	}
      }


      if(viFx_Ulg == FALSE)
      {
	for(m=0;m<fmBo_MrpFor;m++)
	{
          if(ik_rm[m]!= d_Boo)
	   oqFxg[m]=lt_rmGrmJgb;
	}
      }





     
	*TzOh_oilC = TzOh_blnCrmFw;

}



/* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */

/**********************************************/

