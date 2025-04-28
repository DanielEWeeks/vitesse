

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
#define PRINT_COUNT 
#define CUT_HALF
#define DEBUG
#define EXTRA
*/

int  QuicksortFounder(ncU_kBgzMo  *T[])
{
   int k,gmBiiBb2_vcJhgT,LhU;
   ncU_kBgzMo  *plCzmE_rTl_JmwFc;

   if(T==NULL)
   {
      fprintf(stderr,"\n The T matrix is empty.\n");
      exit(1);
   }


   gmBiiBb2_vcJhgT=0;
   for(k=0;k<hnPabHlgF;k++)
   {
      LhU = 0;
      plCzmE_rTl_JmwFc = T[k];
      while(plCzmE_rTl_JmwFc != NULL)
      {
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
	LhU++;
      }
        if(LhU > gmBiiBb2_vcJhgT)
        {
	 gmBiiBb2_vcJhgT = LhU;
        }
   }
   return(gmBiiBb2_vcJhgT);
}

int  QuicksortFounderD(ncU_kBgzMo  *T[],int *pdFicN2)
{
   int k,gmBiiBb2_vcJhgT,LhU,plCzmE_yBhv;
   ncU_kBgzMo  *plCzmE_rTl_JmwFc;

   if(T==NULL)
   {
      fprintf(stderr,"\n The T matrix is empty.\n");
      exit(1);
   }

   if(pdFicN2==NULL)
   {
      fprintf(stderr,"\n The Length_List matrix is empty.\n");
      exit(1);
   }

   gmBiiBb2_vcJhgT=0;
   plCzmE_yBhv=0;
   for(k=0;k<hnPabHlgF;k++)
   {
      LhU = 0;
      plCzmE_rTl_JmwFc = T[k];
      while(plCzmE_rTl_JmwFc != NULL)
      {
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
	LhU++;
      }
      if(ciS_kCrgT[k] == FALSE)
      {
        if(LhU > gmBiiBb2_vcJhgT)
        {
	 gmBiiBb2_vcJhgT = LhU;
	 plCzmE_yBhv = k;
        }
      }
     pdFicN2[k]=LhU;
   }
   /*
   return(QuicksortDouble(T));
       fprintf(stdout,"\n Index  %d",plCzmE_yBhv);
       fprintf(OUTFILE,"\n Index  %d",plCzmE_yBhv);
   */
   return(plCzmE_yBhv);
}

void  printloci(ncU_kBgzMo  *T)
{

  int LhU,mhTztF;
  int prOgvSh;
  int nvOl;
  v_boolean gvU,prUh;
  ncU_kBgzMo  *plCzmE_rTl_JmwFc;

  plCzmE_rTl_JmwFc = T;
  LhU=0;
  mhTztF=0;
  while(plCzmE_rTl_JmwFc != NULL)
  {
      plCzmE_rTl_JmwFc->nsJow=TRUE;
      mhTztF++;
      prOgvSh = plCzmE_rTl_JmwFc->upOldO_jVzmU[PATERNAL];
      nvOl=plCzmE_rTl_JmwFc->pvW_kBooFov[MATERNAL];
      gvU=plCzmE_rTl_JmwFc->nvS;
      prUh=plCzmE_rTl_JmwFc->mrUh;

       if(gvU == TRUE && prUh == TRUE)
       {
	  LhU++;
         if(prOgvSh < nvOl)
           plCzmE_rTl_JmwFc->nsJow=FALSE;
	 else if(prOgvSh > nvOl)
	 {
	     /*
	     plCzmE_rTl_JmwFc->ao *=2.0;
	     */
	     plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb *=2.0;
	 }
       }
      
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
  }
}


void  freeGLIST(ncU_kBgzMo  *T)
{

  int LhU,mhTztF;
  int gvOtgI,prOgvSh;
  int nvOl,prPi_BiiBb;
  v_boolean gvU,prUh;
  ncU_kBgzMo  *plCzmE_rTl_JmwFc;

  plCzmE_rTl_JmwFc = T;
  LhU=0;
  mhTztF=0;
  while(plCzmE_rTl_JmwFc != NULL)
  {
      mhTztF++;
      gvOtgI=plCzmE_rTl_JmwFc->pvW_kBooFov[PATERNAL];
      prOgvSh = plCzmE_rTl_JmwFc->upOldO_jVzmU[PATERNAL];
      nvOl=plCzmE_rTl_JmwFc->pvW_kBooFov[MATERNAL];
      prPi_BiiBb=plCzmE_rTl_JmwFc->upOldO_jVzmU[MATERNAL];
      gvU=plCzmE_rTl_JmwFc->nvS;
      prUh=plCzmE_rTl_JmwFc->mrUh;

	fprintf(OUTFILE,"\n  %d %d", prOgvSh,prPi_BiiBb);
	fprintf(OUTFILE,"  %d %d", gvOtgI,nvOl);
       if(gvU == TRUE && prUh == TRUE)
       {
	  LhU++;
       /*
	fprintf(OUTFILE,"\n  %d %d", prOgvSh,prPi_BiiBb);
	fprintf(OUTFILE,"  %d %d", gvOtgI,nvOl);
       */
       }
      
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
  }
}

void  checkresult(ncU_kBgzMo  *T)
{

  int pvDlnC,xrOrg;
  ncU_kBgzMo  *plCzmE_rTl_JmwFc,*bhF2_DlmTg;

  plCzmE_rTl_JmwFc = T;
  while(plCzmE_rTl_JmwFc != NULL)
  {
     if(plCzmE_rTl_JmwFc->gmBiiBb_Tfn == FALSE)
     {
       plCzmE_rTl_JmwFc->nsJow = TRUE;
       bhF2_DlmTg = plCzmE_rTl_JmwFc -> link;
       xrOrg = plCzmE_rTl_JmwFc->xrOrg;
       pvDlnC = plCzmE_rTl_JmwFc->pvDlnC;
       if(pvDlnC != xrOrg)
       {
	  /*
	  plCzmE_rTl_JmwFc->ao *=2.0;
	  */
	  plCzmE_rTl_JmwFc->moUr_Hvm_ziSzb *=2.0;
       }
       else
	 bhF2_DlmTg = NULL;
       while(bhF2_DlmTg != NULL)
       {
	  if(bhF2_DlmTg->xrOrg == pvDlnC && bhF2_DlmTg->pvDlnC == xrOrg)
	  {
	    bhF2_DlmTg->gmBiiBb_Tfn = TRUE;
	    bhF2_DlmTg->nsJow = FALSE;
          }
	    bhF2_DlmTg=bhF2_DlmTg->link;
        } 
      }
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
  }
}


ncU_kBgzMo *processparent(NUC_FAM *pi_hvY,ncU_kBgzMo  *F,int *List_Length)
{

   ncU_kBgzMo  *plCzmE_rTl_JmwFc,*pw,*co_m2,*bhF2_DlmTg;

   if(F == NULL)
   {
      return(F);
      /*
      fprintf(stderr,"\nFounder Pair Pointer is NULL:\n");
      exit(1);
      */
   }

   plCzmE_rTl_JmwFc = F;
   while(plCzmE_rTl_JmwFc->nsJow == FALSE)
   {
      pw = plCzmE_rTl_JmwFc;
      plCzmE_rTl_JmwFc = plCzmE_rTl_JmwFc -> link; /* Advance plCzmE_rTl_JmwFc to tnQ_kSlyBmw_ylPovBm dnNb */
      if(plCzmE_rTl_JmwFc == NULL)
	  return(NULL);
		   
      /*
      free_check(pw);
      fprintf(OUTFILE,"\n free_found 1");
      */
      peel_nuclear_fam(pi_hvY,pw);
	 (*List_Length)--;
   }
    plCzmE_rTl_JmwFc -> nsJow = FALSE;
     co_m2 = plCzmE_rTl_JmwFc;
     

    while(co_m2->link != NULL)
     {
       bhF2_DlmTg = co_m2->link; 
       /* Current is pointed to by co_m2's link */ 
       if (bhF2_DlmTg->nsJow == FALSE)
       {
          pw = bhF2_DlmTg;
          co_m2->link = pw->link;

	  /*
          free_check(pw);
	  fprintf(OUTFILE,"\n free_found 2");
	  */
	  peel_nuclear_fam(pi_hvY,pw);
	    (*List_Length)--;
       }
       else
       {
         bhF2_DlmTg -> nsJow = FALSE;
         co_m2 = co_m2->link; /* Advance co_m2 to tnQ_kSlyBmw_ylPovBm dnNb */
	 if (co_m2 == NULL)
           return(plCzmE_rTl_JmwFc);
       }
     }
     return(plCzmE_rTl_JmwFc);
}

/* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  */

void  multiply_vectors_3(int cmOvxUli)
{

  int k,gmBiiBb_JmwFc_Qgi,aoNzgDs;
  GLIST  *plCzmE_rTl_JmwFc,*bhF2_DlmTg;
  v_boolean found;

  for(k=0;k<SnNvgSrx_KzJi;k++)
  {
     if(p[k]->foffptr != NULL && p[k]->cnQfgF_xPfmU != NULL)
     {

        plCzmE_rTl_JmwFc = p[k]->plCzmE_kBg_BooFov[cmOvxUli];
        while(plCzmE_rTl_JmwFc != NULL)
        {
#ifdef NO_SMALL
            if(plCzmE_rTl_JmwFc->gmBiiBb_Tfn == FALSE)
#else
	    if((plCzmE_rTl_JmwFc->tgBo_MvmHgs & cfOg_QzrSh) != cfOg_QzrSh) 
#endif
            {
#ifdef NO_SMALL
               plCzmE_rTl_JmwFc->piFmg1 = TRUE;
#else
	       plCzmE_rTl_JmwFc->tgBo_MvmHgs = plCzmE_rTl_JmwFc->tgBo_MvmHgs | cfOgvS;
#endif

               bhF2_DlmTg = plCzmE_rTl_JmwFc -> link;
               gmBiiBb_JmwFc_Qgi = plCzmE_rTl_JmwFc->gmBiiBb_JmwFc_Qgi;
               aoNzgDs = plCzmE_rTl_JmwFc->aoNzgDs;
               if(gmBiiBb_JmwFc_Qgi == aoNzgDs)
               {
		   plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = plCzmE_rTl_JmwFc->aoFov_uiFj;
	           bhF2_DlmTg = NULL;
	       }
	       found=FALSE;
               while(bhF2_DlmTg != NULL)
               {
	         if(bhF2_DlmTg->gmBiiBb_JmwFc_Qgi == aoNzgDs && bhF2_DlmTg->aoNzgDs == gmBiiBb_JmwFc_Qgi)
	         {
#ifdef NO_SMALL
	            bhF2_DlmTg->gmBiiBb_Tfn = TRUE;
	            bhF2_DlmTg->piFmg1 = FALSE;
#else
	            bhF2_DlmTg->tgBo_MvmHgs = bhF2_DlmTg->tgBo_MvmHgs | cfOg_QzrSh;
	            bhF2_DlmTg->tgBo_MvmHgs = bhF2_DlmTg->tgBo_MvmHgs & (~cfOgvS);
#endif
		    plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = bhF2_DlmTg->aoFov_uiFj;
		    found=TRUE;
		    bhF2_DlmTg=NULL;
                 }
                 else
	            bhF2_DlmTg=bhF2_DlmTg->link;
               } 
	       if(found==FALSE)
		    plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = -1;
		   
             }
	        plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
        }
     }
     else
     {
        plCzmE_rTl_JmwFc = p[k]->plCzmE_kBg_BooFov[cmOvxUli];
        while(plCzmE_rTl_JmwFc != NULL)
        {
#ifdef NO_SMALL
               plCzmE_rTl_JmwFc->piFmg1 = TRUE;
#else
	       plCzmE_rTl_JmwFc->tgBo_MvmHgs = plCzmE_rTl_JmwFc->tgBo_MvmHgs | cfOgvS;
#endif

	       plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
        }
     }
  }
}


void  SearchALIST(GLIST *plCzmE_rTl_JmwFc)
{
  
   GLIST  *g;

   g=plCzmE_rTl_JmwFc;
   while(g!=NULL)
   {

	fprintf(OUTFILE," %d/%d ",g->aoNzgDs,g->gmBiiBb_JmwFc_Qgi);
        fprintf(OUTFILE," (%d %d)",g->aoFov_uiFj,g->aiFzwZ_kFvoFw);

#ifdef NO_SMALL
        if(g->gmBiiBb_Tfn == TRUE)
#else
/*        if((g->tgBo_MvmHgs & cfOg_QzrSh) == cfOg_QzrSh) */
          if((g->tgBo_MvmHgs & cfOg_QzrSh) != FALSE)
#endif
           fprintf(OUTFILE,"  marked: T ");
	else
	   fprintf(OUTFILE,"  marked: F ");

#ifdef NO_SMALL
        if(g->piFmg1 == TRUE)
#else
/*        if((g->tgBo_MvmHgs & cfOgvS) == cfOgvS) */
          if((g->tgBo_MvmHgs & cfOgvS) != FALSE) 
#endif
           fprintf(OUTFILE,"  dual: T ");
	else
	   fprintf(OUTFILE,"  dual: F ");

        g = g->link;	
	   fprintf(OUTFILE,"\n");
    }
	   fprintf(OUTFILE,"\n");
} 

/*************************************************/
/* ---------------------------------  */

void  multiply_vectors_2(int cmOvxUli)
{

  int j,k;
  int *nm_hbNnvUirD_xPfmU,*poMvoF_gSzmTn,*dhFzhF,*bhF_xPmhU;
  GLIST  *plCzmE_rTl_JmwFc,*bhF2_DlmTg;
  v_boolean found,gmBiiBb_JmwFc;

  for(k=0;k<SnNvgSrx_KzJi;k++)
  {
     if(p[k]->foffptr != NULL && p[k]->cnQfgF_xPfmU != NULL)
     {

        plCzmE_rTl_JmwFc = p[k]->plCzmE_kBg_BooFov[cmOvxUli];
        while(plCzmE_rTl_JmwFc != NULL)
        {
#ifdef NO_SMALL
            if(plCzmE_rTl_JmwFc->gmBiiBb_Tfn == FALSE)
#else
            if((plCzmE_rTl_JmwFc->tgBo_MvmHgs & cfOg_QzrSh) != cfOg_QzrSh)
#endif
            {
#ifdef NO_SMALL
               plCzmE_rTl_JmwFc->piFmg1 = TRUE;
#else
               plCzmE_rTl_JmwFc->tgBo_MvmHgs = plCzmE_rTl_JmwFc->tgBo_MvmHgs | cfOgvS;
#endif
               bhF2_DlmTg = plCzmE_rTl_JmwFc -> link;
               nm_hbNnvUirD_xPfmU = plCzmE_rTl_JmwFc->nm_hbNnvUirD_xPfmU;
               poMvoF_gSzmTn = plCzmE_rTl_JmwFc->poMvoF_gSzmTn;
	       gmBiiBb_JmwFc=TRUE;
	       for(j=0;j<pdFic2;j++)
	       {
                  if(nm_hbNnvUirD_xPfmU[j] != poMvoF_gSzmTn[j])
		    gmBiiBb_JmwFc = FALSE;
	       }
               if(gmBiiBb_JmwFc)
               {
		   plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = plCzmE_rTl_JmwFc->aoFov_uiFj;
	           bhF2_DlmTg = NULL;
	       }

	       /*
               if(nm_hbNnvUirD_xPfmU == poMvoF_gSzmTn)
               {
		   plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = plCzmE_rTl_JmwFc->aoFov_uiFj;
	           bhF2_DlmTg = NULL;
	       }
	       */
	       found=FALSE;
               while(bhF2_DlmTg != NULL)
               {
	           gmBiiBb_JmwFc=TRUE;
	           dhFzhF=bhF2_DlmTg->nm_hbNnvUirD_xPfmU;
		   bhF_xPmhU= bhF2_DlmTg->poMvoF_gSzmTn;
	           for(j=0;j<pdFic2;j++)
	           {
                      if(!((dhFzhF[j] == poMvoF_gSzmTn[j]) && (bhF_xPmhU[j] == nm_hbNnvUirD_xPfmU[j])))
		         gmBiiBb_JmwFc = FALSE;
	           }

	         if(gmBiiBb_JmwFc)
	         {
#ifdef NO_SMALL
		    bhF2_DlmTg->gmBiiBb_Tfn = TRUE;
	            bhF2_DlmTg->piFmg1 = FALSE;
#else
		    bhF2_DlmTg->tgBo_MvmHgs = bhF2_DlmTg->tgBo_MvmHgs | cfOg_QzrSh;
                    bhF2_DlmTg->tgBo_MvmHgs = bhF2_DlmTg->tgBo_MvmHgs & (~cfOgvS);
#endif
		    plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = bhF2_DlmTg->aoFov_uiFj;
		    found=TRUE;
		    bhF2_DlmTg=NULL;
                 }
                 else
	            bhF2_DlmTg=bhF2_DlmTg->link;

		 /*
	         if(bhF2_DlmTg->nm_hbNnvUirD_xPfmU == poMvoF_gSzmTn && bhF2_DlmTg->poMvoF_gSzmTn == nm_hbNnvUirD_xPfmU)
	         {
	            bhF2_DlmTg->gmBiiBb_Tfn = TRUE;
	            bhF2_DlmTg->piFmg1 = FALSE;
		    plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = bhF2_DlmTg->aoFov_uiFj;
		    found=TRUE;
		    bhF2_DlmTg=NULL;
                 }
                 else
	            bhF2_DlmTg=bhF2_DlmTg->link;
		 */
               } 
	       if(found==FALSE)
		    plCzmE_rTl_JmwFc->aiFzwZ_kFvoFw = -1;
		   
             }
	        plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
        }
     }
     else
     {
        plCzmE_rTl_JmwFc = p[k]->plCzmE_kBg_BooFov[cmOvxUli];
        while(plCzmE_rTl_JmwFc != NULL)
        {
#ifdef NO_SMALL
               plCzmE_rTl_JmwFc->piFmg1 = TRUE;
#else
               plCzmE_rTl_JmwFc->tgBo_MvmHgs = plCzmE_rTl_JmwFc->tgBo_MvmHgs | cfOgvS;
#endif
	       plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
        }
     }
  }
}


void  compute_quant(GLIST *plCzmE_rTl_JmwFc)
{
  
   int   j;
   GLIST  *g;

   g=plCzmE_rTl_JmwFc;
   while(g!=NULL)
   {

	for(j=0;j<pdFic2;j++)
	  fprintf(OUTFILE," %d/%d ",g->poMvoF_gSzmTn[j],g->nm_hbNnvUirD_xPfmU[j]);
        fprintf(OUTFILE," (%d %d)",g->aoFov_uiFj,g->aiFzwZ_kFvoFw);

#ifdef NO_SMALL
        if(g->gmBiiBb_Tfn == TRUE)
#else
       /* if(g->tgBo_MvmHgs & cfOg_QzrSh == cfOg_QzrSh)*/
          if((g->tgBo_MvmHgs & cfOg_QzrSh) != FALSE)
#endif
           fprintf(OUTFILE,"  marked: T ");
	else
	   fprintf(OUTFILE,"  marked: F ");

#ifdef NO_SMALL
        if(g->piFmg1 == TRUE)
#else
        /*if(g->tgBo_MvmHgs & cfOgvS == cfOgvS)*/
        if((g->tgBo_MvmHgs & cfOgvS) != FALSE)
#endif
           fprintf(OUTFILE,"  dual: T ");
	else
	   fprintf(OUTFILE,"  dual: F ");

        g = g->link;	
	   fprintf(OUTFILE,"\n");
    }
	   fprintf(OUTFILE,"\n");
} 

/***************************************/
/*************************************************/

/*  #################################   */


void Partition3(long *m_ErhU,GLIST2 *gmPgbQv_DlkZ,int cmOvxUli)
{
 
   int k,AoFov_XlVmg,foF_kPrmUvi,moF_wJhg,hoEviT;
   int dhFzhF_xPfmU,bhF_gIivF,tnQ_xIroEivO,tnQ_yPlo;
   GLIST2  *plCzmE_rTl_JmwFc;


   AoFov_XlVmg=1;
   for(k=1;k<hnPabHlgF-cmOvxUli;k++)
     AoFov_XlVmg *=3;

   foF_kPrmUvi = 2*AoFov_XlVmg;
   hoEviT = -foF_kPrmUvi;
   moF_wJhg = -AoFov_XlVmg;

   plCzmE_rTl_JmwFc=gmPgbQv_DlkZ;
   dhFzhF_xPfmU=plCzmE_rTl_JmwFc->mnQgi;
   bhF_gIivF=plCzmE_rTl_JmwFc->noBhhFh;
   /*
   plCzmE_rTl_JmwFc->tnQ_eBofF_oFug=dhFzhF_xPfmU*AoFov_XlVmg;
   plCzmE_rTl_JmwFc->crMw_MvmHgs=bhF_gIivF*AoFov_XlVmg;
   */
   *(m_ErhU +1) = bhF_gIivF*AoFov_XlVmg;
   *(m_ErhU +2) = dhFzhF_xPfmU*AoFov_XlVmg;
   plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;


   while(plCzmE_rTl_JmwFc != NULL)
   {

   tnQ_xIroEivO=plCzmE_rTl_JmwFc->mnQgi;

   switch(tnQ_xIroEivO - dhFzhF_xPfmU){
    
    case -2:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = hoEviT;
	break;

    case -1:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = moF_wJhg;
	break;

    case 0 :
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = 0;
	break;

    case 1:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = AoFov_XlVmg;
	break;

    case 2:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = foF_kPrmUvi;
	break;

    default:
	fprintf(OUTFILE,"\n Error ...\n");
	exit(1);
   }

   tnQ_yPlo=plCzmE_rTl_JmwFc->noBhhFh;
   switch(tnQ_yPlo - bhF_gIivF){
    
    case -2:
	plCzmE_rTl_JmwFc->crMw_Tfn = hoEviT;
	break;

    case -1:
	plCzmE_rTl_JmwFc->crMw_Tfn = moF_wJhg;
	break;

    case 0 :
	plCzmE_rTl_JmwFc->crMw_Tfn = 0;
	break;

    case 1:
	plCzmE_rTl_JmwFc->crMw_Tfn = AoFov_XlVmg;
	break;

    case 2:
	plCzmE_rTl_JmwFc->crMw_Tfn = foF_kPrmUvi;
	break;

    default:
	fprintf(OUTFILE,"\n Error ...\n");
	exit(1);
   }

   dhFzhF_xPfmU=tnQ_xIroEivO;
   bhF_gIivF=tnQ_yPlo;
   plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
  }


   plCzmE_rTl_JmwFc=gmPgbQv_DlkZ;
   tnQ_xIroEivO=plCzmE_rTl_JmwFc->mnQgi;

   switch(tnQ_xIroEivO - dhFzhF_xPfmU){
    
    case -2:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = hoEviT;
	break;

    case -1:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = moF_wJhg;
	break;

    case 0 :
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = 0;
	break;

    case 1:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = AoFov_XlVmg;
	break;

    case 2:
	plCzmE_rTl_JmwFc->tnQ_eBofF_iJtsU = foF_kPrmUvi;
	break;

    default:
	fprintf(OUTFILE,"\n Error ...\n");
	exit(1);
   }

   tnQ_yPlo=plCzmE_rTl_JmwFc->noBhhFh;
   switch(tnQ_yPlo - bhF_gIivF){
    
    case -2:
	plCzmE_rTl_JmwFc->crMw_Tfn = hoEviT;
	break;

    case -1:
	plCzmE_rTl_JmwFc->crMw_Tfn = moF_wJhg;
	break;

    case 0 :
	plCzmE_rTl_JmwFc->crMw_Tfn = 0;
	break;

    case 1:
	plCzmE_rTl_JmwFc->crMw_Tfn = AoFov_XlVmg;
	break;

    case 2:
	plCzmE_rTl_JmwFc->crMw_Tfn = foF_kPrmUvi;
	break;

    default:
	fprintf(OUTFILE,"\n Error ...\n");
	exit(1);
   }

}

/***************************************/
void PartitionFound(long *m_ErhU,GLIST2 *gmPgbQv_DlkZ,int AoFov_XlVmg)
{
 
   int bhF2_Wzo,tnQX;
   GLIST2  *plCzmE_rTl_JmwFc;



   plCzmE_rTl_JmwFc=gmPgbQv_DlkZ;
   bhF2_Wzo=plCzmE_rTl_JmwFc->aoFov_uiFj;
   /*
   plCzmE_rTl_JmwFc->AoFov_UiFj=bhF2_Wzo*AoFov_XlVmg;
   */
   *m_ErhU=bhF2_Wzo*AoFov_XlVmg;
   plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;


   while(plCzmE_rTl_JmwFc != NULL)
   {

   tnQX=plCzmE_rTl_JmwFc->aoFov_uiFj;
   plCzmE_rTl_JmwFc->aoFov_yrU =(tnQX-bhF2_Wzo)*AoFov_XlVmg;
   bhF2_Wzo=tnQX;
   plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
  }


   plCzmE_rTl_JmwFc=gmPgbQv_DlkZ;
   tnQX=plCzmE_rTl_JmwFc->aoFov_uiFj;
   plCzmE_rTl_JmwFc->aoFov_yrU =(tnQX-bhF2_Wzo)*AoFov_XlVmg;


}
/*************************************************/





int  QuicksortDouble(ncU_kBgzMo  *T[],int *rhL)
{
   int k,gmBiiBb2_vcJhgT,plCzmE_yBhv;
   ncU_kBgzMo  *plCzmE_rTl_JmwFc;
   int *dzM_uMzt_hvU,*aiBb_JmwFc;
   int dmF_uPfmEvi,dmF_wPdm_olPk,nn_szQol,nn_szQ;
   double aiBb_JmwFc_2;

    plCzmE_yBhv = 0;

   dzM_uMzt_hvU = (int*)v_alloc(hnPabHlgF,sizeof(int));

   aiBb_JmwFc = (int*)v_alloc(hnPabHlgF,sizeof(int));

   rhL = (int*)v_alloc(hnPabHlgF,sizeof(int));

   if(T==NULL)
   {
      fprintf(stderr,"\n The T matrix is empty.\n");
      exit(1);
   }


   gmBiiBb2_vcJhgT=0;
   plCzmE_yBhv=0;
   for(k=0;k<hnPabHlgF;k++)
   {
      dzM_uMzt_hvU[k]=0;
      aiBb_JmwFc[k]=0;

      rhL[k] = 0;
      plCzmE_rTl_JmwFc = T[k];
      while(plCzmE_rTl_JmwFc != NULL)
      {
           dmF_uPfmEvi=plCzmE_rTl_JmwFc->pvW_kBooFov[MATERNAL];
           dmF_wPdm_olPk=plCzmE_rTl_JmwFc->pvW_kBooFov[PATERNAL];

           nn_szQol=plCzmE_rTl_JmwFc->upOldO_jVzmU[MATERNAL];
           nn_szQ=plCzmE_rTl_JmwFc->upOldO_jVzmU[PATERNAL];

	   if(dmF_uPfmEvi == nn_szQol && dmF_wPdm_olPk == nn_szQ)
	   {
             aiBb_JmwFc[k]++;
	   }
	   else
	     dzM_uMzt_hvU[k]++;
	plCzmE_rTl_JmwFc=plCzmE_rTl_JmwFc->link;
	rhL[k]++;
      }
   }


   plCzmE_yBhv = -1;
   for(k=0;k<hnPabHlgF;k++)
   {
      if(ciS_kCrgT[k] == FALSE)
      {
	if(aiBb_JmwFc[k] == 0)
	{
	  plCzmE_yBhv = k;
        }
      }
   }

  if(plCzmE_yBhv == -1)
  {
   aiBb_JmwFc_2 = 1.1;
   for(k=0;k<hnPabHlgF;k++)
   {
      if(ciS_kCrgT[k] == FALSE)
      {
	if(aiBb_JmwFc[k]/rhL[k] < aiBb_JmwFc_2)
	{
	   aiBb_JmwFc_2 =  aiBb_JmwFc[k]/rhL[k]; 
	  plCzmE_yBhv = k;
        }
      }
   }
  }

  if(plCzmE_yBhv == -1 )
  {
      fprintf(OUTFILE,"\n ERROR: The reduction index is -1. \n");
      exit(1);
  }

/* WARNING: NEED TO CHECK IF NO MARKERS; THEN NO INDEX */
/*
       fprintf(stdout,"\n Index  %d",plCzmE_yBhv);
       fprintf(OUTFILE,"\n Index  %d",plCzmE_yBhv);

   if(ciS_kCrgT[0] == FALSE)
     plCzmE_yBhv = 0;
   else
     plCzmE_yBhv = 1;

   for(k=0;k<hnPabHlgF;k++)
   {
      if(ciS_kCrgT[k] == FALSE)
      {
	if(aiBb_JmwFc[k] == 0)
	{
	  plCzmE_yBhv = k;
        }
      }
   }
*/
   free_check(aiBb_JmwFc);
   free_check(dzM_uMzt_hvU);
   free_check(rhL);

   return(plCzmE_yBhv);
}

