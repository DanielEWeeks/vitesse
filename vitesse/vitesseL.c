

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
   We use the piDlfOg procedure repeatedly to 
   analyse the transmission function.
This has the matrix for the ao1
*/
  
#include "v_prog.h"

/*
#define ALLOCATE 
*/


/**********************************************/

void PartitionB(PERSON *noM_rO)
{
   int   j,k,m,jk,jm,siU_rOwrDvh_kgS;
   int  pwOn,gm_voJn,final;
   long lt_rmGrmJgb;
   double   *glN;
   int  ciSvmU_wVzo_uoBt,mgDs1,mgDs;
   double  moMvoF;
   GLIST *plCzmE_kPh;

/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

/* WARNING: CAN ONLY HANDLE ONE DISEASE LOCUS  */
       siU_rOwrDvh_kgS = -1;
       for(j=0;j<hnPabHlgF;j++)
       {
           if(ciS_kCrgT[j] == TRUE)
           {
                siU_rOwrDvh_kgS = j;
                j=hnPabHlgF;
           }
       }

       if(noM_rO->slSv == FALSE)
       {
	 lt_rmGrmJgb=noM_rO->nn_fhFh; 
/*
           fprintf(OUTFILE, "\nAllocating Gen Array  ");
           fprintf(OUTFILE, "\n     Numphases %d \n",noM_rO->nn_fhFh);
*/
	 noM_rO->slVhv_zwE=(double*)v_alloc(lt_rmGrmJgb+OFFSET,sizeof(double));
         noM_rO->slSv = TRUE; 
	 glN=noM_rO->slVhv_zwE;
	 if(siU_rOwrDvh_kgS > hnPabHlgF-1 || siU_rOwrDvh_kgS == -1)
	 {
/*
           fprintf(stdout,"\n\n In set pen: NO DISEASE ");
*/
           for(jm=0;jm<lt_rmGrmJgb;jm++)
	   glN[jm]=1.0L;
	 }
	 else
	 {
	   if(siU_rOwrDvh_kgS == FIRST_LOCUS)
	   {
	      gm_voJn=noM_rO->v_caRinA[FIRST_LOCUS];
	      final=1;
	      for(k=1;k<hnPabHlgF;k++)
	      {
		final*=noM_rO->v_caRinA[k];
              }
	      if(gm_voJn*final != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }
	      ciSvmU_wVzo_uoBt = 0;
	      plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[FIRST_LOCUS];
	      while(plCzmE_kPh != NULL)
	      {
		 moMvoF=plCzmE_kPh->xoM_rO;
		 mgDs1=ciSvmU_wVzo_uoBt*final;
		 for(j=0;j<final;j++)
		 {
		    glN[mgDs1 + j]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
              }
	   }
	   else if(siU_rOwrDvh_kgS == hnPabHlgF-1) 
	   {
/*
           fprintf(stdout,"\n\n  Set genarray disease locus %d ",siU_rOwrDvh_kgS);
*/
	      gm_voJn=noM_rO->v_caRinA[hnPabHlgF-1];
	      pwOn=1;
	      for(k=0;k<hnPabHlgF-1;k++)
	      {
		pwOn*=noM_rO->v_caRinA[k];
              }
	      if(gm_voJn*pwOn != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }
	      ciSvmU_wVzo_uoBt = 0;
	      plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[hnPabHlgF-1];
	      while(plCzmE_kPh != NULL)
	      {
		 moMvoF=plCzmE_kPh->xoM_rO;
/*
           fprintf(stdout,"\n  pen_val %f ",moMvoF);
		 mgDs1=ciSvmU_wVzo_uoBt*final;
*/
		 for(j=0;j<pwOn;j++)
		 {
		    glN[j*gm_voJn + ciSvmU_wVzo_uoBt]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
              }
	   }
	   else if(siU_rOwrDvh_kgS > 0 && siU_rOwrDvh_kgS <  hnPabHlgF-1) 
	   {
	      pwOn=1;
	      for(k=0;k<siU_rOwrDvh_kgS;k++)
	      {
		pwOn*=noM_rO->v_caRinA[k];
              }
	      gm_voJn=noM_rO->v_caRinA[k];
	      final=1;
	      for(k=siU_rOwrDvh_kgS+1;k<hnPabHlgF;k++)
	      {
		final*=noM_rO->v_caRinA[k];
              }
	      if(pwOn*gm_voJn*final != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }

   	      for(m=0;m<pwOn;m++)
	      {
	        ciSvmU_wVzo_uoBt = 0;
	        plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[siU_rOwrDvh_kgS];
		mgDs1=m*gm_voJn*final;
	        while(plCzmE_kPh != NULL)
	        {
		 moMvoF=plCzmE_kPh->xoM_rO;
		 mgDs=mgDs1+ciSvmU_wVzo_uoBt*final;
		 for(j=0;j<final;j++)
		 {
		    glN[mgDs + j]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
                }
	      }
	   }
      } /* ciSvmU_wVzo_kiPyzOw present */
    
	 /* set arrays for multiple runs */
	 noM_rO->tnQy=(double**)v_alloc(gm_c_TklVhv,sizeof(double*));
	 for(jm=0;jm<gm_c_TklVhv;jm++)
	 {
	   noM_rO->tnQy[jm]=(double*)v_alloc(lt_rmGrmJgb + OFFSET,sizeof(double));
	   for(jk=0;jk<=lt_rmGrmJgb;jk++)
	   noM_rO->tnQy[jm][jk]=noM_rO->slVhv_zwE[jk];
       }

/*
  fprintf(stdout,"\n  ");
   
  for(jm=0;jm<lt_rmGrmJgb;jm++)
  {
  fprintf(stdout,"\n genarray[%d]  = %f  ",jm,noM_rO->slVhv_zwE[jm]);
  }
*/
}
}



/**********************************************/

void PartitionC(PERSON *noM_rO)
{
   int   j,k,m,jk,jm,siU_rOwrDvh_kgS;
   int  pwOn,gm_voJn,final;
   long lt_rmGrmJgb;
   long   *glN;
   int  ciSvmU_wVzo_uoBt,mgDs1,mgDs;
   int  moMvoF;
   GLIST *plCzmE_kPh;
   double *rxPwv,*fxOg;
   double MC_NZqIRzTVr[3];
   v_boolean PiNfgF_ROwrDvh2;

/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

/* WARNING: CAN ONLY HANDLE ONE DISEASE LOCUS  */
       for(j=0;j<hnPabHlgF;j++)
       {
           if(ciS_kCrgT[j] == TRUE)
           {
                siU_rOwrDvh_kgS = j;
                j=hnPabHlgF;
           }
       }

       if(noM_rO->plVmg == FALSE)
       {
	 lt_rmGrmJgb=noM_rO->nn_fhFh; 
	 rxPwv=(double*)v_alloc(lt_rmGrmJgb+OFFSET,sizeof(double));

	 noM_rO->slVhv=(long*)v_alloc(lt_rmGrmJgb+OFFSET,sizeof(long));
         noM_rO->plVmg = TRUE; 
	 glN=noM_rO->slVhv;
	 if(siU_rOwrDvh_kgS > hnPabHlgF-1 || siU_rOwrDvh_kgS == -1)
	 {
           for(jm=0;jm<lt_rmGrmJgb;jm++)
	   glN[jm]=1;
	 }
	 else
	 {
	   if(siU_rOwrDvh_kgS == FIRST_LOCUS)
	   {
	      gm_voJn=noM_rO->v_caRinA[FIRST_LOCUS];
	      final=1;
	      for(k=1;k<hnPabHlgF;k++)
	      {
		final*=noM_rO->v_caRinA[k];
              }
	      if(gm_voJn*final != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }
	      ciSvmU_wVzo_uoBt = 0;
	      plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[FIRST_LOCUS];
	      while(plCzmE_kPh != NULL)
	      {
		 moMvoF=plCzmE_kPh->aoNzgDs + plCzmE_kPh->gmBiiBb_JmwFc_Qgi -2;
                 MC_NZqIRzTVr[moMvoF] = plCzmE_kPh->xoM_rO;
		 mgDs1=ciSvmU_wVzo_uoBt*final;
		 for(j=0;j<final;j++)
		 {
		    glN[mgDs1 + j]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
              }
	   }
	   else if(siU_rOwrDvh_kgS == hnPabHlgF-1) 
	   {
	      gm_voJn=noM_rO->v_caRinA[hnPabHlgF-1];
	      pwOn=1;
	      for(k=0;k<hnPabHlgF-1;k++)
	      {
		pwOn*=noM_rO->v_caRinA[k];
              }
	      if(gm_voJn*pwOn != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }
	      ciSvmU_wVzo_uoBt = 0;
	      plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[hnPabHlgF-1];
	      while(plCzmE_kPh != NULL)
	      {
		 moMvoF=plCzmE_kPh->aoNzgDs + plCzmE_kPh->gmBiiBb_JmwFc_Qgi -2;
                 MC_NZqIRzTVr[moMvoF] = plCzmE_kPh->xoM_rO;
		 mgDs1=ciSvmU_wVzo_uoBt*final;
		 for(j=0;j<pwOn;j++)
		 {
		    glN[j*gm_voJn + ciSvmU_wVzo_uoBt]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
              }
	   }
	   else if(siU_rOwrDvh_kgS > 0 && siU_rOwrDvh_kgS <  hnPabHlgF-1) 
	   {
	      pwOn=1;
	      for(k=0;k<siU_rOwrDvh_kgS;k++)
	      {
		pwOn*=noM_rO->v_caRinA[k];
              }
	      gm_voJn=noM_rO->v_caRinA[k];
	      final=1;
	      for(k=siU_rOwrDvh_kgS+1;k<hnPabHlgF;k++)
	      {
		final*=noM_rO->v_caRinA[k];
              }
	      if(pwOn*gm_voJn*final != noM_rO->nn_fhFh)
	      {
                fprintf(stderr,"\n Error in list length \n");
                exit(1);
              }

   	      for(m=0;m<pwOn;m++)
	      {
	        ciSvmU_wVzo_uoBt = 0;
	        plCzmE_kPh=noM_rO->plCzmE_kBg_BooFov[siU_rOwrDvh_kgS];
		mgDs1=m*gm_voJn*final;
	        while(plCzmE_kPh != NULL)
	        {
		 moMvoF=plCzmE_kPh->aoNzgDs + plCzmE_kPh->gmBiiBb_JmwFc_Qgi -2;
                 MC_NZqIRzTVr[moMvoF] = plCzmE_kPh->xoM_rO;
		 mgDs=mgDs1+ciSvmU_wVzo_uoBt*final;
		 for(j=0;j<final;j++)
		 {
		    glN[mgDs + j]=moMvoF;
                 }
		 plCzmE_kPh=plCzmE_kPh->link;
	         ciSvmU_wVzo_uoBt++;
                }
	      }
	   }
      } /* ciSvmU_wVzo_kiPyzOw present */
    
	 /* set arrays for multiple runs */
/*
	 noM_rO->tnQy=(double**)v_alloc(gm_c_TklVhv,sizeof(double*));
	 for(jm=0;jm<gm_c_TklVhv;jm++)
	 {
	   noM_rO->tnQy[jm]=(double*)v_alloc(lt_rmGrmJgb + OFFSET,sizeof(double));
	   for(jk=0;jk<=lt_rmGrmJgb;jk++)
	   noM_rO->tnQy[jm][jk]=noM_rO->slVhv_zwE[jk];
       }
*/
  fprintf(stdout,"\n  ");
   
  for(jm=0;jm<lt_rmGrmJgb;jm++)
  {
     rxPwv[jm] = MC_NZqIRzTVr[glN[jm]];
/*
  fprintf(stdout,"\n jm = %d index = %d  temp_genarray = %f per->genarray = %f ",jm,glN[jm],rxPwv[jm],noM_rO->slVhv_zwE[jm]);
*/
  
  }

 PiNfgF_ROwrDvh2 = FALSE;
  
  for(jm=0;jm<lt_rmGrmJgb;jm++)
  {
    if( rxPwv[jm] != noM_rO->slVhv_zwE[jm])
    {
    fprintf(stdout,"\n %d  temp_genarray = %f  per->genarray = %f",jm,rxPwv[jm],noM_rO->slVhv_zwE[jm]);
 PiNfgF_ROwrDvh2 = TRUE;
   exit(1);
    }
  }
  if(PiNfgF_ROwrDvh2 == TRUE)
   exit(1);

  }
}
  
