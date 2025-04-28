

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
   We use the piDlfOg procedure repeatedly to 
   analyse the transmission function.
*/
  
#include "v_prog.h"

/*
#define CHECKSORT
#define CHECKSORT3
#define GSORT
*/

/*********************************************************/


void Add_Alist(char **A,ncU_kBgzMo **B,int p,int r,int ciSvmU_kPh)
{

 int q;

 if(p<r)
 {
   q=set_genotypes(A,B,p,r,ciSvmU_kPh);
   Add_Alist(A,B,p,q,ciSvmU_kPh);
   Add_Alist(A,B,q+1,r,ciSvmU_kPh);
 }
}


int  set_genotypes(char **A,ncU_kBgzMo **B,int p, int r,int ciSvmU_kPh)
{
    int i,j;
    ncU_kBgzMo  *vxUliT;
    char *x,*lxVh;
    v_boolean PiNfgF_ROwrDvh;

    x=A[p];
    i=p-1;
    j=r+1;

    PiNfgF_ROwrDvh = FALSE;
    while(!PiNfgF_ROwrDvh){
      do{
	j--;
      } while(Quicksort6(A[j],x,ciSvmU_kPh));
      
      do{
	i++;
      } while(Quicksort6(x,A[i],ciSvmU_kPh));
    
     if(i<j)
     {
      lxVh = A[i];
      A[i]=A[j];
      A[j] = lxVh;

      vxUliT = B[i];
      B[i]=B[j];
      B[j]=vxUliT;
    }
    else
     PiNfgF_ROwrDvh = TRUE;
    }
      return j;
}
    
v_boolean Quicksort6(char *v1, char *v2, int ciSvmU_kPh)
{
  int j;

/*
     fprintf(stdout,"\n Entering compare. \n");
  
  fprintf(stdout,"\n V1");
  for(j=0;j<ciSvmU_kPh;j++)
  {
     fprintf(stdout," %d ",v1[j]);
  }

  fprintf(stdout,"\n V2");
  for(j=0;j<ciSvmU_kPh;j++)
  {
     fprintf(stdout," %d ",v2[j]);
  }
*/

  /* is v1 > v2 ? */


  for(j=0;j<ciSvmU_kPh;j++)
  {
    if(v1[j] > v2[j])
    {
      return(TRUE);
    }
    else if(v2[j] > v1[j])
    {
      return(FALSE);
    }
  }

  if(j==ciSvmU_kPh)
  {
    return(FALSE);
  }
  else
  {
     fprintf(stderr,"\n Error in compare. \n");
     exit(1);
  }
    return(FALSE);
}

/* ***************************************  */
    
v_boolean assign_rel(char *v1, char *v2, int ciSvmU_kPh)
{
  int j;

     /*
     fprintf(stdout,"\n Entering equal. \n");
  
  fprintf(stdout,"\n V1");
  for(j=0;j<ciSvmU_kPh;j++)
  {
     fprintf(stdout," %d ",v1[j]);
  }

  fprintf(stdout,"\n V2");
  for(j=0;j<ciSvmU_kPh;j++)
  {
     fprintf(stdout," %d ",v2[j]);
  }
   */

  /* is v1 = v2 ? */


  for(j=0;j<ciSvmU_kPh;j++)
  {
    if(v1[j] != v2[j])
    {
      return(FALSE);
    }
  }

  if(j==ciSvmU_kPh)
  {
    return(TRUE);
  }
  else
  {
     fprintf(stderr,"\n Error in compare. \n");
     exit(1);
  }
    return(TRUE);
}



void setupptrs_nucfam(ncU_kBgzMo *plCzmE_rTl_JmwFc,NUC_FAM *vjVzmU,int cmOvxUli)
{
    
    int k,plCzmE_hUziU,fnJob_hrAv;
    PERSON *tgBo_JmwJxvT;
    GLIST  *A,*B,*C,*D;

    fnJob_hrAv=vjVzmU->fnJob_hrAv;



 for(k=0;k< fnJob_hrAv;k++)
 {
      tgBo_JmwJxvT=vjVzmU->ahXvi[k];

      A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
      plCzmE_hUziU = tgBo_JmwJxvT->v_caRinA[cmOvxUli];
      if(plCzmE_hUziU == 0)
      {
       fprintf(OUTFILE,"\n glength is 0");
       exit(1);
      }
	
      

     


	switch(plCzmE_hUziU){

	   case 1:
              break;

           case 2:
             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     if(A->noBhhFh ==  B->noBhhFh)
	     {
	       if(A->mnQgi >  B->mnQgi)
	       {
                 /* B A */
                  A->link = NULL;
		  B->link=A;
		  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
	       }
             }
	     else if(A->noBhhFh >  B->noBhhFh)
	     {
                 /* B A */
                  A->link = NULL;
		  B->link=A;
		  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
             }
             
             break;

           case 3:
             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     C=B->link;
	     switch(A->noBhhFh + B->noBhhFh + C->noBhhFh)
	     {
		/* sum can only be 2 or 4. If 2 then we have 0/0 0/2 2/? */
		case 2:
	        if(A->noBhhFh == 2)
	        {
                  /* means 2/? is first slot */

	         if(B->mnQgi ==2)
                 {
		   /* means 0/2 in 2nd  and 0/0 in 3rd*/
		  /* C B A */ 
	          A->link=NULL;
	          B->link=A;
	          C->link=B;
                  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=C;
		 }
		 else
		 {
	       
		   /* means 0/0 in 2nd  and 0/2 in 3rd*/
		  /* B C A */ 
	          A->link=NULL;
	          B->link=C;
		  C->link=A;
                  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
                 }
                }
                else
	        {
		  /* 2/? is 2nd or 3rd */

	          if(A->mnQgi == 2)
	          {
		    /* means 0/2 is 1st */
                    if(B->noBhhFh == 0)
		    {
		      /* means 0/0 is 2nd */
		       /* B A C */
     			A->link=C;
			B->link=A;
                        tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
                    }
		    else
		    {
		      /* means 0/0 is 3rd */
		       /* C A B */
     			A->link=B;
			B->link=NULL;
			C->link=A;
                        tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=C;
                    }
                 }
		 else
		 {
		    /* means 0/0 is 1st */
                    if(B->noBhhFh == 2)
		    {
		      /* means 2/? is 2nd */
		       /* A C B*/
     			A->link=C;
			B->link=NULL;
			C->link=B;
                    }
		    /* else in order */
                 }
               }
		break;

		/* now we have 2/0 2/2 0/?  */
		case 4:
	        if(C->noBhhFh == 0)
	        {
                    /* means 0/? is 3rd */
	         if(B->mnQgi ==0)
                 {
		  /* means 2/0 is 2nd */
		  /* C B A */ 
	          A->link=NULL;
	          B->link=A;
	          C->link=B;
                  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=C;
		 }
		 else
		 {
	          /* means 2/2 is 2nd */ 
		  /* C A B */ 
	          A->link=B;
	          B->link=NULL;
		  C->link=A;
                  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=C;
                 }
                }
                else
	        {
		  /* means 0/? is 1st or 2nd */
	          if(C->mnQgi == 0)
	          {
		     /* means 2/0 is 3rd */
                    if(B->noBhhFh == 2)
		    {
		       /* means 2/2 is second */
		       /* A C B */
     			A->link=C;
			B->link=NULL;
			C->link=B;
                    }
		    else
		    {
		       /* means 2/2 is 1st */
		       /* B C A */
     			A->link=NULL;
			B->link=C;
			C->link=A;
                        tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
                    }
                  }
		  else
		  {
		    /* means that 2/2 is 3rd */
                    if(B->noBhhFh == 0)
		    {
		       /* B A C*/
     			A->link=C;
			B->link=A;
                        tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
                    }
                 }
               }
               break;
            }

	     /* check the sort */
             break;


	   /* Will end up with 0/0 0/2 2/0 2/2  */
           case 4:
             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     C=B->link;
	     D=C->link;
	     /* put the two with noBhhFh = 0 in the first two positions */
	     if(A->noBhhFh == 0)
	     {

	       if(C->noBhhFh == 0)
	       {
		  /* A C B D */ 
	          A->link=C;
	          B->link=D;
	          C->link=B;
               }
               else if(D->noBhhFh == 0)
	       {
		  /* A D C B*/ 
	          A->link=D;
	          B->link=NULL;
	          C->link=B;
	          D->link=C;
               }
            }
            else
	    {
	       /* 2/? in 1st place */
	       if(B->noBhhFh == 2)
	       {
	           /* 2/? in 2nd place */
		  /* D C B A */ 
	          A->link=NULL;
	          B->link=A;
	          C->link=B;
	          D->link=C;
                  tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=D;

               }
	       else
	       {
	           /* 0/? in 2nd place */
	          if(C->noBhhFh == 0)
                  {
		     /* C B A D */ 
	             A->link=D;
	             B->link=A;
	             C->link=B;
                     tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=C;
		  }
	          else
	          {
		    /* D B C A */ 
	             A->link=NULL;
	             B->link=C;
	             C->link=A;
	             D->link=B;
                     tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=D;
	          }

               }
             }

              /* Now sort the 1-2 and 3-4 fields separately */
             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     C=B->link;
	     D=C->link;
	     if(A->mnQgi == 2)
	     {
	           /* B A C D */ 
	            A->link=C;
	            B->link=A;
                    tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli]=B;
             }
	    
             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     C=B->link;
	     D=C->link;
	     if(C->mnQgi == 2)
	     {
	  	  /* A B D C */ 
	           B->link=D;
	           C->link=NULL;
	           D->link=C;
            }

	     /* check the sort */
             break;

	     default:
               fprintf(OUTFILE,"\n Error");
	       exit(1);

       }/* switch  */

	
/* Assign MC_RMsH to each sorted configuration to use in the slVhv_kiJli sort */	
	switch(plCzmE_hUziU){

	     case 1:

             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     switch(A->noBhhFh)
	     {
		case 0:
		switch(A->mnQgi)
		{
		  case 0:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=0;
		    break;

		  case 1:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=1;
		    break;

		  case 2:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=2;
		    break;
		   
                  default:
                   fprintf(stderr,"\n Unexpected Recombination \n");
                   fprintf(stderr,"\n (%d,%d)\n",A->noBhhFh,A->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	           exit(1);
                }
                  break;

		case 1:
		switch(A->mnQgi)
		{
		  case 0:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=3;
		    break;

		  case 1:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=4;
		    break;

		  case 2:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=5;
		    break;
		  
                  default:
        first_bit_array(A);
	display_founders3(A);
                   fprintf(stderr,"\n Unexpected Recombination \n");
                   fprintf(stderr,"\n (%d,%d)\n",A->noBhhFh,A->mnQgi);
	           exit(1);
                }
                  break;

		case 2:
		switch(A->mnQgi)
		{
		  case 0:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=6;
		    break;

		  case 1:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=7;
		    break;

		  case 2:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=8;
		    break;
		  
                  default:
                   fprintf(stderr,"\n Unexpected Recombination \n");
                   fprintf(stderr,"\n (%d,%d)\n",A->noBhhFh,A->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	           exit(1);
                }
                  break;

	        default:
                fprintf(stderr,"\n Unexpected Recombination \n");
                fprintf(stderr,"\n (%d,%d)\n",A->noBhhFh,A->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	        exit(1);
              }
	       break;

           /* There are two mg  */
	     case 2:

             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     switch(A->noBhhFh)
	     {
		/* 0/?  ?/? */
		case 0:
		switch(A->mnQgi)
		{
		  /* 0/0  ?/? */
		  case 0:
		    switch(B->noBhhFh)
		    {
		      /* 0/0  0/? */
		      case 0:
		       switch(B->mnQgi)
		       {
		        /* 0/0  0/2 */
			case 2:
	                plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=9;
		        break;

	                default:
                        fprintf(stderr,"\n Unexpected Recombination \n");
                        fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		        fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	                exit(1);
                       }
		       break;
		      /* 0/0  2/? */
		      case 2:
		       switch(B->mnQgi)
		       {
		        /* 0/0  2/0 */
			case 0:
	                plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=10;
		        break;

		        /* 0/0  2/2 */
			case 2:
	                plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=11;
		        break;

	                default:
                        fprintf(stderr,"\n Unexpected Recombination \n");
                        fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		        fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	                exit(1);
                       }
		       break;

	              default:
                      fprintf(stderr,"\n Unexpected Recombination \n");
                      fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		      fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	              exit(1);
                   }
		    break;

		  /* 0/1  ?/? */
		  case 1:
                   switch(B->noBhhFh)
		   {
		    /* 0/1  2/? */
		    case 2:
		    switch(B->mnQgi)
		    {
		     /* 0/1  2/1 */
		     case 1:
		      plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=12;
		      break;

	              default:
                      fprintf(stderr,"\n Unexpected Recombination \n");
		      fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		      fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	              exit(1);
		    }
		    break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);
		  }
                  break;

		/* 0/2  ?/? */
		    case 2:
		    switch(B->mnQgi)
		    {
		     /* 0/2  2/0 */
		     case 0:
		      plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=13;
		      break;

		     /* 0/2  2/2 */
		     case 2:
		      plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=14;
		      break;

	              default:
                      fprintf(stderr,"\n Unexpected Recombination \n");
                      fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		      fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	              exit(1);
		    }
		    break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);

                  }
		  break;

   /* back to first efBo_QzrSh  */

		/* 1/?  ?/? */
		case 1:
                   switch(A->mnQgi)
		   {
		    /* 1/0  ?/? */
		    case 0:
		    switch(B->noBhhFh)
		    {
		     /* 1/0  1/2 */
		     case 1:
		      switch(B->mnQgi)
		      {
			case 2:
		        plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=15;
		        break;

	                default:
                        fprintf(stderr,"\n Unexpected Recombination \n");
                        fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		        fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	                exit(1);
		      }
		      break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);
                   }
		    break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);
                  }
		  break;

	        /* 2/?  ?/? */
		case 2:
                   switch(A->mnQgi)
		   {
		    /* 2/0  ?/? */
		    case 0:
		    switch(B->noBhhFh)
		    {
		     /* 2/0  2/? */
		     case 2:
		      switch(B->mnQgi)
		      {
			case 2:
		        plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=16;
		        break;

	                default:
                        fprintf(stderr,"\n Unexpected Recombination \n");
                        fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		        fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	                exit(1);
		      }
		      break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);
                   }
		    break;

	            default:
                    fprintf(stderr,"\n Unexpected Recombination \n");
                    fprintf(stderr,"\n (%d,%d)",A->noBhhFh,A->mnQgi);
		    fprintf(stderr," (%d,%d)",B->noBhhFh,B->mnQgi);
        first_bit_array(A);
	display_founders3(A);
	            exit(1);
                  }
		  break;
		
	        default:
                    fprintf(OUTFILE,"\n Error");
	            exit(1);
             }
	     break;

	     /* There are three mg  */
             case 3:

             A=tgBo_JmwJxvT->plCzmE_kBg_BooFov[cmOvxUli];
	     B=A->link;
	     C=B->link;
	     switch(A->mnQgi)
	     {
		/* 0/0 ? ? */
		case 0:
		switch(B->noBhhFh)
		{
		  /* 0/0 0/2 ? */
		  case 0:
		  switch(C->mnQgi)
		  {
		     /* 0/0 0/2  2/0 */
		     case 0:
		       plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=17;
		       break;

		     /* 0/0 0/2  2/2 */
		     case 2:
		       plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=18;
		       break;

	             default:
                       fprintf(OUTFILE,"\n Error");
	               exit(1);
		  }
		  break;

		  /* 0/0 2/0 2/2 */
		  case 2:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=19;
		    break;

	          default:
                    fprintf(OUTFILE,"\n Error");
	            exit(1);
		}
                 break;

		/* 0/2 2/0 2/2  */
		case 2:
		    plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=20;
		    break;

	       default:
                  fprintf(OUTFILE,"\n Error");
	          exit(1);
	     }
	      break;

	     case 4:
	       plCzmE_rTl_JmwFc->lt10_UlgBo_Tfn[k]=21;
        }
	
	
     }/* num tgBo_JmwJxvT */
}


