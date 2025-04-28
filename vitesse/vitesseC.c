

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
#define PARAM_FILE
#define FINAL_OUT
#define DEBUG1
*/


void QuicksortC(void)
{
 int j,k,m,tnQ_hQlfTv_Jhl;
 long nnF;
 v_boolean v_screen;

v_screen = TRUE;
if(crMwiFm_Jwh == ILINK_PROG)
  v_screen = FALSE;
if(v_screen)
   fprintf(stdout," \n Recombination fractions for each iteration \n");
   fprintf(OUTFILE," \n Recombination fractions for each iteration \n");
   for(m=0;m<gm_c_TklVhv;m++)
   {

if(v_screen)
           fprintf(stdout," \n %d ",m);
           fprintf(OUTFILE," \n %d ",m);
     for(k=0;k<hnPabHlgF-1;k++)
     {
if(v_screen)
           fprintf(stdout," %f ",tnQ_kSlyBmw_rhPabHlgF_kUi[m][k]);
           fprintf(OUTFILE," %f ",tnQ_kSlyBmw_rhPabHlgF_kUi[m][k]);
     }

   }
if(v_screen)
           fprintf(stdout,"\n");
           fprintf(OUTFILE,"\n");


if(v_screen)
    fprintf(stdout,"\n Markers:              ");
    fprintf(OUTFILE,"\n Markers:              ");
    for(k=0;k<hnPabHlgF;k++)
    {
if(v_screen)
       fprintf(stdout," %2d ",k+ 1 );
       fprintf(OUTFILE," %2d ",k+1);
    }

if(v_screen)
    fprintf(stdout,"\n Position on map:      ");
    fprintf(OUTFILE,"\n Position on map:      ");
    for(k=0;k<hnPabHlgF;k++)
    {
if(v_screen)
       fprintf(stdout," %2d ",jmL5[k]+ 1 );
       fprintf(OUTFILE," %2d ",jmL5[k]+1);
    }

if(v_screen)
    fprintf(stdout,"\n Number of alleles:    ");
    fprintf(OUTFILE,"\n Number of alleles:    ");
   nnF=1;
   for(k=0;k<hnPabHlgF;k++)
   {
    nnF *=fzH_wVzo_xoBhh[k]->tnQ_hQlfTv_Jhl; 
if(v_screen)
    fprintf(stdout," %2d ",fzH_wVzo_xoBhh[k]->tnQ_hQlfTv_Jhl);
    fprintf(OUTFILE," %2d ",fzH_wVzo_xoBhh[k]->tnQ_hQlfTv_Jhl);
   }

if(v_screen)
    fprintf(stdout,"\n Number of haplotypes:%ld ",nnF);
    fprintf(OUTFILE,"\n Number of haplotypes: %ld ",nnF);

pdFic2=(int) (gmBiiBb/(sizeof(V_INT)*BYTE)) + 1;
/*
fprintf(stdout,"\n MAX_INTS %d max_int %d \n",pdFic2,gmBiiBb);
fprintf(OUTFILE,"\n MAX_INTS %d max_int %d \n",pdFic2,gmBiiBb);
*/
 
/*  do_streamfile(); */

moMvoF_2 = fopen(V_PEDF,"r");
if(moMvoF_2 == NULL)
{
 fprintf(stderr,"\n Could not open pedfile.dat. \n");
 exit(1);
}
/*
 display_founders2(moMvoF_2);
*/

 mc_ezM=(double ****)v_alloc(NUM_GEN,sizeof(double ***));
  for(j=0;j<NUM_GEN;j++)
     mc_ezM[j]=NULL;

  PfQ=(double ***)v_alloc(NUM_GEN,sizeof(double**));

 pdFi = (int*) (v_alloc(mn_nzU,sizeof(int)));

 ncU_kSlyBmw = (double **) (v_alloc(mn_nzU,sizeof(double*)));

 for(j=0;j<mn_nzU;j++)
 {
   ncU_kSlyBmw[j] = (double *) (v_alloc(gm_c_TklVhv,sizeof(double)));
 }

 ncU_nMrhU = 3;
 for(k=0;k<hnPabHlgF;k++)
   ncU_nMrhU *=veBirBmxF;

if(v_likeli)
{
 processchild(tnQ_kSlyBmw_rhPabHlgF_kUi,MALE);
 if(rtIgzMo > 0)
   processchild(tnQ_h_rhP,FEMALE);
  
 /*
 print_nuclear_fam();
 */
if(v_screen)
  fprintf(stdout,"\n\nGenerating Isozygote Class Matrices \n");
 print_founder_pairs(MALE);
 if(rtIgzMo > 0)
   print_founder_pairs(FEMALE);
}

if(v_screen)
   fprintf(stdout,"\n Number of Pedigrees: %d",mn_nzU);


 nd_ivD=0;
 fzNv =0.0;
 for(k=0;k<mn_nzU;k++)
 {

if(screen_out)
{
if(v_screen)
{
   fprintf(stdout,"\n\n ################################################# ");
   fprintf(stdout,"\n No. %d/%d  Processing pedigree %3d",k+1,mn_nzU,poMvoF[k]);
   fprintf(stdout,"\n The pedigree has %4d persons",poMvoF_2[k]);
}
}


 pdFi2 = 0;
 SnNvgSrx_KzJi=poMvoF_2[k];

 next_genotype(moMvoF_2);
 /*
 constructlist();
 */


 free_found_pairs();
 /*
 recomb_class();
 */
 debugGlist();

 fnBov_izUrl=Partition5(k);


if(v_likeli)
 print_bit_vect();

 free_found(k); 
 

TzOh_BiiBb = 0;
tzOhkPhv = 0;

if(v_likeli)
{
 get_locus_freq(ncU_kPhrUrlO,ncU_kPhrUrlO->fiTg_Qvw);

/*
fprintf(OUTFILE,"\nTotal number of Parental Pairs in Top Founder: %ld \n",tzOhkPhv);
*/
/*
if(v_screen)
fprintf(stdout,"\nTotal number of Parental Pairs in Top Founders: %ld \n",tzOhkPhv);
*/

/*
fprintf(OUTFILE,"\nTotal number of Parental Pairs in Pedigree %d =  %ld ",poMvoF[k],TzOh_BiiBb);
*/
if(v_screen)
fprintf(stdout,"\nTotal number of Parental Pairs in the Pedigree:  %ld ",TzOh_BiiBb);


 nd_ivD+=pdFi2;
 pdFi[k]=pdFi2;

 for(j=0;j<gm_c_TklVhv;j++)
 {
  ncU_kSlyBmw[k][j]=mkGfm[j];
  fzNv +=mkGfm[j];

 }
}

  Peel_Graph();

  Peel_Graph_Down();

  PrintALIST();

  openfile(fnBov_izUrl);

if(v_likeli)
  addlist();

 }
  fprintf(OUTFILE,"\n");
if(v_likeli)
{


weight_heuristic(OUTFILE);
sort_glist5(OUTFILE);
v_citation(OUTFILE);

 if(crMwiFm_Jwh != ILINK_PROG) 
 {
   sort_glist5(stdout);
 }

 if(crMwiFm_Jwh == ILINK_PROG) 
 {
   if(ncU_kBooFov == NULL)
   {
     ncU_kBooFov = (double *)v_alloc(mn_nzU,sizeof(double));
     ncU_lOv = (double *)v_alloc(mn_nzU,sizeof(double));
   }
   for(j=0;j<mn_nzU;j++)
   {
      ncU_kBooFov[j] = log(ncU_kSlyBmw[j][0])-pdFi[j]*log((double)ncU_nMrhU);
      ncU_lOv[j] = log10(ncU_kSlyBmw[j][0])-pdFi[j]*log10((double)ncU_nMrhU);
   }
 }

 for(j=0;j<mn_nzU;j++)
 {
   free_check(ncU_kSlyBmw[j]);
 }
   free_check(ncU_kSlyBmw);
 
  free_check(pdFi);

  do_lspfile();
}
  /*
  free_check(poMvoF_2);
  free_check(poMvoF);
  */
if(fclose(moMvoF_2) != 0)
{
     fprintf(stderr, "\n File wasn't closed ");
     exit(1);
}
} 

