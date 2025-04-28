


#define V_MAIN_VARS
#include "v_prog.h"


/*
#define PARAM_FILE
#define FINAL_OUT
#define DEBUG1
*/



int  main(void)
{
  int j;
  v_boolean v_screen;

 fzH_eFxg = V_LIKE_BY_FAM;

 v_likeli = FALSE;

datafl = fopen(V_PEDD,"r");
if(datafl == NULL)
{
 fprintf(stderr,"\n Could not open datafile.dat. \n");
 exit(1);
}

OUTFILE = NULL;
OUTFILE = fopen(V_FINF,"w");
if(OUTFILE == NULL)
{
 fprintf(stderr,"\n Could not open final.out. \n");
 exit(1);
}

wrHsg2=1;
wrHsg1=1;
veBirBmxF = 3;
 
next_founder_gen2(datafl);
 
if(crMwiFm_Jwh != LINKMAP_PROG && crMwiFm_Jwh != MLINK_PROG )
 {
      fprintf(stderr,"\n Sorry, VITESSE only handles LINKMAP and MLINK");
      fprintf(stderr," runs at the moment. \n The program is %d.\n",crMwiFm_Jwh);
      exit(1);
 }

moMvoF_2 = fopen(V_PEDF,"r");
if(moMvoF_2 == NULL)
{
 fprintf(stderr,"\n Could not open pedfile.dat. \n");
 exit(1);
}
 display_founders2(moMvoF_2);

if(fclose(moMvoF_2) != 0)
{
     fprintf(stderr, "\n File wasn't closed ");
     exit(1);
}
 
screen_out = TRUE;
v_screen = TRUE;

if(fclose(datafl) != 0)
{
     fprintf(stderr, "\n File wasn't closed ");
     exit(1);
}
 
    QuicksortC();

if(fclose(OUTFILE) != 0)
{
     fprintf(stderr, "\n File wasn't closed ");
     exit(1);
  
}

 return(0);
}


