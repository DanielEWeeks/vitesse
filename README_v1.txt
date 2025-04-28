VITESSE Documentation						December 7, 1995 

		VITESSE V1.0  (c)1995  Jeff O'Connell

***************************************************************
Summary:
* Program:       VITESSE 
                 (means 'speed' in French)
* Function:      Likelihood Calculation on Pedigrees
* Version:       1.0
* Date:          1995.12.06
* Author:        Jeff O'Connell 
* Copyright:     (c) 1995 Jeff O'Connell 
* Collaborators: Daniel E. Weeks
* Language:      C (ANSI standard)
* Distribution:  via anonymous FTP from 
		 watson.hgen.pitt.edu directory /pub/vitesse (US)
* Registration:  via e-mail to jeff@watson.hgen.pitt.edu
* Reference:    "The VITESSE algorithm for rapid exact multilocus 
		 linkage analysis via genotype set-recoding and 
		 fuzzy inheritance", O'Connell JR, Weeks DE,
		 Nature Genetics 11:402-408, December 1995

******************************************************************

--------------------------------------------------------------------
DESCRIPTION:

VITESSE is a software package that computes likelihoods with
the functionality of the LINKMAP and MLINK programs from LINKAGE.

VITESSE uses the novel algorithms of set-recoding and fuzzy inheritance
to reduce the number of genotypes needed for exact computation of the
likelihood, which accelerates the calculation. It also represents
multilocus genotypes locus-by-locus to reduce the memory requirements.

The algorithms in VITESSE were developed and coded by Jeff O'Connell at 
the University of Pittsburgh. Dan Weeks at the University of Pittsburgh 
and the Wellcome Trust Centre for Human Genetics at Oxford University 
collaborated.

-------------------------------------------------------------------
FTP DIRECTIONS:

Here are the instructions for retrieving VITESSE.

ftp watson.hgen.pitt.edu
(login as 'anonymous' with your e-mail address as password)

cd pub/vitesse
get vitesse.tar.Z
 
On your machine:

uncompress vitesse.tar.Z
tar xvf vitesse.tar

All files will appear in the directory ./vitesse.

For DOS users: 

cd pub/vitesse/DOS   

get README.DOS

This file will contain further instructions on how to get 
DOS executables when they are available.  
		   
-----------------------------------------------------------------
PLATFORMS:

VITESSE has been tested on the following platforms:

Sun SPARCStation 10, Solaris 2.3 
Sun SPARCStation 2, SunOS 4.1.3
DEC Alpha, OSF/1.3
HP workstation, HP/UX 9.01
IBM SP2 (parallel RS/6000 CPUs), AIX 3.21
Silicon Graphics Challenge, IRIX 5.2
Silicon Graphics Indigo, IRIX 4.0.5
IBM-PC DOS, BorlandC Compiler (being tested)

--------------------------------------------------------------------
COMPILING:

There are a total of 4 executables:

1.) vitesse - the main program 
2.) pedcheck - determines the pedigrees VITESSE can currently 
	       handle. See below.
3.) cnvrt_sh - converts an lcp-generated shell to run VITESSE. See below.
4.) pedshell - converts an lcp-generated shell to run 'pedcheck'. See below. 

If you have the Make utility, then type 'make all'.

**The default compiler is gcc with optimization -O. Edit the
 'Makefile' to change these options.

If you don't use the Make utility then compile as follows:

1.  gcc -O vitsrc.c -o vitesse -lm
2.  gcc -O pedsrc.c -o pedcheck -lm
3.  gcc -O cnvrt_sh.c -o cnvrt_sh -lm
4.  gcc -O pedshell.c -o pedshell -lm

Replace gcc -O as necessary or desired.

**VITESSE is written in ANSI C and uses new style function
prototypes. I've had trouble compiling with older Sun compilers.

**For DOS users, executables are provided.
------------------------------------------------
RUNNING:

Running VITESSE requires that you be familiar with running the 
LINKAGE/FASTLINK package. If you have not used LINKAGE/FASTLINK you 
can find user manuals and other very helpful information on Jurg Ott's
Linkage Analysis Web Site at Columbia University. 
   http://linkage.cpmc.columbia.edu/ 

**VITESSE is currently set up to run in an lcp-generated script file.

**I have included a program called cnvrt_sh to convert any lcp-generated 
  shell file to another shell file for running VITESSE. Basically it 
  strips out calls to the LINKAGE/FASTLINK programs mlink and linkmap 
  and replaces them with calls to vitesse. For example,

 % cnvrt_sh
 Input file :pedin
 Output file: vpedin
   -> 79 lines processed
 %

**Set any necessary paths to where you keep the vitesse and other 
  executables.

**The converted shell file also has some different file names. The output
  files are the same as for LINKAGE/FASTLINK except that the letter 'v' is 
  appended as a prefix. For example, 'final.out' becomes 'vfinal.out'. 

The reason the output is not the same is so that the user can compare
answers with LINKAGE/FASTLINK by running both shell files.

**The LINKAGE/FASTLINK program UNKNOWN is not used by VITESSE and so it is
  stripped from the shell file.  

---------------------------------------------------------------------------
LIMITATIONS:

Version 1.0 will only handle simple pedigrees. Simple means there are no
loops and there is only one set of parents who are founders. I'm actively
working on general pedigrees without loops. Yes, I know this limitation 
is annoying!!

**If VITESSE encounters a pedigree it cannot handle it issues a message and
exits. I have included yet another shell converter 'pedshell' to assist you
in determining which pedigrees are simple. Use 'pedshell' to convert
an lcp-generated shell file to another shell file containing a call
to the program 'pedcheck'. For example, 
 
 %pedshell 
 Input file :pedin
 Output file: ppedin
   -> 79 lines processed
 %

When you run 'ppedin', a message will be displayed after each 
pedigree that is not simple. Delete those pedigrees from the pedigree 
file (make a backup first) and then run VITESSE.

For example, if you generated 'pedin' using lcp, then convert 'pedin'
to say 'ppedin' and run it. Delete pedigrees, if necessary. Then 
convert 'pedin' to 'vpedin' and run it.

------------------------------------------------------------
ALELLE LUMPING/RECODING:

   NEVER RECODE ALLELES. VITESSE does its own allele recoding.
   Hand recoding may lead to errors and any 'allele lumping'
   will not affect VITESSE's running time. My experience is that
   LINKAGE/FASTLINK does not always handle more than 31 alleles 
   at a locus correctly. VITESSE has no restrictions on the number
   of alleles at a locus.

------------------------------------------------------------
SCREEN OUTPUT:

   The final output from VITESSE appears at the end of the run 
   because the program is optimized for MLINK and LINKMAP runs.
   This means that when a trait locus moves between two markers,
   all thetas are done for that pedigree, instead of doing all 
   pedigrees for one theta. VITESSE will print which pedigree is
   is processing.

   The output during the run is much different than LINKAGE.
   VITESSE prints the state of the calculation to give the user
   an idea of the of how complex the problem is. As each nuclear
   family is peeled from the pedigree, VITESSE prints the id's of
   the parents and children, and the number of Parental Pairs. This
   is the number of valid multilocus genotype pairings for the parents.
   For each pair, a calculation involving the compatible child genotypes
   is done, so this number relates to the complexity of that nuclear 
   family, and thus the pedigree - the more pairs, the longer the 
   calculation. As you run different pedigrees or add markers, you 
   should get a feel for the time complexity of the problems.

   When you reach the last nuclear family, the output looks slightly
   different because VITESSE uses a novel algorithm to exploit a
   special symmetry in this family which can lead up to 2-fold 
   speed up. 

**The time and space complexity of LINKAGE/FASTLINK is associated to 
  the product of the number of alleles of the markers, called Maxhap.
  This constant is actually a false indicator of the complexity of the 
  problem.  The space and time complexity of VITESSE is a function of the 
  number of markers and the number of parental genotypes in the pedigree. 
  Note that Maxhap is irrelevant in VITESSE.

I'm developing a preprocessor program that will quickly give this 
complexity information without doing the entire likelihood calculation. 

--------------------------------------------------------------------------
PEDIGREE/MENDELIAN INCONSISTENCIES:

VITESSE assumes to some extent that the pedigree file is correct and does
not have extensive diagnostic checking. Assuming the pedigree file is 
correct, VITESSE is guaranteed to find any Mendelian inheritance 
inconsistencies in your pedigree. If VITESSE finds an inconsistency it 
exits and gives information on the screen and in the file 'vitesse.dbg' 
to assist you in locating the problem.

--------------------------------------------------------------------------
MEMORY:

The memory requirements are a function of the number of loci and 
number of parental pairs. On most problems I've reached
the time complexity before the memory limit. 

--------------------------------------------------------------------------
COMPARING RESULTS FROM VITESSE AND LINKAGE/FASTLINK.

** IMPORTANT: VITESSE ALWAYS computes the exact likelihood.
   In LINKAGE/FASTLINK, if you use UNKNOWN, then the likelihood may 
   not be preserved because if everyone is untyped then they
   are made homozygous 1/1 with that allele frequency. (And I believe
   there are other things LINKAGE does unrelated to UNKNOWN that may
   alter the likelihood.) In general, my experience is that the more
   complex the pedigree, the more likely the likelihood will be 
   preserved. Although these techniques alter the raw likelihood the 
   LOD scores are preserved.

** IMPORTANT: The log base 10 output of LINKAGE/FASTLINK is not
   always accurate. That is because the programs don't use the log10
   function, but rather divide the natural log function by a fixed-
   length constant. This problem may have been fixed in newer versions of
   Pascal LINKAGE 5.* and is fixed LINKAGE 6.0. Note that the inaccuracy
   is minor, becuase the answers are fine to 2-3 decimal places.
   
   For example, on my Sun Workstation I get inaccurate log10 answers from
   FASTLINK.  To get better accuracy in FASTLINK change the log10_ value 
   in commondefs.h to  2.302585093 or replace it with '(log(10.0))'. 

** IMPORTANT: Since VITESSE is a new program and uses completely different
   algorithms than LINKAGE/FASTLINK, comparing output is advised where 
   possible. In general, the bugs I found while developing VITESSE 
   appeared while running 2-point analyses, so if all the 2-point runs
   agreed, other multipoint runs using those markers agreed.

***THUS, TO COMPARE:
     1. CHANGE THE LOG 10 CONSTANT; OTHERWISE LODS WON'T BE ACCURATE
	AND WON'T MATCH WITH VITESSE.
     2. IN THOSE CASES WHEN LINKAGE/FASTLINK GIVES EXACT LIKELIHOODS 
		                 		   ^^^^^
        THEN VITESSE MUST GIVE THE SAME NATURAL LOG AND BASE 10 LOG 
	VALUES. THE LOD SCORES, HOWEVER, MUST ALWAYS MATCH BETWEEN THE
                    ^^^                       ^^^^^^
	PROGRAMS. 

------------------------------------------------------------------
PUBLICATION CITATION: "The VITESSE algorithm for rapid exact multilocus
linkage analysis via genotype set-recoding and fuzzy inheritance",
O'Connell JR, Weeks DE, Nature Genetics 11:402-408, December 1995

Please reference this if you use VITESSE in any published work.

----------------------------------------------------------------------
Contact info

jeff@watson.hgen.pit.edu
