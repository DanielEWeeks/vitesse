#
#	Makefile for VITESSE 
#



CC =  gcc 

HFILE =  v_prog.h v_var.h v_extern.h	
CFLAGS = -O 
GCCFLAG =  -c 


MODG  = main.o \
	vitesseA.o vitesseB.o  vitesseC.o \
	vitesseD.o vitesseE.o  vitesseF.o \
	vitesseG.o vitesseH.o  vitesseI.o \
	vitesseJ.o vitesseK.o  vitesseL.o \
	vitesseM.o vitesseN.o  vitesseO.o \
	vitesseP.o vitesseQ.o  vitesseR.o \
	vitesseS.o vitesseT.o  vitesseU.o \
	vitesseV.o vitesseW.o  vitesseX.o \
	vitesseY.o vitesseZ.o  vitesse1.o \
	vitesse2.o 
#vitesse3.o   

MODH  = main2.o \
	vitesseA.o vitesseB.o  vitesseC.o \
	vitesseD.o vitesseE.o  vitesseF.o \
	vitesseG.o vitesseH.o  vitesseI.o \
	vitesseJ.o vitesseK.o  vitesseL.o \
	vitesseM.o vitesseN.o  vitesseO.o \
	vitesseP.o vitesseQ.o  vitesseR.o \
	vitesseS.o vitesseT.o  vitesseU.o \
	vitesseV.o vitesseW.o  vitesseX.o \
	vitesseY.o vitesseZ.o  vitesse1.o \
	vitesse2.o    

MODE  = cnvrt_sh.o

all	:	vitesse pedcheck cnvrt_sh pedshell

vitesse	:	$(HFILE) $(MODG)
	$(CC) $(CFLAGS) -o vitesse $(MODG) -lm  

pedcheck	:	$(HFILE) $(MODH)
	$(CC) $(CFLAGS) -o pedcheck $(MODH) -lm  

cnvrt_sh	:	$(HFILE) $(MODE)
	$(CC) $(CFLAGS) -o cnvrt_sh $(MODE) -lm 

cnvrt_sh.o		: cnvrt_sh.c	
	$(CC) $(CFLAGS) $(GCCFLAG) -o cnvrt_sh.o cnvrt_sh.c 

pedshell	:	pedshell.o	
	$(CC) $(CFLAGS) -o pedshell  pedshell.o -lm 


main.o		:	$(HFILE) main.c	
	$(CC) $(CFLAGS) $(GCCFLAG) -o main.o main.c 

main2.o		:	$(HFILE) main2.c	
	$(CC) $(CFLAGS) $(GCCFLAG) -o main2.o main2.c 

vitesseA.o        :       $(HFILE) vitesseA.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseA.o vitesseA.c 

vitesseB.o        :       $(HFILE) vitesseB.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseB.o vitesseB.c 

vitesseC.o        :       $(HFILE) vitesseC.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseC.o vitesseC.c 

vitesseD.o       :       $(HFILE) vitesseD.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseD.o vitesseD.c 

vitesseE.o        :       $(HFILE) vitesseE.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseE.o vitesseE.c 

vitesseF.o        :       $(HFILE) vitesseF.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseF.o vitesseF.c 

vitesseG.o        :       $(HFILE) vitesseG.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseG.o vitesseG.c 

vitesseH.o        :       $(HFILE) vitesseH.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseH.o vitesseH.c 

vitesseI.o        :       $(HFILE) vitesseI.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseI.o vitesseI.c 

vitesseJ.o        :       $(HFILE) vitesseJ.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseJ.o vitesseJ.c 

vitesseK.o        :       $(HFILE) vitesseK.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseK.o vitesseK.c 

vitesseL.o        :       $(HFILE) vitesseL.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseL.o vitesseL.c 

vitesseM.o        :       $(HFILE) vitesseM.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseM.o vitesseM.c 

vitesseN.o        :       $(HFILE) vitesseN.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseN.o vitesseN.c 

vitesseO.o        :       $(HFILE) vitesseO.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseO.o vitesseO.c 

vitesseP.o        :       $(HFILE) vitesseP.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseP.o vitesseP.c 

vitesseQ.o        :       $(HFILE) vitesseQ.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseQ.o vitesseQ.c 

vitesseR.o        :       $(HFILE) vitesseR.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseR.o vitesseR.c 

vitesseS.o        :       $(HFILE) vitesseS.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseS.o vitesseS.c 

vitesseT.o        :       $(HFILE) vitesseT.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseT.o vitesseT.c 

vitesseU.o        :       $(HFILE) vitesseU.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseU.o vitesseU.c 

vitesseV.o        :       $(HFILE) vitesseV.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseV.o vitesseV.c 

vitesseW.o        :       $(HFILE) vitesseW.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseW.o vitesseW.c 

vitesseX.o        :       $(HFILE) vitesseX.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseX.o vitesseX.c 

vitesseY.o        :       $(HFILE) vitesseY.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseY.o vitesseY.c 

vitesseZ.o        :       $(HFILE) vitesseZ.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesseZ.o vitesseZ.c 

vitesse1.o        :       $(HFILE) vitesse1.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesse1.o vitesse1.c 

vitesse2.o        :       $(HFILE) vitesse2.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesse2.o vitesse2.c 

vitesse3.o        :       $(HFILE) vitesse1.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesse3.o vitesse3.c 

vitesse4.o        :       $(HFILE) vitesse1.c 
	$(CC) $(CFLAGS)  $(GCCFLAG) -o vitesse4.o vitesse4.c 

v_ptype.h			:	 
	touch v_ptype.h 

v_var.h			: 
	touch v_var.h  

v_extern.h		: 
	touch v_extern.h 

v_prog.h		:  
	rm -f *.o

clean	:
	rm -f *.o  core


