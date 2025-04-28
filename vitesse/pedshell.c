/* 
This program will convert an lcp-generated shell file for use by
vitesse. The ouput files are given the prefix v, the program 
unknown is eliminated, and mlink and linkmap are renamed to 
vitesse. 

This program is adapted from a program written by Alan Young.

#define SCREEN
#define SCREEN_P
#define PROTS 
#include <sys/ioctl.h>
#include <signal.h>

#include <termio.h>
#ifdef __sparc
#ifdef __SVR4
#define BROKEN_TIOCGWINSZ
#include <sys/termios.h>
#endif
#endif
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define	INT		int
#define REAL		double	/* INTERNAL STORAGE FOR REAL VARIABLES	*/
#define	ASTERISK	'*'
#define	BACKSLASH	'\\'
#define	DOT		'.'
#define	EOL		'\n'
#define	EOS		0
#define	HASH		'#'
#define	LB		'('
#define	NO		'n'
#define	QUOTE		'\"'
#define	SEMICOLON	';'
#define	SINGLE_QUOTE	'\''
#define	SLASH		'/'
#define	SPACE		' '
#define	TAB		'\t'
#define	UNDER_SCORE	'_'
#define	YES		'y'
/*--------------------------------------------------------------------------*/
FILE	*fi, *fo;					/* FILE POINTERS */
/*--------------------------------------------------------------------------*/
char	w[200];		
char	s_new_word[200];	
/*============================================================================*/
typedef	struct	{		
	char	c[35];
	}	WORD;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

WORD	s_keyword[] = {
	"linkmap",
	"mlink"
	 };

WORD s_vit_keyword[] = {
	"pedcheck",
	"pedcheck"
	 };

WORD	keyword[] = {
	"final.bak",
	"final.out",
	"outfile.dat",
	"stream.bak",
	"stream.dat",
	"stream.out",
	"unknown"
	 };

WORD vit_keyword[] = {
	"vfinal.bak",
	"vfinal.out",
	"voutfile.dat",
	"vstream.bak",
	"vstream.dat",
	"vstream.out",
	" "
	 };


/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#define	NUM_KEYWORDS	(sizeof(keyword)/sizeof(WORD))
#define	NUM_S_KEYWORDS	(sizeof(s_keyword)/sizeof(WORD))

/*--------------------------------------------------------------------------*/
			/* DECIDE IF WORD IS TO BE PRESERVED IN CURRENT FORM  */
int is_keyword( char *word )
{	INT	wlo = -1;
	INT	whi = NUM_KEYWORDS;
	INT	m, sc;
	while ( whi - wlo > 1 )
		{m = ( whi + wlo ) / 2;
		if ( ( sc = strcmp( word, keyword[m].c ) ) > 0 )
			{wlo = m;
			}
		else if ( sc < 0 )
			{whi = m;			
			}
		else
			{return(m);	/* WORD IS IDENTIFIED AS KEYWORD */
			}
		}
	return(-1);			/* WORD CAN SAFELY BE REPLACED */
}



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
			/* DECIDE IF WORD IS TO BE PRESERVED IN CURRENT FORM  */
int is_s_keyword( char *word )
{	INT	wlo = -1;
	INT	whi = NUM_S_KEYWORDS;
	INT	m, sc;
	while ( whi - wlo > 1 )
		{m = ( whi + wlo ) / 2;
		if ( ( sc = strcmp( word, s_keyword[m].c ) ) > 0 )
			{wlo = m;
			}
		else if ( sc < 0 )
			{whi = m;			
			}
		else
			{return(m);	/* WORD IS IDENTIFIED AS KEYWORD */
			}
		}
	return(-1);			/* WORD CAN SAFELY BE REPLACED */
}



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
int process_line( char *s )
{	INT	i = 0;
	INT	j, word_len, word_offset;
	char	c;
	int     cont,count;
	int     len1,len2;
	int val;
        int s_word_offset;
        int s_word_offset2;
        int i_temp,v_change,true,false;

	true = 1;
	false = 0;

	cont = 0;
	count = 0;
	s_word_offset = 0;
	s_word_offset2 = 0;
	word_offset = 0;
	START:
	while ( !isalpha( s[i] ) && UNDER_SCORE != s[i] && DOT != s[i] )/* SEEK NEXT WORD */
		{if ( EOS == s[i] )
			{
			s_new_word[s_word_offset] = s[i];
			return(0);		/* REACHED END OF LINE */
			}
			s_new_word[s_word_offset] = s[i];
		i++;
	        s_word_offset2++;
	        s_word_offset++;

		}

#ifdef SCREEN
fprintf( stdout, " s =  %s  s_new =  %s\n",s,s_new_word ); 	
#endif
	count++;
	word_len = 0;
	word_offset = i;
#ifdef SCREEN
			fprintf( stdout, " word offset =  %d  s word offset = %d\n",word_offset,s_word_offset ); 	
#endif
	while ( isalnum( s[i] ) || UNDER_SCORE == s[i] || DOT == s[i] ) /* LOAD NEXT WORD */
		{w[word_len++] = s[i];
		  s_new_word[s_word_offset + i-word_offset ] = s[i];
		  s_word_offset2++;
		  i++;
		}
	if ( word_len > 1  && HASH != ( c = s[word_offset-1] )
			&& SINGLE_QUOTE != c && !isdigit(c)
			&& UNDER_SCORE != w[0] )
		{w[word_len] = EOS;		/* EXAMINE CANDIDATE WORD */
#ifdef SCREEN
			fprintf( stdout, "%s", w ); 	
#endif

		val =is_keyword( w );
#ifdef SCREEN
			fprintf( stdout, "  %d", val ); 	
			fprintf( stdout, "\n");  
#endif
		if (val  > -1 )
        	{
/*  write to a file;  */
#ifdef SCREEN
			fprintf( stdout, "%s", w ); 	
			fprintf( stdout, "\n");  
			fprintf( stdout, "%s", vit_keyword[val].c ); 	
			fprintf( stdout, "\n");  
#endif
			len2 = strlen(vit_keyword[val].c);
			for ( j = 0; j < len2; j++ )
				{s_new_word[j+s_word_offset] = vit_keyword[val].c[j];
				}		
#ifdef SCREEN
		fprintf( stdout, "s_new after write  %s\n", s_new_word ); 	
#endif
			s_word_offset +=len2;
#ifdef SCREEN
	fprintf( stdout, " i =  %d s_word_offset = %d\n",i,s_word_offset); 	
#endif
			s_word_offset2 = s_word_offset2 + len2;
	      	}
		else if(is_s_keyword(w) > -1)
		{
			val = is_s_keyword(w) ;
			v_change = true;
			i_temp = i;
			while(s[i_temp] != EOL)
			{
                            if(s[i_temp] == BACKSLASH)
			      v_change = false;
			     i_temp++;
			}

			if(v_change == true)
	                {	
			/*
			fprintf( stdout, "%s", &s[i] ); 	
			*/
			len2 = strlen(s_vit_keyword[val].c);
			for ( j = 0; j < len2; j++ )
				{s_new_word[j+s_word_offset] = s_vit_keyword[val].c[j];
				}		
			s_word_offset +=len2;
			s_word_offset2 = s_word_offset2 + len2;
                        } 
			else
			s_word_offset +=word_len;

		}
		else
			s_word_offset +=word_len;
			  
		}
		else
		 {
		 s_word_offset2 +=word_len;
		 s_word_offset +=word_len;
		 }
	word_offset = i - word_offset;	
	goto START;
}
/*============================================================================*/
						
void process_file()
{	INT	lines = 0;
	char	text[202],p[100];
	int    first_char;
	int val;
	int count,stmt;
	int val2;
        int word_len,i;

	count = 0;
	while ( EOF != ( first_char = fgetc(fi) ) )
		{lines++;
		text[0] = first_char;
		if ( EOL != text[0] )
		{
			fgets( text+1, 200, fi );
			process_line( text );
		        fprintf( fo, "%s", s_new_word);		/* WRITE LINE */
		}
		else
			{fprintf( fo, "\n" );
			}
		}
	printf("\n-> %d lines processed\n", lines );
}
/*============================================================================*/
					 	/* GET NAMES OF DATA FILES */
int main( int argc, char *argv[] )
{	int	i = -1;
	char	infile[200], outfile[200];	/* INPUT/OUTPUT SOURCE FILES */
	int   chmod();

	if ( argc >= 2 )
		{strcpy( infile, argv[1] );
		}
	else             /* GET FORMAT FILENAME IF NOT GIVEN ON COMMAND LINE */
		{printf("\nInput file :");
		scanf( "%s", infile );
		}
	if ( argc >= 3 )
		{strcpy( outfile, argv[2] );
		if ( 0 == strcmp( outfile, "*" ) )
			{strcpy( outfile, infile );
			while ( EOS != outfile[++i] )	/* REMOVE VERSION */
				{if ( SEMICOLON == outfile[i] )
					{outfile[i] = EOS;
					}
				}
			}
		}
	else             /* GET FORMAT FILENAME IF NOT GIVEN ON COMMAND LINE */
		{printf("\nOutput file: ");
		scanf( "%s", outfile );
		}
	if ( NULL == ( fi = fopen( infile, "r" ) ) )
		{printf("\nFile \"%s\" cannot be opened", infile );
		exit( 0 );
		}
	fo = fopen( outfile, "w" );
	process_file();
	fclose( fi );
	fclose( fo );
	chmod(outfile,0755);
	exit( 0 );
}
/*============================================================================*/
