/* IMT 10Sep95  Declare jump buffer used to recover from exception exits & aborts */
#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32)
#include <setjmp.h>
extern jmp_buf gRecoverToConsole;
#endif /* Macintosh C compilers  and CW Win32 */


#include "stdio.h"
#include "f2c.h"
#define PAUSESIG 15

#ifdef KR_headers
#define Void /* void */
#define Int /* int */
#else
#define Void void
#define Int int
#undef abs
#undef min
#undef max
#include "stdlib.h"
#include "signal1.h"
#ifdef __cplusplus
extern "C" {
#endif
extern int getpid(void), isatty(int), pause(void);
#endif

extern VOID f_exit(Void);

 static VOID
waitpause(Int n)
{	n = n; /* shut up compiler warning */
	return;
	}

 static VOID
#ifdef KR_headers
s_1paus(fin) FILE *fin;
#else
s_1paus(FILE *fin)
#endif
{
	fprintf(stderr,
	"To resume execution, type go.  Other input will terminate the job.\n");
	fflush(stderr);
	if( getc(fin)!='g' || getc(fin)!='o' || getc(fin)!='\n' ) {
		fprintf(stderr, "STOP\n");
#ifdef NO_ONEXIT
		f_exit();
#endif

/* IMT 10Sep95  Use jump buffer to recover to console instead of exit() */
#if defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) 
		longjmp( gRecoverToConsole, 1 );
#else
		exit(0);
#endif /* Macintosh compilers and CW Win32 */
		}
	}

 int
#ifdef KR_headers
s_paus(s, n) char *s; ftnlen n;
#else
s_paus(char *s, ftnlen n)
#endif
{
	fprintf(stderr, "PAUSE ");
	if(n > 0)
		fprintf(stderr, " %.*s", (int)n, s);
	fprintf(stderr, " statement executed\n");
	if( isatty(fileno(stdin)) )
		s_1paus(stdin);
	else {
#ifdef MSDOS
		FILE *fin;
		fin = fopen("con", "r");
		if (!fin) {
			fprintf(stderr, "s_paus: can't open con!\n");
			fflush(stderr);
			exit(1);
			}
		s_1paus(fin);
		fclose(fin);
		
/* IMT 14Aug97  Create appropriate pause mechanism for MacOS */
#elif defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) 
		fprintf( stderr, "To resume execution, click the mouse button.\n" );
		fflush( stderr );
		pause();			/* Special MacOS version */
		
#else
		fprintf(stderr,
		"To resume execution, execute a   kill -%d %d   command\n",
			PAUSESIG, getpid() );
		signal1(PAUSESIG, waitpause);
		fflush(stderr);
		pause();
#endif
		}
	fprintf(stderr, "Execution resumes after PAUSE.\n");
	fflush(stderr);
	return 0; /* NOT REACHED */
#ifdef __cplusplus
	}
#endif
}
