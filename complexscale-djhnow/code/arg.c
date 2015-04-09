

//   BE CAREFUL USING NEXTRAN!!   WITH MPI ALL PROCS MUST BY SYNCED.

#ifdef PGFFLAG

#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

    int myiargc_() {         
      extern int    __argc_save;
      return __argc_save - 1;
    } 

    void mygetarg_(int* i, char** buffer) {          
      extern char **__argv_save;

      printf("Go mygetarg! \n");

      strcpy(*buffer, __argv_save[*i+1]);

    //	printf("Arg %i is %50c \n",*i+1,*buffer);
    } 
#endif



