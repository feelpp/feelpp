#ifndef CL_BASE_CONFIG_H
#define CL_BASE_CONFIG_H

/* CL_GETTIMEOFDAY */
/* Define if you have the gettimeofday() function. */
#cmakedefine HAVE_GETTIMEOFDAY
/* Define if the declaration of gettimeofday() needs dots. */
#cmakedefine GETTIMEOFDAY_DOTS
/* Define as the type of `tzp' in gettimeofday() declaration. */
#cmakedefine GETTIMEOFDAY_TZP_T

/* CL_TIMES_CLOCK */
/* Define if you have the times() function and it returns the real time,
   but don't have the gettimeofday() function. */
#cmakedefine HAVE_TIMES_CLOCK
#endif /* CL_BASE_CONFIG_H */

