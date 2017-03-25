/* 
Unit:
                         1 length = 1 Mpc/h
                         1 mass = 1 M_sun /h
*/

/* Unit transformation constants */
#define SUNTOKG               1.9891e30
#define MPCTOM                3.08567758e22
#define G0                    6.67384e-11
#define PI                    3.14159265358979323

/* ST parameters */
#define Ast                   0.3222
#define pst                   0.3
#define qst                   0.707

/* Mass function transfer fitting param */
#define a_m                   -0.52569
#define b_m                   0.88117            /* Mx = 10^b_m */

/* c-M relation transfer fitting param */
#define a_cm                  -0.27127
#define b_cm                  1.31206            /* Mx = 10^b_cm */

/* 1D Spline related parameters */
#define Mmin                  1                  /* Mmin < 1e5 is enough for Kmax < 100 h/Mpc */
#define Mmax                  1e19               /* Mmax > 1e18 is enough */
#define NI                    1500               /* the interpolation point number */
