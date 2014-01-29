/* LSQR.c
* I removed all the comments from the original 
* file because they were causing trouble with the 
* compiler. The original file can be found at
* http://www.stanford.edu/group/SOL/software/lsqr/c/clsqr2.zip
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "lsqr.h"

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif

#define ZERO   0.0
#define ONE    1.0

static double
d2norm( const double a, const double b )
{
    double scale;
    const double zero = 0.0;

    scale  = fabs( a ) + fabs( b );
    if (scale == zero)
        return zero;
    else
        return scale * sqrt( (a/scale)*(a/scale) + (b/scale)*(b/scale) );
}

static void
dload( const int n, const double alpha, double x[] )
{    
    int i;
    for (i = 0; i < n; i++) x[i] = alpha;
    return;
}

void lsqr( 
          int m,
          int n,
          void (*aprod)(int mode, int m, int n, double x[], double y[],
                        void *UsrWrk),
          double damp,
          void   *UsrWrk,
          double u[],     
          double v[],     
          double w[],     
          double x[],     
          double se[],    
          double atol,
          double btol,
          double conlim,
          int    itnlim,
          FILE   *nout,

          int    *istop_out,
          int    *itn_out,
          double *anorm_out,
          double *acond_out,
          double *rnorm_out,
          double *arnorm_out,
          double *xnorm_out
         )
{


    int
        istop  = 0,
        itn    = 0;
    double
        anorm  = ZERO,
        acond  = ZERO,
        rnorm  = ZERO,
        arnorm = ZERO,
        xnorm  = ZERO;


    const bool
        extra  = false,       
        damped = damp > ZERO,
        wantse = se != NULL;
    int
        i, maxdx, nconv, nstop;
    double
        alfopt, alpha, arnorm0, beta, bnorm,
        cs, cs1, cs2, ctol,
        delta, dknorm, dnorm, dxk, dxmax,
        gamma, gambar, phi, phibar, psi,
        res2, rho, rhobar, rhbar1,
        rhs, rtol, sn, sn1, sn2,
        t, tau, temp, test1, test2, test3,
        theta, t1, t2, t3, xnorm1, z, zbar;
    char
        enter[] = "Enter LSQR.  ",
        exit[]  = "Exit  LSQR.  ",
        msg[6][100] =
        {
            {"The exact solution is  x = 0"},
            {"A solution to Ax = b was found, given atol, btol"},
            {"A least-squares solution was found, given atol"},
            {"A damped least-squares solution was found, given atol"},
            {"Cond(Abar) seems to be too large, given conlim"},
            {"The iteration limit was reached"}
        };

    char fmt_1000[] = 
        " %s        Least-squares solution of  Ax = b\n"
        " The matrix  A  has %7d rows  and %7d columns\n"
        " damp   = %-22.2e    wantse = %10i\n"
        " atol   = %-22.2e    conlim = %10.2e\n"
        " btol   = %-22.2e    itnlim = %10d\n\n";
    char fmt_1200[] =
        "    Itn       x(1)           Function"
        "     Compatible    LS      Norm A   Cond A\n";
    char fmt_1300[] =
        "    Itn       x(1)           Function"
        "     Compatible    LS      Norm Abar   Cond Abar\n";
    char fmt_1400[] =
        "     phi    dknorm  dxk  alfa_opt\n";
    char fmt_1500_extra[] =
        " %6d %16.9e %16.9e %9.2e %9.2e %8.1e %8.1e %8.1e %7.1e %7.1e %7.1e\n";
    char fmt_1500[] =
        " %6d %16.9e %16.9e %9.2e %9.2e %8.1e %8.1e\n";
    char fmt_1550[] =
        " %6d %16.9e %16.9e %9.2e %9.2e\n";
    char fmt_1600[] = 
        "\n";
    char fmt_2000[] =
        "\n"
        " %s       istop  = %-10d      itn    = %-10d\n"
        " %s       anorm  = %11.5e     acond  = %11.5e\n"
        " %s       vnorm  = %11.5e     xnorm  = %11.5e\n"
        " %s       rnorm  = %11.5e     arnorm = %11.5e\n";
    char fmt_2100[] =
        " %s       max dx = %7.1e occured at itn %-9d\n"
        " %s              = %7.1e*xnorm\n";
    char fmt_3000[] =
        " %s       %s\n";


    if (nout != NULL)
        fprintf(nout, fmt_1000,
                enter, m, n, damp, wantse,
                atol, conlim, btol, itnlim);

    itn    =   0;
    istop  =   0;
    nstop  =   0;
    maxdx  =   0;
    ctol   =   ZERO;
    if (conlim > ZERO) ctol = ONE / conlim;
    anorm  =   ZERO;
    acond  =   ZERO;
    dnorm  =   ZERO;
    dxmax  =   ZERO;
    res2   =   ZERO;
    psi    =   ZERO;
    xnorm  =   ZERO;
    xnorm1 =   ZERO;
    cs2    = - ONE;
    sn2    =   ZERO;
    z      =   ZERO;


    dload( n, 0.0, v );
    dload( n, 0.0, x );

    if ( wantse )
        dload( n, 0.0, se );
    
    alpha  =   ZERO;
    beta   =   cblas_dnrm2 ( m, u, 1 );

    if (beta > ZERO) {
        cblas_dscal ( m, (ONE / beta), u, 1 );
        aprod ( 2, m, n, v, u, UsrWrk );
        alpha  =   cblas_dnrm2 ( n, v, 1 );
    }

    if (alpha > ZERO) {
        cblas_dscal ( n, (ONE / alpha), v, 1 );
        cblas_dcopy ( n, v, 1, w, 1 );
    }

    arnorm = arnorm0 = alpha * beta;
    if (arnorm == ZERO) goto goto_800;
    
    rhobar =   alpha;
    phibar =   beta;
    bnorm  =   beta;
    rnorm  =   beta;

    if (nout != NULL) {
        if ( damped ) 
            fprintf(nout, fmt_1300);
        else
            fprintf(nout, fmt_1200);

        test1  = ONE;
        test2  = alpha / beta;
        
        if ( extra ) 
            fprintf(nout, fmt_1400);

        fprintf(nout, fmt_1550, itn, x[0], rnorm, test1, test2);
        fprintf(nout, fmt_1600);
    }


    while (1) {
        itn    = itn + 1;
        

        cblas_dscal ( m, (- alpha), u, 1 );
        aprod ( 1, m, n, v, u, UsrWrk );
        beta   =   cblas_dnrm2 ( m, u, 1 );


        temp   =   d2norm( alpha, beta );
        temp   =   d2norm( temp , damp );
        anorm  =   d2norm( anorm, temp );

        if (beta > ZERO) {
            cblas_dscal ( m, (ONE / beta), u, 1 );
            cblas_dscal ( n, (- beta), v, 1 );
            aprod ( 2, m, n, v, u, UsrWrk );
            alpha  =   cblas_dnrm2 ( n, v, 1 );
            if (alpha > ZERO) {
                cblas_dscal ( n, (ONE / alpha), v, 1 );
            }
        }

        rhbar1 = rhobar;
        if ( damped ) {
            rhbar1 = d2norm( rhobar, damp );
            cs1    = rhobar / rhbar1;
            sn1    = damp   / rhbar1;
            psi    = sn1 * phibar;
            phibar = cs1 * phibar;
        }

        rho    =   d2norm( rhbar1, beta );
        cs     =   rhbar1 / rho;
        sn     =   beta   / rho;
        theta  =   sn * alpha;
        rhobar = - cs * alpha;
        phi    =   cs * phibar;
        phibar =   sn * phibar;
        tau    =   sn * phi;

        t1     =   phi   / rho;
        t2     = - theta / rho;
        t3     =   ONE   / rho;
        dknorm =   ZERO;

        if ( wantse ) {
            for (i = 0; i < n; i++) {
                t      =  w[i];
                x[i]   =  t1*t  +  x[i];
                w[i]   =  t2*t  +  v[i];
                t      = (t3*t)*(t3*t);
                se[i]  =  t     +  se[i];
                dknorm =  t     +  dknorm;
            }
        }
        else {
            for (i = 0; i < n; i++) {
                t      =  w[i];
                x[i]   =  t1*t  +  x[i];
                w[i]   =  t2*t  +  v[i];
                dknorm = (t3*t)*(t3*t)  +  dknorm;
            }
        }

        dknorm = sqrt( dknorm );
        dnorm  = d2norm( dnorm, dknorm );
        dxk    = fabs( phi * dknorm );
        if (dxmax < dxk ) {
            dxmax   =  dxk;
            maxdx   =  itn;
        }

        delta  =   sn2 * rho;
        gambar = - cs2 * rho;
        rhs    =   phi    - delta * z;
        zbar   =   rhs    / gambar;
        xnorm  =   d2norm( xnorm1, zbar  );
        gamma  =   d2norm( gambar, theta );
        cs2    =   gambar / gamma;
        sn2    =   theta  / gamma;
        z      =   rhs    / gamma;
        xnorm1 =   d2norm( xnorm1, z     );

        acond  =   anorm * dnorm;
        res2   =   d2norm( res2 , psi    );
        rnorm  =   d2norm( res2 , phibar );
        arnorm =   alpha * fabs( tau );


        alfopt =   sqrt( rnorm / (dnorm * xnorm) );
        test1  =   rnorm /  bnorm;
        test2  =   ZERO;
        if (rnorm   > ZERO) test2 = arnorm / (anorm * rnorm);
        test3  =   ONE   /  acond;
        t1     =   test1 / (ONE  +  anorm * xnorm / bnorm);
        rtol   =   btol  +  atol *  anorm * xnorm / bnorm;

        t3     =   ONE + test3;
        t2     =   ONE + test2;
        t1     =   ONE + t1;
        if (itn >= itnlim) istop = 5;
        if (t3  <= ONE   ) istop = 4;
        if (t2  <= ONE   ) istop = 2;
        if (t1  <= ONE   ) istop = 1;

        if (test3 <= ctol) istop = 4;
        if (test2 <= atol) istop = 2;
        if (test1 <= rtol) istop = 1; 

        if (nout  == NULL     ) goto goto_600;
        if (n     <= 40       ) goto goto_400;
        if (itn   <= 10       ) goto goto_400;
        if (itn   >= itnlim-10) goto goto_400;
        if (itn % 10 == 0     ) goto goto_400;
        if (test3 <=  2.0*ctol) goto goto_400;
        if (test2 <= 10.0*atol) goto goto_400;
        if (test1 <= 10.0*rtol) goto goto_400;
        if (istop != 0        ) goto goto_400;
        goto goto_600;

    goto_400:
        if ( extra ) {
            fprintf(nout, fmt_1500_extra,
                    itn, x[0], rnorm, test1, test2, anorm,
                    acond, phi, dknorm, dxk, alfopt);
        }
        else {
            fprintf(nout, fmt_1500,
                    itn, x[0], rnorm, test1, test2, anorm, acond);
        }
        if (itn % 10 == 0) fprintf(nout, fmt_1600);

    goto_600:
        if (istop == 0) {
            nstop  = 0;
        }
        else {
            nconv  = 1;
            nstop  = nstop + 1;
            if (nstop < nconv  &&  itn < itnlim) istop = 0;
        }

        if (istop != 0) break;
        
    }

    if ( wantse ) {
        t    =   ONE;
        if (m > n)     t = m - n;
        if ( damped )  t = m;
        t    =   rnorm / sqrt( t );
      
        for (i = 0; i < n; i++)
            se[i]  = t * sqrt( se[i] );
        
    }

 goto_800:
    if (damped  &&  istop == 2) istop = 3;
    if (nout != NULL) {
        fprintf(nout, fmt_2000,
                exit, istop, itn,
                exit, anorm, acond,
                exit, bnorm, xnorm,
                exit, rnorm, arnorm);
        fprintf(nout, fmt_2100,
                exit, dxmax, maxdx,
                exit, dxmax/(xnorm + 1.0e-20));
        fprintf(nout, fmt_3000,
                exit, msg[istop]);
    }

    *istop_out  = istop;
    *itn_out    = itn;
    *anorm_out  = anorm;
    *acond_out  = acond;
    *rnorm_out  = rnorm;
    *arnorm_out = test2;
    *xnorm_out  = xnorm;

    return;
}


