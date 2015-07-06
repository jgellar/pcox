/* copied @ 2015-07-06 from survival version 2.38-3
/* code written by Terry Terneau, copied & modified by Fabian Scheipl
/* modifications marked by "fs"

/* Automatically generated from all.nw using noweb */
#include "survS.h"

SEXP coxcount1(SEXP y2, SEXP strat2) {
    int ntime, nrow;
    int i, j, n;
    int stratastart=0;  /* start row for this strata */
    int nrisk=0;  /* number at risk (=0 to stop -Wall complaint)*/
    double *time, *status;
    int *strata;
    double dtime;
    SEXP rlist, rlistnames, rtime, rn, rindex, rstatus;
    int *rrindex, *rrstatus;
    
    n = nrows(y2);
    time = REAL(y2);
    status = time +n;
    strata = INTEGER(strat2);
    
    /* 
    ** First pass: count the total number of death times (risk sets)
    **  and the total number of rows in the new data set.
    */
    ntime=0; nrow=0;
    for (i=0; i<n; i++) {
        if (strata[i] ==1) nrisk =0;
        nrisk++;
        if (status[i] ==1) {
            ntime++;
            dtime = time[i];
            /* walk across tied times, if any */
            for (j=i+1; j<n && time[j]==dtime && status[j]==1 && strata[j]==0;
                 j++) nrisk++;
            i = j-1;
            nrow += nrisk;
        }
    }
    /*
    **  Allocate memory
    */
    PROTECT(rtime = allocVector(REALSXP, ntime));
    PROTECT(rn = allocVector(INTSXP, ntime));
    PROTECT(rindex=allocVector(INTSXP, nrow));
    PROTECT(rstatus=allocVector(INTSXP,nrow));
    rrindex = INTEGER(rindex);
    rrstatus= INTEGER(rstatus);
    
    /*
    ** Pass 2, fill them in
    */
    ntime=0; 
    for (i=0; i<n; i++) {
        if (strata[i] ==1) stratastart =i;
        if (status[i]==1) {
            dtime = time[i];
            for (j=stratastart; j<i; j++) *rrstatus++=0; /*non-deaths */
            *rrstatus++ =1; /* this death */
            /* tied deaths */
            for(j= i+1; j<n && status[j]==1 && time[j]==dtime  && strata[j]==0;
                j++) *rrstatus++ =1;
            i = j-1;

            REAL(rtime)[ntime] = dtime;
            INTEGER(rn)[ntime] = i +1 -stratastart;
            ntime++;
            for (j=stratastart; j<=i; j++) *rrindex++ = j+1;
            }
    }
    /* return the list */
    PROTECT(rlist = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(rlist, 0, rn);
    SET_VECTOR_ELT(rlist, 1, rtime);
    SET_VECTOR_ELT(rlist, 2, rindex);
    SET_VECTOR_ELT(rlist, 3, rstatus);
    PROTECT(rlistnames = allocVector(STRSXP, 4));
    SET_STRING_ELT(rlistnames, 0, mkChar("nrisk"));
    SET_STRING_ELT(rlistnames, 1, mkChar("time"));
    SET_STRING_ELT(rlistnames, 2, mkChar("index"));
    SET_STRING_ELT(rlistnames, 3, mkChar("status"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(6);
    return(rlist);
}
