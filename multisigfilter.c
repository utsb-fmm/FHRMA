//123456789012345678901234567890123456789012345678901234567890123456789012
//
//     FHR Morphological Analysis Toolbox  Copyright (C) 2018 Samuel Boudet, Faculté de Médecine et Maïeutique,
//     samuel.boudet@gmail.com
//
//     This file is part of FHR Morphological Analysis Toolbox 
//
//     FHR Morphological Analysis Toolbox  is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     FHR Morphological Analysis Toolbox  is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <process.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

void filter(double *fdata, double *data,int samples,double *b,double *a,int order,double *zi);
unsigned __stdcall threadfunc(void *arg);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]) {
    int samples, nchan, order;
    double *fdata, *data, *a, *b, *zi;
    int processMax=2;
    int chancount[1];
    int i;
    void *args[9];
    HANDLE *hThread;
    unsigned *threadID;

    if (nrhs ==5) {

        processMax =((int) (*mxGetPr(prhs[4]) ));
        hThread=malloc(processMax*sizeof(HANDLE));
        threadID=malloc(processMax*sizeof(unsigned));
        

    }else if (nrhs != 4) {
        mexErrMsgTxt("Incorrect number of input arguments");
    }
    samples = ((int) mxGetM(prhs[2]));
    nchan = ((int) mxGetN(prhs[2]));
    order = ((int) mxGetN(prhs[0]));
    
    if ( mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Complex values are not supported");
    }
    if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
        mexErrMsgTxt("Only double values are supported");
    
    b = mxGetPr(prhs[0]);
    a = mxGetPr(prhs[1]);
    data = mxGetPr(prhs[2]);
    zi = mxGetPr(prhs[3]);
    
    plhs[0]=mxCreateDoubleMatrix(samples, nchan, mxREAL);
    fdata = mxGetPr(plhs[0]);
    
    *chancount=-1;
    for(i=0;i<processMax;i++) {
        args[0]=data;args[1]=fdata;args[2]=chancount;args[3]=&nchan;
        args[4]=&samples;args[5]=b;args[6]=a;args[7]=&order;args[8]=zi;
        hThread[i] = (HANDLE)_beginthreadex( NULL, 0, &threadfunc, args, 0, &threadID[i] );
    }
    
    for(i=0;i<processMax;i++) {
        WaitForSingleObject( hThread[i], INFINITE );
        CloseHandle( hThread[i] );
    }
    
    free(hThread);
    free(threadID);
    return;
}

unsigned __stdcall threadfunc(void *arg) {
    void **args=(void **) arg;

    double *data=args[0];
    double *fdata=args[1];
    int *chancount=args[2];
    int nchans=*((int *) args[3]);
    int samples=*((int *) args[4]);
    double *b=args[5];
    double *a=args[6];
    int order=*((int *) args[7]);
    double *zi=args[8];
    int n;
    double g=0;
    while( (n=++(*chancount))<nchans ) {
        filter(&fdata[n*samples],&data[n*samples],samples,b,a,order,zi);
    }
    _endthreadex( 0 );
    return 0;
}

void filter(double *fdata, double *data,int samples,double *b,double *a,int order,double *zi) {
    double *d,*dinit;
    int i,j,nfact;
    nfact=3*(order-1); // Prolongation for minimizing border effect
    d=malloc((samples+nfact)*sizeof(double));
    dinit=malloc(nfact*sizeof(double));
    //*****************First pass *************************
    //Initial values***************************************
    //data[i] (i>=0?data[i]:2*data[0]-data[-i])
    //d[i]= (i>=0?d[i]:dinit[nfact+i])
    
    for(i=0;i<order-1;i++)
        dinit[i]=(2*data[0]-data[nfact])*zi[i];
    for(i=order-1;i<nfact;i++)
        dinit[i]=0;
    
    for(i=-nfact;i<0;i++) {
        for(j=MIN(order-1,i+nfact);j>=1;j--) {
            dinit[nfact+i]=b[j]*(2*data[0]-data[-i+j])+dinit[nfact+i]-a[j]*dinit[nfact+i-j];
        }
        dinit[nfact+i]+=b[0]*(2*data[0]-data[-i]);
    }
    
    for(i=0;i<order-1;i++) {
        d[i]=0;
        for(j=order-1;j>i && j>=1;j--) {
            d[i]=b[j]*(2*data[0]-data[-i+j])+d[i]-a[j]*dinit[nfact+i-j];
        }
        for(j=i;j>=1;j--) {
            d[i]=b[j]*data[i-j]+d[i]-a[j]*d[i-j];
        }
        d[i]+=b[0]*data[i];
    }
    
    
    
    //Normal values***************************************
    for(i=order-1;i<samples;i++) {
        d[i]=0;
        for(j=order-1;j>=1;j--) {
            d[i]=b[j]*data[i-j]+d[i]-a[j]*d[i-j];
        }
        d[i]+=b[0]*data[i];
    }    
    
    //End values********************************************
    //data[i] (i<samples?data[i]:2*data[samples-1]-data[2*samples-2-i])
    for(i=samples;i<samples+nfact;i++) {
        d[i]=0;
        for(j=order-1;j>i-samples && j>=1;j--) {
            d[i]=b[j]*data[i-j]+d[i]-a[j]*d[i-j];
        }
        for(j=MIN(order-1,i-samples);j>=1;j--) {
            d[i]=b[j]*(2*data[samples-1]-data[2*samples-2-i+j])+d[i]-a[j]*d[i-j];
        }
        d[i]+=b[0]*(2*data[samples-1]-data[2*samples-2-i]);
    }
    
    //***********Reverse order************************ 
    //init values ************************************
    for(i=samples+nfact-1;i>=samples;i--) {
        if(samples+nfact-i<order)
            dinit[i-samples]=d[samples+nfact-1]*zi[samples+nfact-i-1];
        else
            dinit[i-samples]=0;
        for(j=MIN(order-1,nfact+samples-i-1);j>=1;j--) {
            dinit[i-samples]=b[j]*d[i+j]+ dinit[i-samples] - a[j]*dinit[i-samples+j];
        }
        dinit[i-samples]+=b[0]*d[i];
    }
    for(i=samples-1;i>samples-order;i--) {
        for(j=order-1;j>=samples-i && j>=1;j--) {
            fdata[i]=b[j]*d[i+j] +fdata[i]- a[j]*dinit[i+j-samples];
        }
        for(j=samples-1-i;j>=1;j--) {
            fdata[i]=b[j]*d[i+j] +fdata[i]- a[j]*fdata[i+j];
        }
        fdata[i]+=b[0]*d[i];
    }
    
    //Normal values***************************************
    for(i=samples-order;i>=0;i--) {
        for(j=order-1;j>=1;j--) {
            fdata[i]=b[j]*d[i+j] +fdata[i]- a[j]*fdata[i+j];
        }
        fdata[i]+=b[0]*d[i];
    }
    free(d);
    free(dinit);
}