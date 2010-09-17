/* This handles the 1e, 2e and DF integrals for all SAPT jobs */

#ifdef HAVEMKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
    
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <time.h>
  
#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "sapt.h"

namespace psi { namespace sapt {

void SAPT::compute_amplitudes()
{
  psio_->open(PSIF_SAPT_AMPS,PSIO_OPEN_NEW);

  if (workflow_.t_arar || workflow_.theta_arar) {
    if (params_.print) {
      fprintf(outfile,"Computing T ARAR Amplitudes...\n");
      fflush(outfile);
    }
    t_arar(workflow_.mp2_opdm,workflow_.theta_ar_ar,workflow_.g_arar);
  }

  if (workflow_.t_bsbs || workflow_.theta_bsbs) {  
    if (params_.print) {
      fprintf(outfile,"Computing T BSBS Amplitudes...\n");
      fflush(outfile);
    }
    t_bsbs(workflow_.mp2_opdm,workflow_.theta_bs_bs,workflow_.g_bsbs);
  }

  if (params_.nat_orbs) {
    if (params_.print) {
      fprintf(outfile,"Computing MP2 Natural Orbitals...\n");
      fflush(outfile);
    }
    natural_orbitalify("AA MP2 OPDM","RR MP2 OPDM",calc_info_.evalsA,
      calc_info_.CA,calc_info_.noccA,calc_info_.nvirA,'A');
    natural_orbitalify("BB MP2 OPDM","SS MP2 OPDM",calc_info_.evalsB,
      calc_info_.CB,calc_info_.noccB,calc_info_.nvirB,'B');
  }

  if (workflow_.t_arbs || workflow_.t_bsar) {
    if (params_.print) {
      fprintf(outfile,"Computing T ARBS Amplitudes...\n");
      fflush(outfile);
    }
    t_arbs(workflow_.t_bsar);
  }

  if (workflow_.t2_arar || workflow_.theta2_arar) {
    if (params_.print) {
      fprintf(outfile,"Computing T2 ARAR Amplitudes...\n");
      fflush(outfile);
    }
    t2_arar(workflow_.theta2_ar_ar);
  }

  if (workflow_.t2_bsbs || workflow_.theta2_bsbs) {
    if (params_.print) {
      fprintf(outfile,"Computing T2 BSBS Amplitudes...\n");
      fflush(outfile);
    }
    t2_bsbs(workflow_.theta2_bs_bs);
  }

  if (workflow_.Y2_ar || workflow_.t_ar) {
    if (params_.print) {
      fprintf(outfile,"Computing Y2 AR Amplitudes...\n");
      fflush(outfile);
    }
    Y2("Y2 AR Amplitudes","T AR Amplitudes","RR MP2 OPDM","AA MP2 OPDM",
      "Theta(AR) AR",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",calc_info_.evalsA,calc_info_.noccA,calc_info_.nvirA);
  }

  if (workflow_.Y2_bs || workflow_.t_bs) {
    if (params_.print) {
      fprintf(outfile,"Computing Y2 BS Amplitudes...\n");
      fflush(outfile);
    }
    Y2("Y2 BS Amplitudes","T BS Amplitudes","SS MP2 OPDM","BB MP2 OPDM",
      "Theta(BS) BS",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",calc_info_.evalsB,calc_info_.noccB,calc_info_.nvirB);
  }

  if (workflow_.gt_ar_arbs) {
    if (params_.print) {
      fprintf(outfile,"Computing gARAR x tARBS...\n");
      fflush(outfile);
    }
    g_arar();
  }

  if (workflow_.gt_bs_arbs) {
    if (params_.print) {
      fprintf(outfile,"Computing gBSBS x tARBS...\n");
      fflush(outfile);
    }
    g_bsbs();  
  }

  fprintf(outfile,"\n");
}

void SAPT::t_arar(int opdm, int theta, int garar)
{

  double **B_p_AR = get_AR_ints(1);
  double **tARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nri,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_AR[0][0]),calc_info_.nrio,0.0,&(tARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA);

  if (garar) {

    double **gARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
      calc_info_.noccA*calc_info_.nvirA);

    for (int a=0, ar=0; a < calc_info_.noccA; a++) {
    for (int r=0; r < calc_info_.nvirA; r++, ar++) {
      for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
      for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
        int aar = aa*calc_info_.nvirA + r;
        int arr = a*calc_info_.nvirA + rr;
        gARAR[ar][aarr] = 2.0*tARAR[ar][aarr] - tARAR[arr][aar];
      }}
    }}

    write_IJKL(gARAR,PSIF_SAPT_AMPS,"g ARAR Integrals",calc_info_.noccA*
      calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  }

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
    for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsA[aa]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsA[rr+calc_info_.noccA];
      tARAR[ar][aarr] /= denom;
    }}
  }}

  double **thetaARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
    for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
      int aar = aa*calc_info_.nvirA + r;
      int arr = a*calc_info_.nvirA + rr;
      thetaARAR[ar][aarr] = 2.0*tARAR[ar][aarr] - tARAR[arr][aar];
    }}
  }}

  if (opdm) {

    double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);
    double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);
  
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccA*
      calc_info_.nvirA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
      calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirA,&(tARAR[0][0]),
      calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirA,0.0,&(xAA[0][0]),
      calc_info_.noccA);
    
    C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccA*
      calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
      calc_info_.nvirA,&(tARAR[0][0]),calc_info_.nvirA,0.0,&(xRR[0][0]),
      calc_info_.nvirA);
  
    write_IJKL(xAA,PSIF_SAPT_AMPS,"AA MP2 OPDM",calc_info_.noccA,
      calc_info_.noccA);
    write_IJKL(xRR,PSIF_SAPT_AMPS,"RR MP2 OPDM",calc_info_.nvirA,
      calc_info_.nvirA);

  }

  if (theta) {

    double **tAR_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
      calc_info_.nrio);

    C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
      calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
      calc_info_.noccA*calc_info_.nvirA,&(B_p_AR[0][0]),calc_info_.nrio,
      0.0,&(tAR_p_AR[0][0]),calc_info_.nrio);

    write_IJKL(tAR_p_AR,PSIF_SAPT_AMPS,"Theta(AR) AR",calc_info_.noccA*
      calc_info_.nvirA,calc_info_.nrio);

  }

  free_block(B_p_AR);

  write_IJKL(tARAR,PSIF_SAPT_AMPS,"T ARAR Amplitudes",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);
  write_IJKL(thetaARAR,PSIF_SAPT_AMPS,"T ARAR Antisym Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

}

void SAPT::t_bsbs(int opdm, int theta, int gbsbs)
{

  double **B_p_BS = get_BS_ints(1);
  double **tBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,1.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(tBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  if (gbsbs) {

    double **gBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
      calc_info_.noccB*calc_info_.nvirB);

    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
      for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
        int bbs = bb*calc_info_.nvirB + s;
        int bss = b*calc_info_.nvirB + ss;
        gBSBS[bs][bbss] = 2.0*tBSBS[bs][bbss] - tBSBS[bss][bbs];
      }}
    }}

    write_IJKL(gBSBS,PSIF_SAPT_AMPS,"g BSBS Integrals",calc_info_.noccB*
      calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);
  
  }

  for (int b=0, bs=0; b < calc_info_.noccB; b++) {
  for (int s=0; s < calc_info_.nvirB; s++, bs++) {
    for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
    for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
      double denom = calc_info_.evalsB[b]+calc_info_.evalsB[bb]-
        calc_info_.evalsB[s+calc_info_.noccB]-
        calc_info_.evalsB[ss+calc_info_.noccB];
      tBSBS[bs][bbss] /= denom;
    }}
  }}

  double **thetaBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);

  for (int b=0, bs=0; b < calc_info_.noccB; b++) {
  for (int s=0; s < calc_info_.nvirB; s++, bs++) {
    for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
    for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
      int bbs = bb*calc_info_.nvirB + s;
      int bss = b*calc_info_.nvirB + ss;
      thetaBSBS[bs][bbss] = 2.0*tBSBS[bs][bbss] - tBSBS[bss][bbs];
    }}
  }}

  if (opdm) {

    double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);
    double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);
  
    C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.noccB*
      calc_info_.nvirB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
      calc_info_.noccB*calc_info_.nvirB*calc_info_.nvirB,&(tBSBS[0][0]),
      calc_info_.noccB*calc_info_.nvirB*calc_info_.nvirB,0.0,&(xBB[0][0]),
      calc_info_.noccB);
    
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccB*
      calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
      calc_info_.nvirB,&(tBSBS[0][0]),calc_info_.nvirB,0.0,&(xSS[0][0]),
      calc_info_.nvirB);
  
    write_IJKL(xBB,PSIF_SAPT_AMPS,"BB MP2 OPDM",calc_info_.noccB,
      calc_info_.noccB);
    write_IJKL(xSS,PSIF_SAPT_AMPS,"SS MP2 OPDM",calc_info_.nvirB,
      calc_info_.nvirB);

  }

  if (theta) {

    double **tBS_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
      calc_info_.nrio);

    C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
      calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
      calc_info_.noccB*calc_info_.nvirB,&(B_p_BS[0][0]),calc_info_.nrio,
      0.0,&(tBS_p_BS[0][0]),calc_info_.nrio);

    write_IJKL(tBS_p_BS,PSIF_SAPT_AMPS,"Theta(BS) BS",calc_info_.noccB*
      calc_info_.nvirB,calc_info_.nrio);

  }

  free_block(B_p_BS);

  write_IJKL(tBSBS,PSIF_SAPT_AMPS,"T BSBS Amplitudes",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);
  write_IJKL(thetaBSBS,PSIF_SAPT_AMPS,"T BSBS Antisym Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

}

void SAPT::t_arbs(int bsar)
{
  double **B_p_AR = get_AR_ints(1);
  double **B_p_BS = get_BS_ints(1);
  double **tARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(tARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int b=0, bs=0; b < calc_info_.noccB; b++) {
    for (int s=0; s < calc_info_.nvirB; s++, bs++) {
      double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
        calc_info_.evalsA[r+calc_info_.noccA]-
        calc_info_.evalsB[s+calc_info_.noccB];
      tARBS[ar][bs] /= denom;
    }}
  }} 

  double **tAR_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  C_DGEMM('T','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
    calc_info_.noccA*calc_info_.nvirA,1.0,&(tARBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(B_p_AR[0][0]),calc_info_.nrio,
    0.0,&(tAR_p_BS[0][0]),calc_info_.nrio);

  write_IJKL(tAR_p_BS,PSIF_SAPT_AMPS,"T(AR) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);

  double **tBS_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
    calc_info_.noccB*calc_info_.nvirB,1.0,&(tARBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(B_p_BS[0][0]),calc_info_.nrio,
    0.0,&(tBS_p_AR[0][0]),calc_info_.nrio);

  write_IJKL(tBS_p_AR,PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);

  free_block(B_p_AR);
  free_block(B_p_BS);

  if (bsar) {
    double **tBSAR = block_matrix(calc_info_.noccB*calc_info_.nvirB,
      calc_info_.noccA*calc_info_.nvirA);

    for (int a=0, ar=0; a < calc_info_.noccA; a++) {
    for (int r=0; r < calc_info_.nvirA; r++, ar++) {
      for (int b=0, bs=0; b < calc_info_.noccB; b++) {
      for (int s=0; s < calc_info_.nvirB; s++, bs++) {
        tBSAR[bs][ar] = tARBS[ar][bs];
      }}
    }}

    write_IJKL(tBSAR,PSIF_SAPT_AMPS,"T BSAR Amplitudes",calc_info_.noccB*
      calc_info_.nvirB,calc_info_.noccA*calc_info_.nvirA);
  }

  write_IJKL(tARBS,PSIF_SAPT_AMPS,"T ARBS Amplitudes",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
}

void SAPT::Y2(char *Y2_out, char *T_out, char *VV_opdm, char *OO_opdm, 
  char *theta_OV, int dfnum, char *OO_label, char *OV_label, char *VV_label, 
  double *evals, int nocc, int nvir)
{
  double **Y2 = block_matrix(nocc,nvir);

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,VV_opdm,nvir,nvir);
  double **xAA = read_IJKL(PSIF_SAPT_AMPS,OO_opdm,nocc,nocc);
  
  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,theta_OV,nocc*nvir,
    calc_info_.nrio);

  double **B_p_RR = get_DF_ints(dfnum,VV_label,nvir*nvir);
  double **B_p_AA = get_DF_ints(dfnum,OO_label,nocc*nocc);

  C_DGEMM('N','T',nocc,nvir,nvir*calc_info_.nri,1.0,&(T_p_AR[0][0]),
    nvir*calc_info_.nrio,&(B_p_RR[0][0]),nvir*calc_info_.nrio,0.0,&(Y2[0][0]),
    nvir);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','T',nocc,nvir,calc_info_.nri,-1.0,&(B_p_AA[a*nocc][0]),
      calc_info_.nrio,&(T_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(Y2[0][0]),
      nvir);
  }

  free_block(T_p_AR);

  double **T = block_matrix(nocc,nvir);

  C_DCOPY(nocc*nvir,Y2[0],1,T[0],1);

  for (int a=0; a<nocc; a++) {
    for (int r=0; r<nvir; r++) {
      double denom = evals[a] - evals[r+nocc];
      T[a][r] /= denom;
  }}

  write_IJKL(T,PSIF_SAPT_AMPS,T_out,nocc,nvir);

  double **B_p_AR = get_DF_ints(dfnum,OV_label,nocc*nvir);
  double **C_p_AR = block_matrix(nocc*nvir,calc_info_.nrio);
  double *C_p = init_array(calc_info_.nrio);

  C_DGEMV('t',nvir*nvir,calc_info_.nrio,2.0,&(B_p_RR[0][0]),calc_info_.nrio,
    &(xRR[0][0]),1,0.0,C_p,1);

  C_DGEMV('t',nocc*nocc,calc_info_.nrio,-2.0,&(B_p_AA[0][0]),calc_info_.nrio,
    &(xAA[0][0]),1,1.0,C_p,1);

  C_DGEMV('n',nocc*nvir,calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,
    C_p,1,1.0,&(Y2[0][0]),1);

  free(C_p);

  for (int a=0; a<nocc; a++) {
    C_DGEMM('N','N',nvir,calc_info_.nrio,nvir,1.0,&(xRR[0][0]),nvir,
      &(B_p_AR[a*nvir][0]),calc_info_.nrio,0.0,&(C_p_AR[a*nvir][0]),
      calc_info_.nrio);
  }

  C_DGEMM('N','T',nocc,nvir,nvir*calc_info_.nrio,-1.0,&(C_p_AR[0][0]),
    nvir*calc_info_.nrio,&(B_p_RR[0][0]),nvir*calc_info_.nrio,1.0,&(Y2[0][0]),
    nvir);

  free_block(C_p_AR);
  double **C_p_AA = block_matrix(nocc*nocc,calc_info_.nrio);

  C_DGEMM('N','N',nocc,nocc*calc_info_.nrio,nocc,1.0,&(xAA[0][0]),nocc,
    &(B_p_AA[0][0]),nocc*calc_info_.nrio,0.0,&(C_p_AA[0][0]),nocc*
    calc_info_.nrio);

  for (int a=0; a<nocc; a++) { 
    C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,1.0,&(C_p_AA[a*nocc][0]),
      calc_info_.nrio,&(B_p_AR[a*nvir][0]),calc_info_.nrio,1.0,&(Y2[0][0]),
      nvir);
  }

  free_block(xAA);
  free_block(xRR);
  free_block(C_p_AA);
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);

  write_IJKL(Y2,PSIF_SAPT_AMPS,Y2_out,nocc,nvir);
}

void SAPT::t2_arar(int theta)
{
  double *t2ARAR = t2_solver("T ARAR Amplitudes","Theta(AR) AR",
    PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals","RR RI Integrals",
    calc_info_.evalsA,calc_info_.noccA,calc_info_.nvirA);

  double **thetaARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);

  for (int a=0, ar=0; a < calc_info_.noccA; a++) {
  for (int r=0; r < calc_info_.nvirA; r++, ar++) {
    for (int aa=0, aarr=0; aa < calc_info_.noccA; aa++) {
    for (int rr=0; rr < calc_info_.nvirA; rr++, aarr++) {
      int aar = aa*calc_info_.nvirA + r;
      int arr = a*calc_info_.nvirA + rr;
      int araarr = ar*calc_info_.noccA*calc_info_.nvirA+aarr;
      int aararr = aar*calc_info_.noccA*calc_info_.nvirA+arr;
      thetaARAR[ar][aarr] = 2.0*t2ARAR[araarr] - t2ARAR[aararr];
    }}
  }}

  if (theta) {

    double **B_p_AR = get_AR_ints(1);
    double **tAR_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
      calc_info_.nrio);
    
    C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,
      calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
      calc_info_.noccA*calc_info_.nvirA,&(B_p_AR[0][0]),calc_info_.nrio,
      0.0,&(tAR_p_AR[0][0]),calc_info_.nrio);
    
    write_IJKL(tAR_p_AR,PSIF_SAPT_AMPS,"Theta(2)(AR) AR",calc_info_.noccA*
      calc_info_.nvirA,calc_info_.nrio);
    free_block(B_p_AR);
  
  }

  psio_->write_entry(PSIF_SAPT_AMPS,"T2 ARAR Amplitudes",(char *) &(t2ARAR[0]),
    sizeof(double)*calc_info_.noccA*calc_info_.nvirA*calc_info_.noccA*
    calc_info_.nvirA);
  write_IJKL(thetaARAR,PSIF_SAPT_AMPS,"T2 ARAR Antisym Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);
  free(t2ARAR);
}

void SAPT::t2_bsbs(int theta)
{
  double *t2BSBS = t2_solver("T BSBS Amplitudes","Theta(BS) BS",
    PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals","SS RI Integrals",
    calc_info_.evalsB,calc_info_.noccB,calc_info_.nvirB);

  double **thetaBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);

  for (int b=0, bs=0; b < calc_info_.noccB; b++) {
  for (int s=0; s < calc_info_.nvirB; s++, bs++) {
    for (int bb=0, bbss=0; bb < calc_info_.noccB; bb++) {
    for (int ss=0; ss < calc_info_.nvirB; ss++, bbss++) {
      int bbs = bb*calc_info_.nvirB + s;
      int bss = b*calc_info_.nvirB + ss;
      int bsbbss = bs*calc_info_.noccB*calc_info_.nvirB+bbss;
      int bbsbss = bbs*calc_info_.noccB*calc_info_.nvirB+bss;
      thetaBSBS[bs][bbss] = 2.0*t2BSBS[bsbbss] - t2BSBS[bbsbss];
    }}
  }}

  if (theta) {

    double **B_p_BS = get_BS_ints(1);
    double **tBS_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
      calc_info_.nrio);
    
    C_DGEMM('N','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
      calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
      calc_info_.noccB*calc_info_.nvirB,&(B_p_BS[0][0]),calc_info_.nrio,
      0.0,&(tBS_p_BS[0][0]),calc_info_.nrio);
    
    write_IJKL(tBS_p_BS,PSIF_SAPT_AMPS,"Theta(2)(BS) BS",calc_info_.noccB*
      calc_info_.nvirB,calc_info_.nrio);
    free_block(B_p_BS);
  
  }

  psio_->write_entry(PSIF_SAPT_AMPS,"T2 BSBS Amplitudes",(char *) &(t2BSBS[0]),
    sizeof(double)*calc_info_.noccB*calc_info_.nvirB*calc_info_.noccB*
    calc_info_.nvirB);
  write_IJKL(thetaBSBS,PSIF_SAPT_AMPS,"T2 BSBS Antisym Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);
  free(t2BSBS);

}

double *SAPT::t2_solver(char *T_amps, char *theta_OV, int dfnum,
  char *OO_label, char *OV_label, char *VV_label, double *evals, int nocc,
  int nvir)
{
  double *t2ARAR = init_array(nocc*nvir*nocc*nvir);

  double **B_p_AA = get_DF_ints(dfnum,OO_label,nocc*nocc);
  double **B_p_RR = get_DF_ints(dfnum,VV_label,nvir*nvir);

  double **OVOV = block_matrix(nocc*nvir,nocc*nvir);

  for (int a=0,ar=0; a<nocc; a++) {
    for (int r=0; r<nvir; r++,ar++) {
      C_DGEMM('N','T',nocc,nvir,calc_info_.nrio,1.0,&(B_p_AA[a*nocc][0]),
        calc_info_.nrio,&(B_p_RR[r*nvir][0]),calc_info_.nrio,0.0,
        &(OVOV[ar][0]),nvir);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double *tOVOV = init_array(nocc*nvir*nocc*nvir);
  psio_->read_entry(PSIF_SAPT_AMPS,T_amps,(char *) &(tOVOV[0]),
    sizeof(double)*nocc*nvir*nocc*nvir);

  double *X = init_array(nvir);

  for(int a=0; a<nocc; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<nvir; r++) {
      int ara1 = (a*nvir+r)*nocc*nvir+a1*nvir;
      int a1ra = (a1*nvir+r)*nocc*nvir+a*nvir;
      C_DCOPY(nvir,&(tOVOV[ara1]),1,X,1);
      C_DCOPY(nvir,&(tOVOV[a1ra]),1,&(tOVOV[ara1]),1);
      C_DCOPY(nvir,X,1,&(tOVOV[a1ra]),1);
  }}}

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,nocc*nvir,-1.0,&(OVOV[0][0]),nocc*nvir,
          &(tOVOV[0]),nocc*nvir,1.0,&(t2ARAR[0]),nocc*nvir);

  for(int a=0; a<nocc; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<nvir; r++) {
      int ara1 = (a*nvir+r)*nocc*nvir+a1*nvir;
      int a1ra = (a1*nvir+r)*nocc*nvir+a*nvir;
      C_DCOPY(nvir,&(tOVOV[ara1]),1,X,1);
      C_DCOPY(nvir,&(tOVOV[a1ra]),1,&(tOVOV[ara1]),1);
      C_DCOPY(nvir,X,1,&(tOVOV[a1ra]),1);
  }}}

  for(int a=0; a<nocc; a++) {
  for(int a1=0; a1<=a; a1++) {
    for(int r=0; r<nvir; r++) {
      int ara1 = (a*nvir+r)*nocc*nvir+a1*nvir;
      int a1ra = (a1*nvir+r)*nocc*nvir+a*nvir;
      C_DCOPY(nvir,&(t2ARAR[ara1]),1,X,1);
      C_DCOPY(nvir,&(t2ARAR[a1ra]),1,&(t2ARAR[ara1]),1);
      C_DCOPY(nvir,X,1,&(t2ARAR[a1ra]),1);
  }}}

  free(X);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,nocc*nvir,-1.0,&(OVOV[0][0]),nocc*nvir,
          &(tOVOV[0]),nocc*nvir,1.0,&(t2ARAR[0]),nocc*nvir);

  free_block(OVOV);
 
  double **B_p_OV = get_DF_ints(dfnum,OV_label,nocc*nvir);
  double **T_p_OV = read_IJKL(PSIF_SAPT_AMPS,theta_OV,nocc*nvir,
    calc_info_.nrio);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,calc_info_.nri,1.0,&(B_p_OV[0][0]),
          calc_info_.nrio,&(T_p_OV[0][0]),calc_info_.nrio,1.0,&(t2ARAR[0]),
          nocc*nvir);

  free_block(B_p_OV);
  free_block(T_p_OV);

  for(int ar=0; ar<nocc*nvir; ar++) {
    for(int a1r1=0; a1r1<ar; a1r1++) {
      int ara1r1 = ar*nocc*nvir + a1r1;
      int a1r1ar = a1r1*nocc*nvir + ar;
      double tval = t2ARAR[ara1r1] + t2ARAR[a1r1ar];
      t2ARAR[a1r1ar] = tval;
      t2ARAR[ara1r1] = tval;
  }}
    
  C_DSCAL(nocc*nvir,2.0,&(t2ARAR[0]),nocc*nvir+1);

  double **xVO = block_matrix(nvir,nocc);

  for(int a=0; a<nocc; a++) {
    for(int r=0; r<nvir; r++) {
      C_DCOPY(nocc*nvir,&(tOVOV[a*nvir*nocc*nvir+r]),nvir,xVO[0],1);
      for(int a1=0; a1<nocc; a1++) {
        int aa1 = a*nocc+a1;
        int r1r = r;
        int aa1r1r = aa1*nvir*nvir+r1r;
        C_DCOPY(nvir,&(xVO[0][a1]),nocc,&(tOVOV[aa1r1r]),nvir);
  }}}

  for(int a=0; a<nocc; a++) {
    for(int r=0; r<nvir; r++) {
      C_DCOPY(nocc*nvir,&(t2ARAR[a*nvir*nocc*nvir+r]),nvir,xVO[0],1);
      for(int a1=0; a1<nocc; a1++) {
        int aa1 = a*nocc+a1;
        int r1r = r;
        int aa1r1r = aa1*nvir*nvir+r1r;
        C_DCOPY(nvir,&(xVO[0][a1]),nocc,&(t2ARAR[aa1r1r]),nvir);
  }}}

  free_block(xVO);

  B_p_AA = get_DF_ints(dfnum,OO_label,nocc*nocc);
  double **OOOO = block_matrix(nocc*nocc,nocc*nocc);

  for (int a=0,aaa=0; a<nocc; a++) {
    for (int aa=0; aa<nocc; aa++,aaa++) {
      C_DGEMM('N','T',nocc,nocc,calc_info_.nrio,1.0,&(B_p_AA[a*nocc][0]),
        calc_info_.nrio,&(B_p_AA[aa*nocc][0]),calc_info_.nrio,0.0,
        &(OOOO[aaa][0]),nocc);
  }}

  free_block(B_p_AA);

  C_DGEMM('N','N',nocc*nocc,nvir*nvir,nocc*nocc,1.0,&(OOOO[0][0]),nocc*nocc,
          &(tOVOV[0]),nvir*nvir,1.0,&(t2ARAR[0]),nvir*nvir);

  free_block(OOOO);

  B_p_RR = get_DF_ints(dfnum,VV_label,nvir*nvir);
  double **X_RR = block_matrix(nvir*nvir,nvir);

  for (int r=0; r < nvir; r++) {
    C_DGEMM('N','T',nvir,nvir*nvir,calc_info_.nrio,1.0,&(B_p_RR[r*nvir][0]),
            calc_info_.nrio,&(B_p_RR[0][0]),calc_info_.nrio,0.0,&(X_RR[0][0]),
            nvir*nvir);
    C_DGEMM('N','T',nocc*nocc,nvir*nvir,nvir,1.0,&(tOVOV[r*nvir]),
            nvir*nvir,&(X_RR[0][0]),nvir,1.0,&(t2ARAR[0]),nvir*nvir);
  }

  free(tOVOV);
  free_block(B_p_RR);
  free_block(X_RR);

  double **xOV = block_matrix(nocc,nvir);

  for(int a=0; a<nocc; a++) {
    for(int r=0; r<nvir; r++) {
      C_DCOPY(nocc*nvir,&(t2ARAR[a*nocc*nvir*nvir+r]),nvir,xOV[0],1);
      for(int a1=0; a1<nocc; a1++) {
        for(int r1=0; r1<nvir; r1++) {
          int ar1 = a*nvir+r1;
          int a1r = a1*nvir+r;
          int ar1a1r = ar1*nocc*nvir+a1r;
          double denom = evals[a]+evals[a1]-evals[r+nocc]-evals[r1+nocc];
          t2ARAR[ar1a1r] = xOV[a1][r1]/denom;
  }}}}

  free_block(xOV);

  return(t2ARAR);
}

void SAPT::g_arar()
{
  double **B_p_AR = get_AR_ints(0);
  double **gARAR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccA*calc_info_.nvirA);
  
  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio,2.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(B_p_AR[0][0]),calc_info_.nrio,0.0,&(gARAR[0][0]),calc_info_.noccA*
    calc_info_.nvirA);

  free_block(B_p_AR);

  double **B_p_AA = get_AA_ints(0);
  double **B_p_RR = get_RR_ints(0);

  for(int a=0,ar=0; a<calc_info_.noccA; a++) {
    for(int r=0; r<calc_info_.nvirA; r++,ar++) {
      C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.nrio,-1.0,
        &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,
        &(B_p_RR[r*calc_info_.nvirA][0]),calc_info_.nrio,1.0,
        &(gARAR[ar][0]),calc_info_.nvirA);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  double **gtAR_ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.noccA*calc_info_.nvirA,1.0,&(gARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA,&(tARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,0.0,&(gtAR_ARBS[0][0]),calc_info_.noccB*calc_info_.nvirB);

  free_block(gARAR);
  free_block(tARBS);

  write_IJKL(gtAR_ARBS,PSIF_SAPT_AMPS,"GxT(AR) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
}

void SAPT::g_bsbs()
{
  double **B_p_BS = get_BS_ints(0);
  double **gBSBS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.noccB*calc_info_.nvirB);
  
  C_DGEMM('N','T',calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio,2.0,&(B_p_BS[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(gBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_BS);

  double **B_p_BB = get_BB_ints(0);
  double **B_p_SS = get_SS_ints(0);

  for(int b=0,bs=0; b<calc_info_.noccB; b++) {
    for(int s=0; s<calc_info_.nvirB; s++,bs++) {
      C_DGEMM('N','T',calc_info_.noccB,calc_info_.nvirB,calc_info_.nrio,-1.0,
        &(B_p_BB[b*calc_info_.noccB][0]),calc_info_.nrio,
        &(B_p_SS[s*calc_info_.nvirB][0]),calc_info_.nrio,1.0,
        &(gBSBS[bs][0]),calc_info_.nvirB);
  }}

  free_block(B_p_BB);
  free_block(B_p_SS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  double **gtBS_ARBS = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB,1.0,&(tARBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB,&(gBSBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB,0.0,&(gtBS_ARBS[0][0]),calc_info_.noccB*calc_info_.nvirB);

  free_block(gBSBS);
  free_block(tARBS);

  write_IJKL(gtBS_ARBS,PSIF_SAPT_AMPS,"GxT(BS) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
}

void SAPT::natural_orbitalify(char *OO_opdm, char *VV_opdm, double *evals, 
  double **scfvec, int occ, int vir, char monomer)
{
  double **P = block_matrix(calc_info_.nmo,calc_info_.nmo);
  double **xAA = read_IJKL(PSIF_SAPT_AMPS,OO_opdm,occ,occ);
  double **xRR = read_IJKL(PSIF_SAPT_AMPS,VV_opdm,vir,vir);

  for (int i=0; i<occ; i++) // Put in electrons
    P[i][i] = 2.0;

  for (int i=0; i<occ; i++)
    C_DAXPY(occ,-2.0,xAA[i],1,P[i],1);

  for (int i=0; i<vir; i++)
    C_DAXPY(vir,2.0,xRR[i],1,&(P[i+occ][occ]),1);

  free_block(xAA);
  free_block(xRR);

  double *occnum = init_array(calc_info_.nmo);
  double **nat_orbs_MO = block_matrix(calc_info_.nmo,calc_info_.nmo);

  sq_rsp(calc_info_.nmo,calc_info_.nmo,P,occnum,3,nat_orbs_MO,1.0e-14);

  int num_no_vir = 0;

  double virt_elec = 0.0;
  double temp_elec = 0.0;

  for (int i=occ; i<calc_info_.nmo; i++) {
    if (occnum[i] > params_.occ_cutoff) {
      num_no_vir++;
    }
    else break;
  }

  if (params_.print) 
    fprintf(outfile,"  Monomer %c: %d virtual orbitals dropped\n",monomer,
          vir-num_no_vir);

  double **Fock_MO = block_matrix(calc_info_.nmo,calc_info_.nmo);

  for (int i=0; i<calc_info_.nmo; i++) {
    Fock_MO[i][i] = evals[i];
  }

  double **tempmat = block_matrix(occ+num_no_vir,calc_info_.nmo);
  double **Fock_NO = block_matrix(occ+num_no_vir,occ+num_no_vir);

  C_DGEMM('T','N',occ+num_no_vir,calc_info_.nmo,calc_info_.nmo,1.0,
          &(nat_orbs_MO[0][0]),calc_info_.nmo,&(Fock_MO[0][0]),calc_info_.nmo,
          0.0,&(tempmat[0][0]),calc_info_.nmo);
  C_DGEMM('N','N',occ+num_no_vir,occ+num_no_vir,calc_info_.nmo,1.0,
          &(tempmat[0][0]),calc_info_.nmo,&(nat_orbs_MO[0][0]),calc_info_.nmo,
          0.0,&(Fock_NO[0][0]),occ+num_no_vir);

  double *epsilon = init_array(occ+num_no_vir);
  double **X = block_matrix(occ+num_no_vir,occ+num_no_vir);
  sq_rsp(occ+num_no_vir,occ+num_no_vir,Fock_NO,epsilon,1,X,1.0e-14);

  double **MO_MVO = block_matrix(calc_info_.nmo,occ+num_no_vir);

  C_DGEMM('N','N',calc_info_.nmo,occ+num_no_vir,occ+num_no_vir,1.0,
          &(nat_orbs_MO[0][0]),calc_info_.nmo,&(X[0][0]),occ+num_no_vir,
          0.0,&(MO_MVO[0][0]),occ+num_no_vir);

  if (monomer == 'A') {
    no_info_.CA = block_matrix(calc_info_.nmo,occ+num_no_vir);
    no_info_.evalsA = init_array(occ+num_no_vir);
    no_info_.nvirA = num_no_vir;

    C_DCOPY(calc_info_.nmo*(occ+num_no_vir),&(MO_MVO[0][0]),1,
      &(no_info_.CA[0][0]),1);
    C_DCOPY(occ+num_no_vir,epsilon,1,no_info_.evalsA,1);
  }

  if (monomer == 'B') {
    no_info_.CB = block_matrix(calc_info_.nmo,occ+num_no_vir);
    no_info_.evalsB = init_array(occ+num_no_vir);
    no_info_.nvirB = num_no_vir;

    C_DCOPY(calc_info_.nmo*(occ+num_no_vir),&(MO_MVO[0][0]),1,
      &(no_info_.CB[0][0]),1);
    C_DCOPY(occ+num_no_vir,epsilon,1,no_info_.evalsB,1);
  }

  free(epsilon);
  free(occnum);
  free_block(P);
  free_block(nat_orbs_MO);
  free_block(tempmat);
  free_block(Fock_MO);
  free_block(Fock_NO);
  free_block(X);
  free_block(MO_MVO);

}

}}