/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_Iij(void);
void uhf_sf_Iij(void);

void Iij(void)
{
  if(params.ref == 0) rhf_sf_Iij();
  else if(params.ref == 2) uhf_sf_Iij();
}

void rhf_sf_Iij(void)
{
  dpdfile2 I, F, D;
  dpdbuf4 G, Aints, Dints, Cints, Eints;

  /* I'IJ <-- sum_K fIK (DJK + DKJ) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  dpd_file2_scm(&I, 0.0);

  dpd_file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'IJ <-- 2 fIJ */
  dpd_file2_axpy(&F, &I, 2.0, 0);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'ij <-- sum_k fik (Djk + Dkj) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'ij <-- 2 fij */
  dpd_file2_axpy(&F, &I, 2.0, 0);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);
 
  /* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Aints, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);
 
  /* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
    
  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);
    
  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Cints, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
  dpd_file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

}

void uhf_sf_Iij(void)
{

}

}} /* End namespaces */