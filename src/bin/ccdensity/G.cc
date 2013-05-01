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
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <string>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void G_build(void) {
  dpdbuf4 LIJAB, Lijab, LiJaB, LIjAb, LijAB, LIJab;
  dpdbuf4 tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab;
  dpdfile2 GAE, Gae, GMI, Gmi;
  int G_irr;
  G_irr = params.G_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&GMI, PSIF_CC_GLG, G_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, PSIF_CC_GLG, G_irr, 0, 0, "Gmi");

    /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
    dpd_buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&LIJAB, PSIF_CC_GLG, G_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&tIJAB, &LIJAB, &GMI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&LIJAB);  

    /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
    dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&LIjAb, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&LIjAb);

    /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
    dpd_buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&Lijab, PSIF_CC_GLG, G_irr, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&tijab, &Lijab, &Gmi, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&Lijab); 

    /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
    dpd_buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&LiJaB, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&tiJaB, &LiJaB, &Gmi, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&tiJaB);
    dpd_buf4_close(&LiJaB);

    dpd_file2_close(&Gmi);
    dpd_file2_close(&GMI);


  
    dpd_file2_init(&GAE, PSIF_CC_GLG, G_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, PSIF_CC_GLG, G_irr, 1, 1, "Gae");

    /* T2(IJ,AB) * L2(IJ,EB) --> G(A,E) */
    dpd_buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&LIJAB, PSIF_CC_GLG, G_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&LIJAB, &tIJAB, &GAE, 2, 2, -1.0, 0.0);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&LIJAB);

    /* T2(Ij,Ab) * L2(Ij,Eb) --> G(A,E) */
    dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&LIjAb, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&LIjAb, &tIjAb, &GAE, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&LIjAb);

    /* T2(ij,ab) * L2(ij,eb) --> G(a,e) */
    dpd_buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&Lijab, PSIF_CC_GLG, G_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&Lijab, &tijab, &Gae, 2, 2, -1.0, 0.0);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&Lijab);

    /* T2(iJ,aB) * L2(iJ,eB) --> G(a,e) */
    dpd_buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&LiJaB, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&LiJaB, &tiJaB, &Gae, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&tiJaB);
    dpd_buf4_close(&LiJaB);

    dpd_file2_close(&GAE);
    dpd_file2_close(&Gae);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&GMI, PSIF_CC_GLG, G_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, PSIF_CC_GLG, G_irr, 2, 2, "Gmi");

    /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
    dpd_buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&LIJAB, PSIF_CC_GLG, G_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&tIJAB, &LIJAB, &GMI, 0, 0, 1, 0);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&LIJAB);  

    /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
    dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&LIjAb, PSIF_CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1, 1);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&LIjAb);

    /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
    dpd_buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&Lijab, PSIF_CC_GLG, G_irr, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&tijab, &Lijab, &Gmi, 0, 0, 1, 0);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&Lijab); 

    /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
    dpd_buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&LiJaB, PSIF_CC_GLG, G_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&tiJaB, &LiJaB, &Gmi, 0, 0, 1, 1);
    dpd_buf4_close(&tiJaB);
    dpd_buf4_close(&LiJaB);

    dpd_file2_close(&Gmi);
    dpd_file2_close(&GMI);


  
    dpd_file2_init(&GAE, PSIF_CC_GLG, G_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, PSIF_CC_GLG, G_irr, 3, 3, "Gae");

    /* T2(JI,BA) * L2(JI,BE) --> G(A,E) */
    dpd_buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&LIJAB, PSIF_CC_GLG, G_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&LIJAB, &tIJAB, &GAE, 3, 3, -1, 0);
    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&LIJAB);

    /* T2(jI,bA) * L2(jI,bE) --> G(A,E) */
    dpd_buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&LIjAb, PSIF_CC_GLG, G_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&LIjAb, &tIjAb, &GAE, 3, 3, -1, 1);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&LIjAb);

    /* T2(ji,ba) * L2(ji,be) --> G(a,e) */
    dpd_buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&Lijab, PSIF_CC_GLG, G_irr, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&Lijab, &tijab, &Gae, 3, 3, -1, 0);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&Lijab);

    /* T2(Ji,Ba) * L2(Ji,Be) --> G(a,e) */
    dpd_buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&LiJaB, PSIF_CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&LiJaB, &tiJaB, &Gae, 3, 3, -1, 1);
    dpd_buf4_close(&tiJaB);
    dpd_buf4_close(&LiJaB);

    dpd_file2_close(&GAE);
    dpd_file2_close(&Gae);

  }

  return;
}



}} // namespace psi::ccdensity
