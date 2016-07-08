/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup ccresponse
    \brief Compute the three tensors needed for Raman Optical Activity.

    ROA requires the following polarizability tensors:
      (1) electric-dipole/electric-dipole;
      (2) electric-dipole/electric-quadrupole; and
      (3) electric-dipole/magnetic-dipole.

  -TDC, August 2009
*/

// Remember to edit this comment later on

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <vector>
#include <psi4-dec.h>
#include <libmints/mints.h>
#include <physconst.h>
#include <liboptions/liboptions.h>
#include "libparallel/ParallelPrinter.h"

//#include "oldphysconst.h"
//#include "mass.h"

using namespace boost;
using namespace psi;

 int levi(int a, int b, int c);
 double tensor_mean(SharedMatrix alpha);
 double beta_alpha2(SharedMatrix alpha);
 double beta_G2(SharedMatrix alpha, SharedMatrix G);
 double beta_A2(SharedMatrix alpha, double ***A, double omega);

 double raman_linear(double alpha, double beta2);
 double depolar_linear(double alpha, double beta2);
 double raman_circular(double alpha, double beta2);
 double depolar_circular(double alpha, double beta2);
 void rs(int nm, int n, double **array, double *e_vals, int matz,
         double **e_vecs, double toler);

double to_delta_x_sq(double nu, double theta);

// uncomment the ones needed as and when required

namespace psi { namespace ccresponse {

void print_tensor_der(boost::shared_ptr<OutFile> myfile, std::vector<SharedMatrix> my_tensor_list);

void rotation_vibave_cartesian(boost::shared_ptr<Molecule> molecule,
    Options &options,
    SharedMatrix d2_alpha_dx2)
{
  // printlevel option to control verbose output //
  int print = options.get_int("PRINT");
  // Get the geometry at eq point//
  int natom = molecule->natom();
  SharedMatrix geom(new Matrix(natom,3));
  geom->copy(molecule->geometry());

  // Reading in the Hessian //
  FILE* hessian;
  FILE* dipole_moment;
  hessian=fopen("file15.dat","r");
  SharedMatrix F(new Matrix(natom*3,natom*3));
  double Fval;
  SharedMatrix M(new Matrix(natom*3,natom*3));
  for(int i=0; i < (3*natom); i++)
  {
    for(int j=0; j < (3*natom); j++)
    {
      int statusvalue=fscanf(hessian,"%lf",&F->pointer()[i][j]);
    }
  }
  fclose(hessian);

  //make a copy of the original geometry//
  SharedMatrix geom_orig(new Matrix(natom,3));
  geom_orig->copy(geom);

  // Translate Molecule to Center of Mass //
  if(print >= 1)  {
      outfile->Printf("\tInput coordinates (bohr):\n");
    molecule->geometry().print();
    }
  molecule->move_to_com();
  geom_orig->copy(molecule->geometry()); //


    // Mass-weighting the coordinates //
  outfile->Printf("\tAtomic Masses for zpvc to optical rotation:\n");
  for(int i=0; i < natom; i++) outfile->Printf("\t%d %12.8f\n", i, molecule->mass(i));
  outfile->Printf("\n");

  for(int i=0; i<natom; i++){
    //massi[i] = molecule->mass(i);
    for(int j=0; j<3; j++){
      //geom->set(i,j, geom_orig->get(i,j) * sqrt(massi[i]));
      geom->set(i,j, geom_orig->get(i,j) * sqrt(molecule->mass(i)));
    }
  }
  if(print >= 1)  {
    outfile->Printf("\tMass-Weighted coordinates relative to the center of mass:\n");
    geom->print();
  }
    //delete[] massi;

    //Generating the inertia tensor

    SharedMatrix I(new Matrix(3,3));
    I->copy(molecule->inertia_tensor());
    //I = molecule->inertia_tensor();
  if(print >= 2)  {
    outfile->Printf("\tMoment of Inertia Tensor:\n");
    I->print();
  }

  // Diagonalizing the inertia tensor //

  SharedMatrix Ievecs(new Matrix("Inertia Eigenvectors",3,3));
  SharedVector Ievals(new Vector("Inertia Eigenvalues",3));

  I->diagonalize(Ievecs,Ievals);

  // Constructing I-inverse matrix //

  SharedMatrix Iinv(new Matrix("I-Inverse",3,3));
  SharedMatrix Itmp(new Matrix(3,3));

  Iinv->zero();
  for(int i=0;i<3;i++)
  {
    Iinv->set(i,i,(1.0/Ievals->get(i)));
  }

  Itmp->gemm(0,1,1.0,Iinv,Ievecs,0.0);
  Iinv->gemm(0,0,1.0,Ievecs,Itmp,0.0);

  if(print >= 2)  {
    outfile->Printf("\tInertia Tensor EigenData\n\n");
    Ievecs->print();
    Ievals->print();
    Iinv->print();
  }

  // Generating the 6 pure rotation and translation vectors //
  SharedMatrix P(new Matrix(natom*3,natom*3));

  double total_mass=0.0;
  for(int i=0;i<natom;i++) {
    total_mass+=molecule->mass(i);
  }

  for(int i=0; i < natom*3; i++) {
    int icart = i % 3;
    int iatom = i/3;
    int imass = molecule->mass(iatom);

    P->set(i,i,1.0);

    for(int j=0; j < natom*3; j++) {
      int jcart = j % 3;
      int jatom = j/3;
      int jmass = molecule->mass(jatom);

      P->add(i,j,-1.0*sqrt(imass*jmass)/total_mass*(icart==jcart));

      for(int a=0; a < 3; a++){
        for(int b=0; b < 3; b++){
          for(int c=0; c < 3; c++){
            for(int d=0; d < 3; d++){
              P->add(i,j,-1.0*levi(a,b,icart)*geom->get(iatom,b)*Iinv->get(a,c)*levi(c,d,jcart)*geom->get(jatom,d));
            }
          }
        }
      }
    }
  }

    // Generate mass-weighted Hessian matrix [Eh/(bohr^2 amu)] //
  M->zero();
  SharedMatrix T(new Matrix(natom*3,natom*3));
  for(int i=0; i < natom; i++){
    for(int j=0; j < 3; j++){
      M->set((i*3+j),(i*3+j),1/sqrt((molecule->mass(i))/pc_au2amu));
    }
  }

  T->gemm(0,0,1.0,M,F,0.0);
  F->gemm(0,0,1.0,T,M,0.0);

  // Project out rotational and translational degrees of freedom from mass-weighted Hessian //
  T->zero();
  T->gemm(0,0,1.0,F,P,0.0);
  F->gemm(0,0,1.0,P,T,0.0);
  if(print >= 2)  {
    outfile->Printf("\tProjected, Mass-Weighted Hessian:\n");
    F->print();
  }

  SharedMatrix Fevecs(new Matrix(3*natom,3*natom));
  SharedVector Fevals(new Vector("Feigenval",3*natom));
  // Transformer from 3n Cartesian basis to basis of normal modes
  // with translations-rotations projected out also called the interal
  // coordinate system
  SharedMatrix Lx(new Matrix("Normal Internal Transform Matrix",3*natom,3*natom));
  //The reduced mass each normal mode
  SharedVector redmass(new Vector("ReducedMass",3*natom));
  double norm=0.0;

  // Diagonalize projected mass-weighted Hessian //
  F->diagonalize(Fevecs,Fevals);
  if(print >= 3){
    Fevals->print();
    Fevecs->print();
  }

  // Mass Weight the normal modes to get the transformer //
  Lx->gemm(0,0,1.0,M,Fevecs,0.0);
  if(print >= 2){
    Lx->print();
  }

  // Normalize the reduced masses for each mode//
  for(int i=0; i < 3*natom; i++) {
    norm = 0.0;
    for(int j=0; j < 3*natom; j++){
      norm += Lx->get(j,i)*Lx->get(j,i)/pc_au2amu;
    }
    if(norm > 1e-3) {
      redmass->set(i,1.0/norm);
    }
  }

  redmass->print();

  // From here use Lx to transform the 2nd derivatives to the normal //
  // Internal coordinate system the d^2[alpha]/dQ^2 from the papers  //
  // then compute the corrections using the expressions  from those sources //

  SharedMatrix deriv_normal(new Matrix("Derivatives in Normal Mode", natom*3,natom*3));
  deriv_normal->gemm(0,0,1.0,d2_alpha_dx2,Lx,0.0);
  deriv_normal->gemm(0,0,1.0,Lx->transpose(),deriv_normal,0.0);

  //=================================================
  outfile->Printf("Double derivative in cartesian coordinates\n");
  d2_alpha_dx2->print();
  //=================================================

  //=================================================
  outfile->Printf("Double derivative in normal coordinates\n");
  deriv_normal->print();
  //=================================================

  double km_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg * pc_au2amu);
  double m_convert = 1.0/(2.0 * pc_pi * pc_c);

  SharedVector freq(new Vector("Freq",(3*natom)));
  for (int i=0; i<3*natom; ++i){
    if(Fevals->get(i) < 0.0)
      freq->set(i, m_convert*sqrt(-km_convert*Fevals->get(i)));
    else
      freq->set(i, m_convert*sqrt(km_convert*Fevals->get(i)));
  }

  outfile->Printf("\n\t     Harmonic Freq.   Red. Mass \n");
  outfile->Printf("\t        (m-1)            (amu)    \n");
  outfile->Printf("\t---------------------------------------------------------------------------------------------------\n");
  for(int i=0; i < 3*natom; ++i)
  {
    if(Fevals->get(i) < 0.0)
      outfile->Printf("\t  %3d  %10.3fi    %7.4f\n",
       i, freq->get(i), redmass->get(i));
    else
      outfile->Printf("\t  %3d  %10.3f     %7.4f\n",
       i, freq->get(i), redmass->get(i));
  }
  outfile->Printf("\t---------------------------------------------------------------------------------------------------\n");

  double theta = 300.0;//???

  //firest six frequencies are useless
  SharedVector delta_x_sq(new Vector(3*natom));
  for (int i=6; i<3*natom; ++i){
    delta_x_sq->set(i, to_delta_x_sq(freq->get(i), theta));///???
  }

  //=================================================
  outfile->Printf("delta_x_sq\n");
  delta_x_sq->print();
  //=================================================

  double correction = 0.0;
  for (int i=6; i<3*natom; ++i){
    correction += deriv_normal->get(i,i)*delta_x_sq->get(i);
  }
  correction /= 2.0;

  outfile->Printf("\n\nFinal Result here:\n%15.12f\n", correction);
  printf("\n\nFinal Result here:\n%15.12f\n", correction);

}

}} // namespace psi::ccresponse

double to_delta_x_sq(double nu, double theta){
  double x1 = (pc_h*pc_c*nu)/(2.0*pc_kb*theta);
  double x2 = pc_h/(8*pc_pi*pc_pi*pc_c*nu * tanh(x1)); //???

  return x2;
}

/* The Levi-Civitas evaluator */

//int levi(int a, int b, int c)
//{
//  int val=0;
//  int x=0, y=1, z=2;
//
//  if(a==x && b==y && c==z) val=1;
//  else if(a==y && b==z && c==x) val=1;
//  else if(a==z && b==x && c==y) val=1;
//  else if(a==x && b==z && c==y) val=-1;
//  else if(a==y && b==x && c==z) val=-1;
//  else if(a==z && b==y && c==x) val=-1;
//  else val=0;
//
//  return val;
//}
//
///* compute mean of a property tensor: alpha = 1/3 alpha_ii */
//double tensor_mean(SharedMatrix alpha)
//{
//  double mean=0.0;
//  int i;
//  for(i=0; i < 3; i++)
//    mean += alpha->get(i,i);
//  mean /= 3.0;
//  return mean;
//}
//
///* compute beta(alpha)^2 = 1/2 [ 3 * alpha_ij*alpha_ij - alpha_ii*alpha_jj */
//double beta_alpha2(SharedMatrix alpha)
//{
//  double value = 0.0;
//  int i,j;
//  for(i=0; i < 3; i++)
//    for(j=0; j < 3; j++)
//      value += 0.5*(3.0*alpha->get(i,j)*alpha->get(i,j) - alpha->get(i,i)*alpha->get(j,j));
//
//  return value;
//}
//
///* compute beta(G')^2 = 1/2[3 * alpha_ij*G'_ij - alpha_ii*G'_jj */
//double beta_G2(SharedMatrix alpha, SharedMatrix G)
//{
//  double value = 0.0;
//  int i,j;
//  for(i=0; i < 3; i++)
//    for(j=0; j < 3; j++)
//      value += 0.5*(3.0*alpha->get(i,j)*G->get(i,j) - alpha->get(i,i)*G->get(j,j));
//  return value;
//}
//
///* compute beta(A)^2 = 1/2 omega * alpha_ij epsilon_ikl * A_klj */
//double beta_A2(SharedMatrix alpha, double ***A, double omega)
//{
//  double value=0.0;
//  int i,j,k,l;
//  for(i=0; i < 3; i++)
//    for(j=0; j < 3; j++)
//      for(k=0; k < 3; k++)
//        for(l=0; l < 3; l++)
//          value += 0.5 * omega * alpha->get(i,j) * levi(i,k,l) * A[k][l][j];
//
//  return value;
//}
//
///* compute Raman intensity for linearly polarized radiation:
//    A = 45 (alpha^2) + 4 (beta^2)
//*/
//double raman_linear(double alpha, double beta2)
//{
//  double value = 0.0;
//  value = 45.0 * alpha * alpha + 4.0 * beta2;
//  return value;
//}
//
///* compute Raman depolarization ratio for 90-degree scattering of linearly
//   polarized radiation:
//
//  ratio = [ 3 * beta(alpha)^2)/(45 * alpha^2 + 4 * beta(alpha)^2) ]
//*/
//double depolar_linear(double alpha, double beta2)
//{
//  double numer, denom;
//
//  numer = 3.0 * beta2;
//  denom = (45.0 * alpha * alpha) + (4.0 * beta2);
//
//  if(denom > 1e-6) return numer/denom;
//  else return 0.0;
//}
//
///* compute Raman intensity for circularly polarized radiation:
//    A = 45 (alpha^2) + 7 (beta^2);
//*/
//double raman_circular(double alpha, double beta2)
//{
//  double value = 0.0;
//  value = 45.0 * alpha * alpha + 7.0 * beta2;
//  return value;
//}
//
//// compute Raman depolarization ratio for 90-degree scattering of circularly polarized radiation://
//// ratio = [ 6 * beta(alpha)^2)/(45 * alpha^2 + 7 * beta(alpha)^2) ] //
//
//double depolar_circular(double alpha, double beta2)
//{
//  double numer, denom;
//
//  numer = 6.0 * beta2;
//  denom = (45.0 * alpha * alpha) + (7.0 * beta2);
//
//  if(denom > 1e-6) return numer/denom;
//  else return 0.0;
//}
