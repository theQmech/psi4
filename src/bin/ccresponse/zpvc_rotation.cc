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

// uncomment the ones needed as and when required

namespace psi { namespace ccresponse {

void print_tensor_der(boost::shared_ptr<OutFile> myfile, std::vector<SharedMatrix> my_tensor_list);

void zpvc_rotation(boost::shared_ptr<Molecule> molecule,
    Options &options,
    double step,
    std::vector<SharedMatrix> G,
    SharedMatrix G0)
{
    double mstep = options.get_double("DISP_SIZE");
    //-> This is troublesome, need to decide if option should be set as global or local-"FINDIF"
    outfile->Printf("STEPSIZE from Options Object = %lf bohr\n\n",mstep);
    outfile->Printf("STEPSIZE passed from driver and used = %lf bohr\n\n",step);

    int print = options.get_int("PRINT");
    outfile->Printf("Print Level = %d\n", print);
    int i,j,k;
    int a,b,c,d;

    int count, nomega;
    double omega;
    count = options["OMEGA"].size();
    if(count == 0) { // Assume 0.0 E_h for field energy
      nomega = 1;
      params.omega = init_array(1);
      params.omega[0] = 0.0;
      omega = params.omega[0];
    }
    else if(count == 1) { // Assume E_h for field energy and read value
      params.nomega = 1;
      params.omega = init_array(1);
      params.omega[0] = options["OMEGA"][0].to_double();
      omega = params.omega[0];
    }
    else if(count == 2) {
      params.nomega = count-1;
      params.omega = init_array(params.nomega);
      std::string units = options["OMEGA"][count-1].to_string();
      for(i=0; i < (count-1); i++) {
        params.omega[i] = options["OMEGA"][i].to_double();
        if(units == "HZ" || units == "Hz" || units == "hz")
          params.omega[i] *= pc_h / pc_hartree2J;
        else if(units == "AU" || units == "Au" || units == "au") continue; // do nothing
        else if(units == "NM" || units == "nm")
          params.omega[i] = (pc_c*pc_h*1e9)/(params.omega[i]*pc_hartree2J);
        else if(units == "EV" || units == "ev" || units == "eV")
          params.omega[i] /= pc_hartree2ev;
        else
          throw PsiException("Error in unit for input field frequencies, should be au, Hz, nm, or eV", __FILE__,__LINE__);
      }
      omega = params.omega[0];
    }
    else if (count > 2) {
      throw PsiException("ROA Scattering only working for one wavelength at a time", __FILE__,__LINE__);
    }
    // Print the Wavelength
    outfile->Printf("\t Wavelength (in au): = %20.12f\n\n", omega);

  int natom = molecule->natom();
  double alpha_eq = tensor_mean(G0);
  std::vector<double> alpha_diffs;

  // Iterate through list of G' tensors, set elements of alpha_difs to the
  // numerator of the finite differences
  for(std::vector<SharedMatrix>::iterator G_it= G.begin(); G_it!=G.end(); ++G_it)
  {
    double alpha_p = tensor_mean(*G_it);
    ++G_it;
    double alpha_m = tensor_mean(*G_it);
    alpha_diffs.push_back(alpha_p - 2*alpha_eq + alpha_m);
  }
  // construct a psi4 vector of 2nd derivatives //
  SharedVector d2_alpha_dx2(new Vector("Cartesian 2nd Derivatives [alpha]",3*natom));
  d2_alpha_dx2->set(alpha_diffs.data());
  //divide each by step^2
  d2_alpha_dx2->scale(1/(step*step));
  std::vector<char> int2xyz={'x','y','z'};
  outfile->Printf("\n     Cartesian 2nd Derivatives\n");
  outfile->Printf(   "===================================\n");
  for(int a = 0; a < natom*3; a++){
    outfile->Printf( "d^2[alpha]/d%d_{%c}^{2}:   %10.5f\n",
        int2xyz[a%3],
        a/3,
        d2_alpha_dx2->get(a));
  }


    SharedMatrix geom(new Matrix(natom,3));

    // Reading in the Hessian //
    FILE* hessian;
    FILE* dipole_moment;
    hessian=fopen("file15.dat","r");
    SharedMatrix F(new Matrix(natom*3,natom*3));
    double Fval;
    SharedMatrix M(new Matrix(natom*3,natom*3));
    for(i=0; i < (3*natom); i++)
    {
      for(j=0; j < (3*natom); j++)
      {
        int statusvalue=fscanf(hessian,"%lf",&F->pointer()[i][j]);
      }
    }
    fclose(hessian);

    geom->copy(molecule->geometry());
    SharedMatrix geom_orig(new Matrix(natom,3));
    geom_orig->copy(geom);

	// Translate Molecule to Center of Mass //
	if(print >= 1)  {
      outfile->Printf("\tInput coordinates (bohr):\n");
	  molecule->geometry().print();
    }
	molecule->move_to_com();
	geom_orig->copy(molecule->geometry());


    // Mass-weighting the co-ordinates //
    //double massi[natom];
    outfile->Printf("\tAtomic Masses for Raman Computation:\n");
    for(i=0; i < natom; i++) outfile->Printf("\t%d %12.8f\n", i, molecule->mass(i));
    outfile->Printf("\n");

    for(i=0; i<natom; i++)
    {
        //massi[i] = molecule->mass(i);
        for(j=0; j<3; j++)
        {
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

/*
    for(i=0; i < 3; i++)
    {
      for(j=0; j <= i; j++)
      {
        if(i==j)
          for(k=0; k < natom; k++)
            I->add(i,j,(geom->get(k,(i+1)%3)*geom->get(k,(i+1)%3) + geom->get(k,(i+2)%3)*geom->get(k,(i+2)%3)));
        else
        {
          for(k=0; k < natom; k++)
            I->add(i,j,-1.0 * (geom->get(k,i)*geom->get(k,j)));
        I->set(j,i,I->get(i,j));
        }
      }
    }
    //I->print();
*/

    // Diagonalizing the inertia tensor //

    SharedMatrix Ievecs(new Matrix("Inertia Eigenvectors",3,3));
    SharedVector Ievals(new Vector("Inertia Eigenvalues",3));

    I->diagonalize(Ievecs,Ievals);
	//rs(3,3,I->pointer(),Ievals->pointer(),1,Ievecs->pointer(),1e-12);

    // Constructing I-inverse matrix //

    SharedMatrix Iinv(new Matrix("I-Inverse",3,3));
    SharedMatrix Itmp(new Matrix(3,3));

    Iinv->zero();
    for(i=0;i<3;i++)
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
    int icart,jcart,iatom,jatom;
    double imass,jmass,total_mass;

    total_mass=0.0;
    for(i=0;i<natom;i++)
    {
      //total_mass+=an2mass[(int)zvals[i]];
      total_mass+=molecule->mass(i);
    }

    for(i=0; i < natom*3; i++)
    {
      icart = i % 3;
      iatom = i/3;
      //imass = an2mass[(int)zvals[iatom]];
      imass = molecule->mass(iatom);

      P->set(i,i,1.0);

      for(j=0; j < natom*3; j++)
      {
        jcart = j % 3;
        jatom = j/3;
        //jmass = an2mass[(int)zvals[jatom]];
        jmass = molecule->mass(jatom);

        P->add(i,j,-1.0*sqrt(imass*jmass)/total_mass*(icart==jcart));

        for(a=0; a < 3; a++)
          for(b=0; b < 3; b++)
            for(c=0; c < 3; c++)
              for(d=0; d < 3; d++)
              {
                P->add(i,j,-1.0*levi(a,b,icart)*geom->get(iatom,b)*Iinv->get(a,c)*levi(c,d,jcart)*geom->get(jatom,d));
              }

      }
    }

    // Generate mass-weighted Hessian matrix [Eh/(bohr^2 amu)] //
    M->zero();
    SharedMatrix T(new Matrix(natom*3,natom*3));
    for(i=0; i < natom; i++)
    {
      for(j=0; j < 3; j++)
      {
        //M->set((i*3+j),(i*3+j),1/sqrt(an2mass[(int)zvals[i]]/_au2amu));
        M->set((i*3+j),(i*3+j),1/sqrt((molecule->mass(i))/pc_au2amu));
      }
    }
    //printf("Mass-Weighting Matrix (for Hessian):\n");
	//M->print(stdout);

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

    // Diagonalize projected mass-weighted Hessian //
    SharedMatrix Fevecs(new Matrix(3*natom,3*natom));
    SharedVector Fevals(new Vector("Feigenval",3*natom));
    SharedMatrix Lx(new Matrix("Normal Transform Matrix",3*natom,3*natom));
    SharedVector redmass(new Vector("ReducedMass",3*natom));
    double norm=0.0;

    F->diagonalize(Fevecs,Fevals);
	//rs(3*natom,3*natom,F->pointer(),Fevals->pointer(),1,Fevecs->pointer(),1e-12);
	if(print >= 3)  {
	  Fevals->print();
	  Fevecs->print();
    }

    Lx->gemm(0,0,1.0,M,Fevecs,0.0);
	if(print >= 2)  {
	  Lx->print();
    }

    for(i=0; i < 3*natom; i++)
    {
       norm = 0.0;
       for(j=0; j < 3*natom; j++)
       {
         norm += Lx->get(j,i)*Lx->get(j,i)/pc_au2amu;
       }
       if(norm > 1e-3)
       {
       redmass->set(i,1.0/norm);
       }
    }

  redmass->print();



}

}} // namespace psi::ccresponse


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
