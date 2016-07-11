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

/*! \file util.cc
    \ingroup OPTKING
    \brief miscellaneous
*/

#include "findif.h"

#include <physconst.h>

namespace psi { namespace findif {
namespace {
  int levi(int a, int b, int c){
    int val =0;
    int x =0, y=1,z=2;
    if(a==x && b==y && c==z) val=1;
    else if(a==y && b==z && c==x) val=1;
    else if(a==z && b==x && c==y) val=1;
    else if(a==x && b==z && c==y) val=-1;
    else if(a==y && b==x && c==z) val=-1;
    else if(a==z && b==y && c==x) val=-1;
    else val=0;

    return val;
    }
  }//anonymous namespace

bool ascending(const VIBRATION *vib1, const VIBRATION *vib2) {
  if (vib1->km < vib2->km)
    return true;
  else
    return false;
}

// function to print out (frequencies and normal modes) vector of vibrations
void print_vibrations(boost::shared_ptr<Molecule> mol, std::vector<VIBRATION *> modes) {

  char **irrep_lbls = mol->irrep_labels();
  int Natom = mol->natom();

  // compute harmonic frequencies, +/- in wavenumbers
  /* Convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  /* v = 1/(2 pi c) sqrt( eval ) */
  const double k_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg);
  const double cm_convert = 1.0/(2.0 * pc_pi * pc_c * 100.0);

  for (int i=0; i<modes.size(); ++i) {
    if(modes[i]->km < 0.0)
      modes[i]->cm = -1*cm_convert * sqrt(-k_convert * modes[i]->km);
    else
      modes[i]->cm =    cm_convert * sqrt( k_convert * modes[i]->km);
  }

  // Sort modes by increasing eigenvalues.
  sort(modes.begin(), modes.end(), ascending);

  // Print out frequencies and irreps to output file.
  outfile->Printf( "\n\t  Irrep      Harmonic Frequency   \n");
  outfile->Printf(   "\t                  (cm-1)          \n");
  outfile->Printf(   "\t-----------------------------------------------\n");

  for(int i=0; i<modes.size(); ++i) {
    if(modes[i]->cm < 0.0)
      outfile->Printf( "\t  %5s   %15.4fi \n", irrep_lbls[modes[i]->irrep], -modes[i]->cm);
    else
      outfile->Printf( "\t  %5s   %15.4f  \n", irrep_lbls[modes[i]->irrep], modes[i]->cm);
  }

  outfile->Printf(   "\t-----------------------------------------------\n");


  // Return list of frequencies to wavefunction object.
  boost::shared_ptr<Vector> freq_vector(new Vector(modes.size()));
  for (int i=0; i<modes.size(); ++i)
    freq_vector->set(i, modes[i]->cm);

  // Reture list of normal modes to wavefunction object.
  boost::shared_ptr<Vector> nm_vector(new Vector(3*Natom*modes.size()));
  int count = 0;
  for (int i=0; i<modes.size(); ++i) {
    freq_vector->set(i, modes[i]->cm);
    for (int a=0; a<Natom; ++a) {
        for (int xyz=0; xyz<3; xyz++) {
            nm_vector->set(count, modes[i]->lx[3*a+xyz]);
            count++;
        }
    }
  }

  Process::environment.set_frequencies(freq_vector);

  double sum = 0.0;
  for (int a=0; a<Natom; ++a)
     sum += mol->mass(a);

  // print out normal modes in format that WebMO likes
  outfile->Printf( "\n\tNormal Modes (non-mass-weighted).\n");
  outfile->Printf( "\tMolecular mass is %10.5f amu.\n", sum);
  outfile->Printf( "\tFrequencies in cm^-1; force constants in au.\n");

  for(int i=0; i<modes.size(); ++i) { // print descending order
    if (fabs(cm_convert * sqrt(k_convert * fabs(modes[i]->km))) < 5.0) continue;
    outfile->Printf("\n");
    if (modes[i]->km < 0.0)
      outfile->Printf( "   Frequency:      %8.2fi\n", cm_convert * sqrt(-k_convert * modes[i]->km));
    else
      outfile->Printf( "   Frequency:      %8.2f\n", cm_convert * sqrt(k_convert * modes[i]->km));

    outfile->Printf(   "   Force constant: %8.4f\n", modes[i]->km);

    //outfile->Printf(   "   IR Intensity: %8.2f\n", irint[i]*ir_prefactor);

    outfile->Printf( "\t     X       Y       Z           mass\t\n");
    for (int a=0; a<Natom; a++) {
      outfile->Printf( "  %s \t", mol->symbol(a).c_str() );

      for (int xyz=0; xyz<3; ++xyz)
        outfile->Printf( "%8.3f", modes[i]->lx[3*a+xyz]);

      outfile->Printf("%15.6f", mol->mass(a));

      outfile->Printf( "\n");
    }
  }

  // awkward, but need nirrep to free labels
  int Nirrep = mol->point_group()->char_table().nirrep();

  for (int i=0; i<Nirrep; ++i)
    free(irrep_lbls[i]);
  free(irrep_lbls);
}

// displaces from a reference geometry: geom += salclist[salc_i] * disp_i * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacement is DX/sqrt(mass)
void displace_cart(boost::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor));

  int nc = salclist[salc_i].ncomponent();

  for (int c=0; c<nc; ++c) {
    int a          = salclist[salc_i].component(c).atom;
    int xyz        = salclist[salc_i].component(c).xyz;
    double coef    = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// displaces from a reference geometry.
// geom += salclist[salc_i] * disp_i * disp_size + salclist[salc_j] * disp_j * disp_size
// disp_size is in mass-weighted coordinates; cartesian displacement is DX/sqrt(mass)
void displace_cart(boost::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int salc_j, int disp_factor_i, int disp_factor_j, double disp_size) {

  geom->set_name("Coord: " + to_string(salc_i) + ", Disp: " + to_string(disp_factor_i)
    + "Coord: " + to_string(salc_j) + ", Disp: " + to_string(disp_factor_j));

  int a, xyz;
  double coef;

  for (int c=0; c<salclist[salc_i].ncomponent(); ++c) {
    a    = salclist[salc_i].component(c).atom;
    xyz  = salclist[salc_i].component(c).xyz;
    coef = salclist[salc_i].component(c).coef;

    geom->add(0, a, xyz, disp_factor_i * disp_size * coef / sqrt(mol->mass(a)));
  }

  for (int c=0; c<salclist[salc_j].ncomponent(); ++c) {
    a    = salclist[salc_j].component(c).atom;
    xyz  = salclist[salc_j].component(c).xyz;
    coef = salclist[salc_j].component(c).coef;

    geom->add(0, a, xyz, disp_factor_j * disp_size * coef / sqrt(mol->mass(a)));
  }

  return;
}

// it's assumed columns are cartesian dimensions
void mass_weight_columns_plus_one_half(boost::shared_ptr<Molecule> mol, SharedMatrix B) {
  double u;

  for (int col=0; col<B->ncol(); ++col) {
    u = sqrt(mol->mass(col/3));
    for (int row=0; row<B->nrow(); ++row)
      B->set(row, col, B->get(row,col) * u);
  }
}

void displace_atom(SharedMatrix geom, const int atom, const int coord, const int sign, const double disp_size) {

  geom->add(0, atom, coord, sign * disp_size);

  return;
}

std::vector<SharedMatrix> atomic_displacements(boost::shared_ptr<Molecule> mol, Options &options)
{

  // This is the size in bohr because geometry is in bohr at this point
  // This equals 0.1 angstrom displacement
  double disp_size = options.get_double("DISP_SIZE");

  int natom = mol->natom();

  // Geometry seems to be in bohr at this point
  Matrix ref_geom_temp = mol->geometry();
  SharedMatrix ref_geom(ref_geom_temp.clone());

  std::vector< SharedMatrix > disp_geoms;

  // Generate displacements
  for(int atom=0; atom < natom; ++atom) {
    for(int coord=0; coord < 3; ++coord) {
      // plus displacement
      SharedMatrix p_geom(ref_geom->clone());
      displace_atom(p_geom, atom, coord, +1, disp_size);
      disp_geoms.push_back(p_geom);
      // minus displacement
      SharedMatrix m_geom(ref_geom->clone());
      displace_atom(m_geom, atom, coord, -1, disp_size);
      disp_geoms.push_back(m_geom);
    }
  }

  // put reference geometry in list
  // disp_geoms.push_back(ref_geom);

  return disp_geoms;

}

/*
 * mixed_atomic_displacements:
 *  Generates displaced geometries by displacing two of the 3n Cartesian coordinates
 *  at a  time.
 *  For every atom N1 the coordinate x1 is displaced +/- directions.
 *  For each of those displaced geometries all atoms N2 <= N1 and all
 *  coordinates {x2: x,y,z if N2<N1; x2 <x1 if N2=N1} are displaced in +/-
 *  directions.
 *
 *  ONLY the mixed displacement geometries are retuned, the atomic_displacement
 *  function that is used here should be used to return the single
 *  displacements.
 */

std::vector<SharedMatrix> mixed_atomic_displacements(
  boost::shared_ptr<Molecule> mol,
  Options &options)
{
  double disp_size = options.get_double("DISP_SIZE");
  int natom = mol->natom();

  std::vector<SharedMatrix> mixed_disp_geoms;
  std::vector<SharedMatrix> disp_geoms = atomic_displacements(mol,options);
  for(int x =0,disp_idx=0; x < natom*3; ++x)
  {
    SharedMatrix disp_x_p = disp_geoms[disp_idx++];
    for(int y =0; y <x; y++)
    {
      int atom = y/3;
      int cord = y%3;
      SharedMatrix disp_x_p_y_p(disp_x_p->clone());
      SharedMatrix disp_x_p_y_m(disp_x_p->clone());

      displace_atom(disp_x_p_y_p, atom,cord, +1,disp_size);
      mixed_disp_geoms.push_back(disp_x_p_y_p);
      displace_atom(disp_x_p_y_m, atom,cord, -1,disp_size);
      mixed_disp_geoms.push_back(disp_x_p_y_m);

    }

    SharedMatrix disp_x_m = disp_geoms[disp_idx++];
    for(int y =0; y <x; y++)
    {
      int atom = y/3;
      int cord = y%3;
      SharedMatrix disp_x_p_y_p(disp_x_p->clone());
      SharedMatrix disp_x_p_y_m(disp_x_p->clone());
      SharedMatrix disp_x_m_y_p(disp_x_m->clone());
      SharedMatrix disp_x_m_y_m(disp_x_m->clone());


      displace_atom(disp_x_m_y_p, atom,cord, +1,disp_size);
      mixed_disp_geoms.push_back(disp_x_m_y_p);
      displace_atom(disp_x_m_y_m, atom,cord, -1,disp_size);
      mixed_disp_geoms.push_back(disp_x_m_y_m);
    }

  }

  return mixed_disp_geoms;



}

/*
 *normal_mode_displacements
 *  Return a list of geometries displaced along normal modes by
 *  options.get_double("DISP_SIZE") in +/- directions.
 */

std::vector<SharedMatrix> normal_mode_displacement_vectors(
  boost::shared_ptr<Molecule> mol,
  Options &options,
  SharedMatrix F //Hessian matrix, read + set by python wrapper
  )
{
  int natom = mol->natom();
  SharedMatrix geom(new Matrix(natom,3));
  SharedMatrix geom_orig(new Matrix(natom,3));
  SharedMatrix M(new Matrix(natom*3,natom*3));
  M->zero();
  mol->move_to_com();
  geom_orig->copy(mol->geometry());
  geom_orig->print();
  geom->zero();
  //Mass-Weighting
  double total_mass=0.00;
  for(int i = 0; i < natom; i++){
    total_mass += mol->mass(i);
    for(int j = 0; j<3; j++){
      geom->set(i,j, geom_orig->get(i,j)*sqrt(mol->mass(i)));
      M->set((i*3+j),(i*3+j),1/sqrt(mol->mass(i)));
    }
  }
  outfile->Printf("TOTAL MASS == %lf (amu)\n", total_mass);
  // Diagonalize MOI tensor
  SharedMatrix I(new Matrix("Inertia Tensor",3,3));
  SharedMatrix Ievect(new Matrix("Inertia Tensor Eigenvectors",3,3));
  SharedVector Ieval(new Vector("Inertia Tensor Eigenvalues",3));
  I->copy(mol->inertia_tensor());
  I->diagonalize(Ievect,Ieval);

  // Build I inverse
  SharedMatrix Iinv(new Matrix("I^{-1}",3,3));
  SharedMatrix Itmp(new Matrix(3,3));
  Iinv->zero();
  for(int i =0; i < 3; i++){
    Iinv->set(i,i,(1.0/Ieval->get(i)));
  }
  Itmp->gemm(0,1,1.0,Iinv,Ievect,0.0);
  Iinv->gemm(0,0,1.0,Ievect,Itmp,0.0);

  // rotation-translation projector
  SharedMatrix P(new Matrix(natom*3,natom*3));
  //Generate rotation-translation vectors
  for(int i = 0; i < natom*3; i++){
    int icart = i%3;
    int iatom = i/3;
    int imass = mol->mass(iatom);

    P->set(i,i,1.0);

    for(int j = 0; j < natom*3; j++){
      int jcart = j%3;
      int jatom = j/3;
      int jmass = mol->mass(jatom);

      P->add(i,j,-1.0*sqrt(imass*jmass)/total_mass*(icart==jcart));

      for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
          for(int c = 0; c < 3; c++){
            for(int d = 0; d < 3; d++){
              P->add(i,j,
                -1.0*levi(a,b,icart)*geom->get(iatom,b)*Iinv->get(a,c)* levi(c,d,jcart)*geom->get(jatom,d));
            }
          }
        }
      }
    }
  }
  F->print("check-Hessian.test");
  P->print("check-Proj.test");
  SharedMatrix T(new Matrix(3*natom,3*natom));
  // Mass Weight F --> FM
  T->zero();
  T->gemm(0,0,1.0,M,F,0.0);
  F->gemm(0,0,1.0,T,M,0.0);

  //Project out rotations and translations
  T->zero();
  T->gemm(0,0,1.0,F,P,0.0);
  F->gemm(0,0,1.0,P,T,0.0);

  // Projected Hessian Evects (normal modes)
  SharedMatrix Fevec(new Matrix(natom*3,natom*3));
  // vibrational freq
  SharedVector freq(new Vector(natom*3));

  /* Convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  const double km_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg);
  /* v = 1/(2 pi c) sqrt( eval ) */
  const double cm_convert = 1.0/(2.0 * pc_pi * pc_c * 100.0);
  //diagonalize MW hessian
  F->diagonalize(Fevec,freq);
  // compute normal mode displacement vectors
  SharedMatrix Lx(new Matrix(natom*3,natom*3));
  T->zero();
  T->gemm(0,0,1.0,P,Fevec,0.0);
  Lx->gemm(0,0,1.0,M,T,0.0);
  SharedVector redmass( new Vector(natom*3) );
  double norm = 0.0;
  redmass->zero();
  for(int i = 0; i < 3*natom; i++){
    norm = 0.0;
    for(int j = 0; j < 3*natom; j++){
      norm += Lx->get(j,i)*Lx->get(j,i);
    }

    if(norm > 1e-3){
      redmass->set(i,1.0/norm);
      outfile->Printf("Reduced mass of mode %9.4f = %20.12f\n",cm_convert*(sqrt(km_convert*freq->get(i))),redmass->get(i));
    }
    norm=sqrt(norm);
    if(norm > 1e-3){
      for(int j = 0; j< 3*natom; j++){
        double normVal = Lx->get(j,i) * 1/norm;
        Lx->set(j,i,normVal);
      }
    }
  }
  //create displacement "vectors"
  std::vector<SharedMatrix> disp_vects;
  for(int i = 6; i < (3*natom); i++ ){
    SharedMatrix d_vec(new Matrix(natom,3));
    for(int j = 0; j<natom*3; j++){
      int jcart = j%3;
      int jatom = j/3;
      d_vec->set(jatom,jcart,Lx->get(j,i) );
    }
    disp_vects.push_back(d_vec);
  }

  return disp_vects;
}


/*
 *normal_mode_rms_amp_displacements
 *  Returns a vector of std::pair (displacement, size). Each normal mode
 *  has its own displacement amplitude which depends on temperature and the
 *  vibrational frequency.
 */
std::vector<std::pair<SharedMatrix,double>> normal_mode_rms_amp_displacements(
  boost::shared_ptr<Molecule> mol,
  Options &options,
  SharedMatrix F)
{
  int natom = mol->natom();
  SharedMatrix geom(new Matrix(natom,3));
  SharedMatrix geom_orig(new Matrix(natom,3));
  SharedMatrix M(new Matrix(natom*3,natom*3));
  M->zero();
  geom_orig->copy(mol->geometry());
  outfile->Printf("\nORIGNAL GEOMETRY\n");
  geom_orig->print();
  outfile->Printf("\nAfter COM GEOMETRY\n");
  mol->move_to_com();
  geom_orig->copy(mol->geometry());
  geom_orig->print();
  geom->zero();
  //Mass-Weighting
  double total_mass=0.00;
  for(int i = 0; i < natom; i++){
    total_mass += mol->mass(i);
    for(int j = 0; j<3; j++){
      geom->set(i,j, geom_orig->get(i,j)*sqrt(mol->mass(i)));
      M->set((i*3+j),(i*3+j),1/sqrt((mol->mass(i))/pc_au2amu));
    }
  }
  outfile->Printf("!!!! In Normal_mod_rms_amp_displacements\n");
  outfile->Printf("--->Total Mass = %lf\n",total_mass);

  // Diagonalize MOI tensor
  SharedMatrix I(new Matrix("Inertia Tensor",3,3));
  SharedMatrix Ievect(new Matrix("Inertia Tensor Eigenvectors",3,3));
  SharedVector Ieval(new Vector("Inertia Tensor Eigenvalues",3));
  I->copy(mol->inertia_tensor());
  I->diagonalize(Ievect,Ieval);

  // Build I inverse
  SharedMatrix Iinv(new Matrix("I^{-1}",3,3));
  SharedMatrix Itmp(new Matrix(3,3));
  Iinv->zero();
  for(int i =0; i < 3; i++){
    Iinv->set(i,i,(1.0/Ieval->get(i)));
  }
  Itmp->gemm(0,1,1.0,Iinv,Ievect,0.0);
  Iinv->gemm(0,0,1.0,Ievect,Itmp,0.0);

  // rotation-translation projector
  SharedMatrix P(new Matrix(natom*3,natom*3));
  //Generate rotation-translation vectors
  for(int i = 0; i < natom*3; i++){
    int icart = i%3;
    int iatom = i/3;
    int imass = mol->mass(iatom);

    P->set(i,i,1.0);

    for(int j = 0; j < natom*3; j++){
      int jcart = j%3;
      int jatom = j/3;
      int jmass = mol->mass(jatom);

      P->add(i,j,-1.0*sqrt(imass*jmass)/total_mass*(icart==jcart));

      for(int a = 0; a < 3; a++){
        for(int b = 0; b < 3; b++){
          for(int c = 0; c < 3; c++){
            for(int d = 0; d < 3; d++){
              P->add(i,j,
                -1.0*levi(a,b,icart)*
                geom->get(iatom,b)*Iinv->get(a,c)*
                levi(c,d,jcart)*geom->get(jatom,d)
                );
            }
          }
        }
      }
    }
  }
  F->print("check-Hessian.test");
  P->print("check-Proj.test");
  SharedMatrix T(new Matrix(3*natom,3*natom));
  // Mass Weight F --> FM
  T->zero();
  T->gemm(0,0,1.0,M,F,0.0);
  F->gemm(0,0,1.0,T,M,0.0);

  //Project out rotations and translations
  T->zero();
  T->gemm(0,0,1.0,F,P,0.0);
  F->gemm(0,0,1.0,P,T,0.0);

  // Diagonalize Fm --> Q(modes),v(freq)
  // Projected Hessian Evects (normal modes)
  SharedMatrix Fevec(new Matrix(natom*3,natom*3));
  // vibrational freq
  SharedVector freq(new Vector(natom*3));

  /* Convert evals from H/(kg bohr^2) to J/(kg m^2) = 1/s^2 */
  const double km_convert = pc_hartree2J/(pc_bohr2m * pc_bohr2m * pc_amu2kg * pc_au2amu);
  /* v = 1/(2 pi c) sqrt( eval ) */
  const double cm_convert = 1.0/(2.0 * pc_pi * pc_c * 100.0);
  //diagonalize MW hessian
  F->diagonalize(Fevec,freq);
  // compute normal mode displacement vectors
  SharedMatrix Lx(new Matrix(natom*3,natom*3));
  T->zero();
  T->gemm(0,0,1.0,P,Fevec,0.0);
  Lx->gemm(0,0,1.0,M,T,0.0);
  SharedVector redmass( new Vector(natom*3) );
  double norm = 0.0;
  redmass->zero();
  for(int i = 0; i < 3*natom; i++){
    norm = 0.0;
    for(int j = 0; j < 3*natom; j++){
      norm += Lx->get(j,i)*Lx->get(j,i);
    }

    if(norm > 1e-3){
      redmass->set(i,1.0/norm);
    }
    norm=sqrt(norm);
    if(norm > 1e-3){
      for(int j = 0; j< 3*natom; j++){
        double normVal = Lx->get(j,i) * 1/norm;
        Lx->set(j,i,normVal);
      }
    }
  }
  Lx->print("Lx.test");
  SharedVector freq_cm(new Vector(natom*3-6));
  for(int i = 6; i < (3*natom); i++ ){
    double omega = cm_convert*(sqrt(km_convert*freq->get(i)));
    freq_cm->set(i-6, omega); 
  }
  //used later while printing final output in zpvc-correction
  freq_cm->print("freq.dat");
  //create displacement "vectors" and compute displacement amplitude for each
  //wiberg's constants
  double c1 = 16.8576;
  double c2 = 0.719384;
  //double temp = options.get_double("T");
  double temp = 300.00;
  std::vector<std::pair<SharedMatrix,double>> disp_amp_pairs;
  for(int i = 6; i < (3*natom); i++ ){
    // if(freq->get(i) > 1e-5){
      double omega = cm_convert*(sqrt(km_convert*freq->get(i)));
      outfile->Printf("%f\n", omega);
      double amp = c2*omega/temp;
      amp = cosh(amp)/sinh(amp);
      amp = sqrt((c1/omega)*amp);
      amp = amp/pc_bohr2angstroms;
      SharedMatrix d_vec(new Matrix(natom,3));
      for(int j = 0; j<natom*3; j++){
        int jcart = j%3;
        int jatom = j/3;
        d_vec->set(jatom,jcart,Lx->get(j,i));
      }
      disp_amp_pairs.push_back(std::make_pair(d_vec,amp));
    // }
  }

  return disp_amp_pairs;
}
}}//namespace psi::findif
