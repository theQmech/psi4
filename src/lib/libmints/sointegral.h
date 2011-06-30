#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

#include "onebody.h"
#include "twobody.h"
#include "basisset.h"
#include "integral.h"
#include "sobasis.h"
#include "gshell.h"
#include "petitelist.h"
#include "wavefunction.h"
#include "cdsalclist.h"
#include "dcd.h"

#include <libqt/qt.h>
#include <vector>

#define DebugPrint 0

#if DebugPrint
#define dprintf(...) fprintf(outfile, __VA_ARGS__)
#else
#define dprintf(...)
#endif

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class Matrix;

template<typename T>
void swap_index(T& a, T& b) {
    T temp;
    temp = b;
    b = a;
    a = temp;
}

#define SWAP_INDEX(a, b) swap_index(a ## abs, b ## abs); swap_index(a ## rel, b ## rel); swap_index(a ## irrep, b ## irrep);

class OneBodySOInt
{
protected:
    boost::shared_ptr<OneBodyAOInt> ob_;
    const IntegralFactory* integral_;
    int deriv_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;

    void common_init();

public:
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>& , const boost::shared_ptr<IntegralFactory> &);
    OneBodySOInt(const boost::shared_ptr<OneBodyAOInt>& , const IntegralFactory*);
    virtual ~OneBodySOInt();

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;

    /**
     * Computes a one-electron integral matrix. Only works for symmetric operators
     * (multipole operators will not work).
     *
     * \param result Where the integrals are going.
     */
    virtual void compute(boost::shared_ptr<Matrix> result);

    /**
     * Computes one-electron integral matrices. Should be able to handle multipole operators
     *
     * \param results Where the integrals are going.
     */
    virtual void compute(std::vector<boost::shared_ptr<Matrix> > results);

    /**
     * Computes one-electron integral derivative matrices.
     *
     * \param result Where the integral derivatives are going.
     * \param cdsalcs The Cartesian displacement SALCs that you are interested in.
     */
    virtual void compute_deriv1(std::vector<boost::shared_ptr<Matrix> > result,
                                const CdSalcList& cdsalcs);
};

// Only include the following function if Doxygen is running to generate appropriate
// documentation.
#ifdef DOXYGEN
class TwoBodySOIntFunctor
{
public:
    void operator()(int pirrep, int pso, int qirrep, int qso, int rirrep, int rso, int sirrep, int sso, double value);
};
#endif

class TwoBodySOInt
{
protected:
    boost::shared_ptr<TwoBodyAOInt> tb_;
    boost::shared_ptr<IntegralFactory> integral_;

    boost::shared_ptr<SOBasisSet> b1_;
    boost::shared_ptr<SOBasisSet> b2_;
    boost::shared_ptr<SOBasisSet> b3_;
    boost::shared_ptr<SOBasisSet> b4_;

    size_t size_;
    double *buffer_;
    double **deriv_;

    int iirrepoff_[8], jirrepoff_[8], kirrepoff_[8], lirrepoff_[8];

    boost::shared_ptr<PetiteList> petite1_;
    boost::shared_ptr<PetiteList> petite2_;
    boost::shared_ptr<PetiteList> petite3_;
    boost::shared_ptr<PetiteList> petite4_;

    boost::shared_ptr<PointGroup> pg_;

    boost::shared_ptr<DCD> dcd_;

    bool only_totally_symmetric_;

    const CdSalcList* cdsalcs_;

    template<typename TwoBodySOIntFunctor>
    void provide_IJKL(int, int, int, int, TwoBodySOIntFunctor& body);

    template<typename TwoBodySOIntFunctor>
    void provide_IJKL_deriv1(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body);

    void common_init();
public:
    TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt>&,
                 const boost::shared_ptr<IntegralFactory>&);
    TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt>& aoint,
                 const boost::shared_ptr<IntegralFactory>& intfac,
                 const CdSalcList& cdsalcs);
    virtual ~TwoBodySOInt();

    bool only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(bool ots) { only_totally_symmetric_ = ots; }

    boost::shared_ptr<SOBasisSet> basis() const;
    boost::shared_ptr<SOBasisSet> basis1() const;
    boost::shared_ptr<SOBasisSet> basis2() const;
    boost::shared_ptr<SOBasisSet> basis3() const;
    boost::shared_ptr<SOBasisSet> basis4() const;

    const double *buffer() const { return buffer_; }

    // Normal integrals
    template<typename TwoBodySOIntFunctor>
    void compute_shell(const SOShellCombinationsIterator& shellIter, TwoBodySOIntFunctor& body) {
        compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s(),
                      body);
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell(int, int, int, int, TwoBodySOIntFunctor& body);

    // Derivative integrals
    template<typename TwoBodySOIntFunctor>
    void compute_shell_deriv1(const SOShellCombinationsIterator& shellIter, TwoBodySOIntFunctor& body) {
        compute_shell_deriv1(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s(),
                      body);
    }

    template<typename TwoBodySOIntFunctor>
    void compute_shell_deriv1(int, int, int, int, TwoBodySOIntFunctor& body);
};

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell(int uish, int ujsh, int uksh, int ulsh, TwoBodySOIntFunctor& body)
{
#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell overall");
#endif

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell setup");
#endif // MINTS_TIMER

    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->sotrans(uish);
    const SOTransform &t2 = b2_->sotrans(ujsh);
    const SOTransform &t3 = b3_->sotrans(uksh);
    const SOTransform &t4 = b4_->sotrans(ulsh);

    const int nso1 = b1_->nfunction(uish);
    const int nso2 = b2_->nfunction(ujsh);
    const int nso3 = b3_->nfunction(uksh);
    const int nso4 = b4_->nfunction(ulsh);
    const size_t nso = nso1*nso2*nso3*nso4;

    const int nao1 = b1_->naofunction(uish);
    const int nao2 = b2_->naofunction(ujsh);
    const int nao3 = b3_->naofunction(uksh);
    const int nao4 = b4_->naofunction(ulsh);
    const size_t nao = nao1*nao2*nao3*nao4;

    const int iatom = tb_->basis1()->shell(t1.aoshell[0].aoshell)->ncenter();
    const int jatom = tb_->basis2()->shell(t2.aoshell[0].aoshell)->ncenter();
    const int katom = tb_->basis3()->shell(t3.aoshell[0].aoshell)->ncenter();
    const int latom = tb_->basis4()->shell(t4.aoshell[0].aoshell)->ncenter();

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell zero buffer");
#endif // MINTS_TIMER

    memset(buffer_, 0, nso1*nso2*nso3*nso4*sizeof(double));

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell zero buffer");
#endif // MINTS_TIMER

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell setup");
#endif // MINTS_TIMER

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::compute_shell full shell transform");
#endif // MINTS_TIMER

    // Get the atomic stablizer (the first symmetry operation that maps the atom
    // onto itself.

    // These 3 sections are not shell specific so we can just use petite1_
    const unsigned short istablizer = petite1_->stablizer(iatom);
    const unsigned short jstablizer = petite1_->stablizer(jatom);
    const unsigned short kstablizer = petite1_->stablizer(katom);
    const unsigned short lstablizer = petite1_->stablizer(latom);

    const int istabdense = dcd_->bits_to_dense_numbering(istablizer);
    const int jstabdense = dcd_->bits_to_dense_numbering(jstablizer);
    const int kstabdense = dcd_->bits_to_dense_numbering(kstablizer);
    const int lstabdense = dcd_->bits_to_dense_numbering(lstablizer);

    const int ijstablizer = dcd_->intersection(istabdense, jstabdense);
    const int klstablizer = dcd_->intersection(kstabdense, lstabdense);
    const int ijklstablizer = dcd_->intersection(ijstablizer, klstablizer);

    const int* R_list = dcd_->dcr(istabdense, jstabdense);
    const int* S_list = dcd_->dcr(kstabdense, lstabdense);
    const int* T_list = dcd_->dcr(ijstablizer, klstablizer);

    const int R_size = R_list[0];
    const int S_size = S_list[0];
    const int T_size = T_list[0];

    // Check with Andy on this:
    int lambda_T = petite1_->nirrep() / dcd_->subgroup_dimensions(ijklstablizer);

    std::vector<int> sj_arr, sk_arr, sl_arr;

    //fprintf(outfile, "for (%d %d | %d %d) need to compute:\n", uish, ujsh, uksh, ulsh);
    int si = petite1_->unique_shell_map(uish, 0);
    const int siatom = tb_->basis1()->shell(si)->ncenter();

    for (int ij=1; ij <= R_size; ++ij) {
        int sj = petite2_->unique_shell_map(ujsh, R_list[ij]);
        const int sjatom = tb_->basis2()->shell(sj)->ncenter();

        for (int ijkl=1; ijkl <= T_size; ++ijkl) {
            int sk = petite3_->unique_shell_map(uksh, T_list[ijkl]);
            int llsh = petite4_->unique_shell_map(ulsh, T_list[ijkl]);
            const int skatom = tb_->basis3()->shell(sk)->ncenter();

            for (int kl=1; kl <= S_size; ++kl) {
                int sl = petite4_->shell_map(llsh, S_list[kl]);
                const int slatom = tb_->basis4()->shell(sl)->ncenter();

                // Check AM
                int total_am = tb_->basis1()->shell(si)->am() +
                        tb_->basis2()->shell(sj)->am() +
                        tb_->basis3()->shell(sk)->am() +
                        tb_->basis4()->shell(sl)->am();

//                fprintf(outfile, "\ttotal_am = %d atoms: %d %d %d %d\n", total_am, siatom, sjatom, skatom, slatom);

                if (!(total_am % 2) ||
                        (siatom != sjatom) ||
                        (sjatom != skatom) ||
                        (skatom != slatom)) {
                    sj_arr.push_back(sj);
                    sk_arr.push_back(sk);
                    sl_arr.push_back(sl);
                }
            }
        }
    }

//    fprintf(outfile, "\tlambda_T: %d\n", lambda_T);
//    fprintf(outfile, "\tgroup = %d, R_list %d S_list %d T_list %d\n", group, R_list, S_list, T_list);
//    fprintf(outfile, "\tistablizer: %d\n", istabdense);
//    fprintf(outfile, "\tjstablizer: %d\n", jstabdense);
//    fprintf(outfile, "\tkstablizer: %d\n", kstabdense);
//    fprintf(outfile, "\tlstablizer: %d\n", lstabdense);
//    fprintf(outfile, "\tijstablizer: %d\n", ijstablizer);
//    fprintf(outfile, "\tklstablizer: %d\n", klstablizer);
//    fprintf(outfile, "\tR.size = %d\n", R_size);
//    for (int i=1; i<=R_size; ++i)
//        fprintf(outfile, "\t%d\n", R_list[i]);
//    fprintf(outfile, "\tS.size = %d\n", S_size);
//    for (int i=1; i<=S_size; ++i)
//        fprintf(outfile, "\t%d\n", S_list[i]);
//    fprintf(outfile, "\tT.size = %d\n", T_size);
//    for (int i=1; i<=T_size; ++i)
//        fprintf(outfile, "\t%d\n", T_list[i]);
//    fprintf(outfile, "\tR_list: ");
//    petite1_->print_group(R_list);
//    fprintf(outfile, "\tS_list: ");
//    petite1_->print_group(S_list);
//    fprintf(outfile, "\tT_list: ");
//    petite1_->print_group(T_list);
//    for (int i=0; i<sj_arr.size(); ++i) {
//        fprintf(outfile, "\t(%d %d | %d %d)\n", si, sj_arr[i], sk_arr[i], sl_arr[i]);
//    }
//    fflush(outfile);

    // Compute integral using si, sj_arr, sk_arr, sl_arr
    // Loop over unique quartets
    const AOTransform& s1 = b1_->aotrans(si);

    for (int n=0; n<sj_arr.size(); ++n) {
        int sj = sj_arr[n];
        int sk = sk_arr[n];
        int sl = sl_arr[n];

        const AOTransform& s2 = b2_->aotrans(sj);
        const AOTransform& s3 = b3_->aotrans(sk);
        const AOTransform& s4 = b4_->aotrans(sl);

        int ns1so = s1.soshell.size();
        int ns2so = s2.soshell.size();
        int ns3so = s3.soshell.size();
        int ns4so = s4.soshell.size();

        // Compute this unique AO shell
        tb_->compute_shell(si, sj, sk, sl);

        for (int itr=0; itr<ns1so; itr++) {
            const AOTransformFunction &ifunc = s1.soshell[itr];
            double icoef = ifunc.coef;
            int iaofunc = ifunc.aofunc;
            int isofunc = b1_->function_offset_within_shell(uish,
                                                            ifunc.irrep)
                    + ifunc.sofunc;
            int iaooff = iaofunc;
            int isooff = isofunc;

            for (int jtr=0; jtr<ns2so; jtr++) {
                const AOTransformFunction &jfunc = s2.soshell[jtr];
                double jcoef = jfunc.coef * icoef;
                int jaofunc = jfunc.aofunc;
                int jsofunc = b2_->function_offset_within_shell(ujsh,
                                                                jfunc.irrep)
                        + jfunc.sofunc;
                int jaooff = iaooff*nao2 + jaofunc;
                int jsooff = isooff*nso2 + jsofunc;

                for (int ktr=0; ktr<ns3so; ktr++) {
                    const AOTransformFunction &kfunc = s3.soshell[ktr];
                    double kcoef = kfunc.coef * jcoef;
                    int kaofunc = kfunc.aofunc;
                    int ksofunc = b3_->function_offset_within_shell(uksh,
                                                                    kfunc.irrep)
                            + kfunc.sofunc;
                    int kaooff = jaooff*nao3 + kaofunc;
                    int ksooff = jsooff*nso3 + ksofunc;

                    for (int ltr=0; ltr<ns4so; ltr++) {
                        const AOTransformFunction &lfunc = s4.soshell[ltr];
                        double lcoef = lfunc.coef * kcoef;
                        int laofunc = lfunc.aofunc;
                        int lsofunc = b4_->function_offset_within_shell(ulsh,
                                                                        lfunc.irrep)
                                + lfunc.sofunc;
                        int laooff = kaooff*nao4 + laofunc;
                        int lsooff = ksooff*nso4 + lsofunc;
                        // If you're doing the two-stage SO integral uncomment the next line
//                        fprintf(outfile, "\t\tbuffer_[%d] += %d * %f * %f symms: %d %d %d %d\n", lsooff, lambda_T, lcoef, aobuff[laooff],
//                                ifunc.irrep, jfunc.irrep, kfunc.irrep, lfunc.irrep);

                        if ((ifunc.irrep ^ jfunc.irrep) == (kfunc.irrep ^ lfunc.irrep)) {
//                            fprintf(outfile, "\t\tadded\n");
                            buffer_[lsooff] += lambda_T * lcoef * aobuff[laooff];
                        }
                    }
                }
            }
        }
    }

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell full shell transform");
#endif // MINTS_TIMER

    provide_IJKL(uish, ujsh, uksh, ulsh, body);

#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::compute_shell overall");
#endif
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL overall");
#endif

    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->sotrans(ish);
    const SOTransform &t2 = b2_->sotrans(jsh);
    const SOTransform &t3 = b3_->sotrans(ksh);
    const SOTransform &t4 = b4_->sotrans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int n1 = b1_->nfunction(ish);
    int n2 = b2_->nfunction(jsh);
    int n3 = b3_->nfunction(ksh);
    int n4 = b4_->nfunction(lsh);

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    int itr, itrfunc;
    int jtr, jtrfunc;
    int ktr, ktrfunc;
    int ltr, ltrfunc;

    for (itr=0, itrfunc=0; itr<n1; itr++, itrfunc++) {

        int ifunc = b1_->function(ish) + itr;
        int isym = b1_->irrep(ifunc);
        int irel = b1_->function_within_irrep(ifunc);
        int iabs = iirrepoff_[isym] + irel;
        int isooff = itr;

        for (jtr=0; jtr<n2; jtr++) {

            int jfunc = b2_->function(jsh) + jtr;
            int jsym = b2_->irrep(jfunc);
            int jrel = b2_->function_within_irrep(jfunc);
            int jabs = jirrepoff_[jsym] + jrel;
            int jsooff = isooff*nso2 + jtr;

            for (ktr=0; ktr<n3; ktr++) {

                int kfunc = b3_->function(ksh) + ktr;
                int ksym = b3_->irrep(kfunc);
                int krel = b3_->function_within_irrep(kfunc);
                int kabs = kirrepoff_[ksym] + krel;
                int ksooff = jsooff*nso3 + ktr;

                for (ltr=0; ltr<n4; ltr++) {

                    int lfunc = b4_->function(lsh) + ltr;
                    int lsym = b4_->irrep(lfunc);
                    int lrel = b4_->function_within_irrep(lfunc);
                    int labs = lirrepoff_[lsym] + lrel;
                    int lsooff = ksooff*nso4 + ltr;

                    int iiabs = iabs;
                    int jjabs = jabs;
                    int kkabs = kabs;
                    int llabs = labs;

                    int iiirrep = isym;
                    int jjirrep = jsym;
                    int kkirrep = ksym;
                    int llirrep = lsym;

                    int iirel = irel;
                    int jjrel = jrel;
                    int kkrel = krel;
                    int llrel = lrel;

                    if (fabs(buffer_[lsooff]) > 1.0e-14) {
                        if (ish == jsh) {
                            if (iabs < jabs)
                            continue;

                            if (ksh == lsh) {
                                if (kabs < labs)
                                continue;
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    if (ish == ksh)   // IIII case
                                    continue;
                                    else {            // IIJJ case
                                        SWAP_INDEX(ii, kk);
                                        SWAP_INDEX(jj, ll);
                                    }
                                }
                            }
                            else{                     // IIJK case
                                if (labs > kabs) {
                                    SWAP_INDEX(kk, ll);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                        }
                        else {
                            if (ksh == lsh) {         // IJKK case
                                if (kabs < labs)
                                continue;
                                if (iabs < jabs) {
                                    SWAP_INDEX(ii, jj);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                            else {                   // IJIJ case
                                if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs))
                                continue;
                                // IJKL case
                                if (iabs < jabs) {
                                    SWAP_INDEX(ii, jj);
                                }
                                if (kabs < labs) {
                                    SWAP_INDEX(kk, ll);
                                }
                                if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                        }

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL functor");
#endif
                        // func off/on
                        body(iiabs, jjabs, kkabs, llabs,
                             iiirrep, iirel,
                             jjirrep, jjrel,
                             kkirrep, kkrel,
                             llirrep, llrel,
                             buffer_[lsooff]);
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL functor");
#endif
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL overall");
#endif
}

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::compute_shell_deriv1(int uish, int ujsh, int uksh, int ulsh, TwoBodySOIntFunctor& body)
{
    dprintf("usi %d usj %d usk %d usk %d\n", uish, ujsh, uksh, ulsh);
    const double *aobuffer = tb_->buffer();

    const SOTransform &t1 = b1_->sotrans(uish);
    const SOTransform &t2 = b2_->sotrans(ujsh);
    const SOTransform &t3 = b3_->sotrans(uksh);
    const SOTransform &t4 = b4_->sotrans(ulsh);

    const int nso1 = b1_->nfunction(uish);
    const int nso2 = b2_->nfunction(ujsh);
    const int nso3 = b3_->nfunction(uksh);
    const int nso4 = b4_->nfunction(ulsh);
    const size_t nso = nso1*nso2*nso3*nso4;

    const int nao1 = b1_->naofunction(uish);
    const int nao2 = b2_->naofunction(ujsh);
    const int nao3 = b3_->naofunction(uksh);
    const int nao4 = b4_->naofunction(ulsh);
    const size_t nao = nao1*nao2*nao3*nao4;

    const int iatom = tb_->basis1()->shell(t1.aoshell[0].aoshell)->ncenter();
    const int jatom = tb_->basis2()->shell(t2.aoshell[0].aoshell)->ncenter();
    const int katom = tb_->basis3()->shell(t3.aoshell[0].aoshell)->ncenter();
    const int latom = tb_->basis4()->shell(t4.aoshell[0].aoshell)->ncenter();

    // These 3 sections are not shell specific so we can just use petite1_
    const unsigned short istablizer = petite1_->stablizer(iatom);
    const unsigned short jstablizer = petite1_->stablizer(jatom);
    const unsigned short kstablizer = petite1_->stablizer(katom);
    const unsigned short lstablizer = petite1_->stablizer(latom);

    const int istabdense = dcd_->bits_to_dense_numbering(istablizer);
    const int jstabdense = dcd_->bits_to_dense_numbering(jstablizer);
    const int kstabdense = dcd_->bits_to_dense_numbering(kstablizer);
    const int lstabdense = dcd_->bits_to_dense_numbering(lstablizer);

    const int ijstablizer = dcd_->intersection(istabdense, jstabdense);
    const int klstablizer = dcd_->intersection(kstabdense, lstabdense);
    const int ijklstablizer = dcd_->intersection(ijstablizer, klstablizer);

    const int* R_list = dcd_->dcr(istabdense, jstabdense);
    const int* S_list = dcd_->dcr(kstabdense, lstabdense);
    const int* T_list = dcd_->dcr(ijstablizer, klstablizer);

    const int R_size = R_list[0];
    const int S_size = S_list[0];
    const int T_size = T_list[0];
    int lambda_T = petite1_->nirrep() / dcd_->subgroup_dimensions(ijklstablizer);

    std::vector<int> sj_arr, sk_arr, sl_arr;

    int si = petite1_->unique_shell_map(uish, 0);
    const int siatom = tb_->basis1()->shell(si)->ncenter();

    for (int ij=1; ij <= R_size; ++ij) {
        int sj = petite2_->unique_shell_map(ujsh, R_list[ij]);
        const int sjatom = tb_->basis2()->shell(sj)->ncenter();

        for (int ijkl=1; ijkl <= T_size; ++ijkl) {
            int sk = petite3_->unique_shell_map(uksh, T_list[ijkl]);
            int llsh = petite4_->unique_shell_map(ulsh, T_list[ijkl]);
            const int skatom = tb_->basis3()->shell(sk)->ncenter();

            for (int kl=1; kl <= S_size; ++kl) {
                int sl = petite4_->shell_map(llsh, S_list[kl]);
                const int slatom = tb_->basis4()->shell(sl)->ncenter();

                // Check AM
                int total_am = tb_->basis1()->shell(si)->am() +
                        tb_->basis2()->shell(sj)->am() +
                        tb_->basis3()->shell(sk)->am() +
                        tb_->basis4()->shell(sl)->am();

                if (!(total_am % 2) ||
                        (siatom != sjatom) ||
                        (sjatom != skatom) ||
                        (skatom != slatom)) {
                    sj_arr.push_back(sj);
                    sk_arr.push_back(sk);
                    sl_arr.push_back(sl);
                }
            }
        }
    }

    // Obtain SALC transformation objects
    // This probably won't work. I'll probably need
    // the SALC for the petite list of shells.
    // If so, this will move into the for loop below.

    // Compute integral using si, sj_arr, sk_arr, sl_arr
    // Loop over unique quartets
    const AOTransform& s1 = b1_->aotrans(si);
    const CdSalcWRTAtom& c1 = cdsalcs_->atom_salc(siatom);

    // Zero out SALC memory
    for (int i=0; i<cdsalcs_->ncd(); ++i)
        memset(deriv_[i], 0, sizeof(double)*nso);

    for (int n=0; n<sj_arr.size(); ++n) {
        int sj = sj_arr[n];
        int sk = sk_arr[n];
        int sl = sl_arr[n];

        dprintf("si %d sj %d sk %d sl %d\n", si, sj, sk, sl);

        const AOTransform& s2 = b2_->aotrans(sj);
        const AOTransform& s3 = b3_->aotrans(sk);
        const AOTransform& s4 = b4_->aotrans(sl);

        const CdSalcWRTAtom& c2 = cdsalcs_->atom_salc(tb_->basis2()->shell(sj)->ncenter());
        const CdSalcWRTAtom& c3 = cdsalcs_->atom_salc(tb_->basis3()->shell(sk)->ncenter());
        const CdSalcWRTAtom& c4 = cdsalcs_->atom_salc(tb_->basis4()->shell(sl)->ncenter());

        int ns1so = s1.soshell.size();
        int ns2so = s2.soshell.size();
        int ns3so = s3.soshell.size();
        int ns4so = s4.soshell.size();

        // Compute this unique AO shell
        tb_->compute_shell_deriv1(si, sj, sk, sl);

        // First loop over SO transformation
        for (int itr=0; itr<ns1so; itr++) {
            const AOTransformFunction &ifunc = s1.soshell[itr];
            double icoef = ifunc.coef;
            int iaofunc = ifunc.aofunc;
            int isofunc = b1_->function_offset_within_shell(uish,
                                                            ifunc.irrep)
                    + ifunc.sofunc;
            int iaooff = iaofunc;
            int isooff = isofunc;

            for (int jtr=0; jtr<ns2so; jtr++) {
                const AOTransformFunction &jfunc = s2.soshell[jtr];
                double jcoef = jfunc.coef * icoef;
                int jaofunc = jfunc.aofunc;
                int jsofunc = b2_->function_offset_within_shell(ujsh,
                                                                jfunc.irrep)
                        + jfunc.sofunc;
                int jaooff = iaooff*nao2 + jaofunc;
                int jsooff = isooff*nso2 + jsofunc;

                for (int ktr=0; ktr<ns3so; ktr++) {
                    const AOTransformFunction &kfunc = s3.soshell[ktr];
                    double kcoef = kfunc.coef * jcoef;
                    int kaofunc = kfunc.aofunc;
                    int ksofunc = b3_->function_offset_within_shell(uksh,
                                                                    kfunc.irrep)
                            + kfunc.sofunc;
                    int kaooff = jaooff*nao3 + kaofunc;
                    int ksooff = jsooff*nso3 + ksofunc;

                    for (int ltr=0; ltr<ns4so; ltr++) {
                        const AOTransformFunction &lfunc = s4.soshell[ltr];
                        double lcoef = lfunc.coef * kcoef;
                        int laofunc = lfunc.aofunc;
                        int lsofunc = b4_->function_offset_within_shell(ulsh,
                                                                        lfunc.irrep)
                                + lfunc.sofunc;
                        int laooff = kaooff*nao4 + laofunc;
                        int lsooff = ksooff*nso4 + lsofunc;

                        int total_symmetry = ifunc.irrep ^ jfunc.irrep ^ kfunc.irrep ^ lfunc.irrep;

                        // If we're only interested in totally symmetric derivatives, skip all others.
                        if (only_totally_symmetric_ == true && total_symmetry != 0)
                            continue;

                        // OK, the integral at 12 aobuff[laooff] needs
                        // to have lambda_T * lcoef * salcCoef applied to
                        // it and saved to salc index
                        // Special case is center B which must be determined
                        // from centers A, C, and D.
                        double A[3], B[3], C[3], D[3];

                        A[0] = aobuffer[0*nao+laooff]; A[1] = aobuffer[1*nao+laooff]; A[2] = aobuffer[2*nao+laooff];
                        C[0] = aobuffer[3*nao+laooff]; C[1] = aobuffer[4*nao+laooff]; C[2] = aobuffer[5*nao+laooff];
                        D[0] = aobuffer[6*nao+laooff]; D[1] = aobuffer[7*nao+laooff]; D[2] = aobuffer[8*nao+laooff];

                        // Use translational invariance to determine B
                        B[0] = -(A[0] + C[0] + D[0]);  B[1] = -(A[1] + C[1] + D[1]);  B[2] = -(A[2] + C[2] + D[2]);

                        A[0] *= lambda_T * lcoef;
                        A[1] *= lambda_T * lcoef;
                        A[2] *= lambda_T * lcoef;
                        B[0] *= lambda_T * lcoef;
                        B[1] *= lambda_T * lcoef;
                        B[2] *= lambda_T * lcoef;
                        C[0] *= lambda_T * lcoef;
                        C[1] *= lambda_T * lcoef;
                        C[2] *= lambda_T * lcoef;
                        D[0] *= lambda_T * lcoef;
                        D[1] *= lambda_T * lcoef;
                        D[2] *= lambda_T * lcoef;

                        dprintf("so' derivatives: A[0] %+lf A[1] %+lf A[2] %+lf\n"
                                         "                 B[0] %+lf B[1] %+lf B[2] %+lf\n"
                                         "                 C[0] %+lf C[1] %+lf C[2] %+lf\n"
                                         "                 D[0] %+lf D[1] %+lf D[2] %+lf\n"
                                         "lsooff: %d\n"
                                         "iirrep: %d jirrep: %d kirrep: %d lirrep: %d combined: %d\n",
                                A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2], D[0], D[1], D[2], lsooff,
                                ifunc.irrep, jfunc.irrep, kfunc.irrep, lfunc.irrep,
                                ifunc.irrep ^ jfunc.irrep ^ kfunc.irrep ^ lfunc.irrep);

                        // For each center apply the so transform and salc coef.
                        // Ax
                        for (int nx=0; nx<c1.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = c1.x(nx);
                            double temp = element.coef * A[0];
                            dprintf("Ax SALC#%d pfac %lf, A[0] %lf, contr %lf\n", element.salc, element.coef, A[0], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Ay
                        for (int ny=0; ny<c1.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = c1.y(ny);
                            double temp = element.coef * A[1];
                            dprintf("Ay SALC#%d pfac %lf, A[1] %lf, contr %lf\n", element.salc, element.coef, A[1], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Az
                        for (int nz=0; nz<c1.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = c1.z(nz);
                            double temp = element.coef * A[2];
                            dprintf("Az SALC#%d pfac %lf, A[2] %lf, contr %lf\n", element.salc, element.coef, A[2], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Bx
                        for (int nx=0; nx<c2.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = c2.x(nx);
                            double temp = element.coef * B[0];
                            dprintf("Bx SALC#%d pfac %lf, B[0] %lf, contr %lf\n", element.salc, element.coef, B[0], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // By
                        for (int ny=0; ny<c2.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = c2.y(ny);
                            double temp = element.coef * B[1];
                            dprintf("By SALC#%d pfac %lf, B[1] %lf, contr %lf\n", element.salc, element.coef, B[1], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Bz
                        for (int nz=0; nz<c2.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = c2.z(nz);
                            double temp = element.coef * B[2];
                            dprintf("Bz SALC#%d pfac %lf, B[2] %lf, contr %lf\n", element.salc, element.coef, B[2], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Cx
                        for (int nx=0; nx<c3.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = c3.x(nx);
                            double temp = element.coef * C[0];
                            dprintf("Cx SALC#%d pfac %lf, C[0] %lf, contr %lf\n", element.salc, element.coef, C[0], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Cy
                        for (int ny=0; ny<c3.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = c3.y(ny);
                            double temp = element.coef * C[1];
                            dprintf("Cy SALC#%d pfac %lf, C[1] %lf, contr %lf\n", element.salc, element.coef, C[1], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Cz
                        for (int nz=0; nz<c3.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = c3.z(nz);
                            double temp = element.coef * C[2];
                            dprintf("Cz SALC#%d pfac %lf, C[2] %lf, contr %lf\n", element.salc, element.coef, C[2], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Dx
                        for (int nx=0; nx<c4.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = c4.x(nx);
                            double temp = element.coef * D[0];
                            dprintf("Dx SALC#%d pfac %lf, D[0] %lf, contr %lf\n", element.salc, element.coef, D[0], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Dy
                        for (int ny=0; ny<c4.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = c4.y(ny);
                            double temp = element.coef * D[1];
                            dprintf("Dy SALC#%d pfac %lf, D[1] %lf, contr %lf\n", element.salc, element.coef, D[1], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        // Dz
                        for (int nz=0; nz<c4.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = c4.z(nz);
                            double temp = element.coef * D[2];
                            dprintf("Dz SALC#%d pfac %lf, D[2] %lf, contr %lf\n", element.salc, element.coef, D[2], temp);
                            if (total_symmetry == element.irrep && fabs(temp) > 1.0e-10)
                                deriv_[element.salc][lsooff] += temp;
                            dprintf(" val: %lf\n", deriv_[element.salc][lsooff]);
                        }

                        dprintf("\n");
                    }
                }
            }
        }
    }

    provide_IJKL_deriv1(uish, ujsh, uksh, ulsh, body);
} // function

template<typename TwoBodySOIntFunctor>
void TwoBodySOInt::provide_IJKL_deriv1(int ish, int jsh, int ksh, int lsh, TwoBodySOIntFunctor& body)
{
#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL overall");
#endif

    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->sotrans(ish);
    const SOTransform &t2 = b2_->sotrans(jsh);
    const SOTransform &t3 = b3_->sotrans(ksh);
    const SOTransform &t4 = b4_->sotrans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int n1 = b1_->nfunction(ish);
    int n2 = b2_->nfunction(jsh);
    int n3 = b3_->nfunction(ksh);
    int n4 = b4_->nfunction(lsh);

    const SOTransformShell &s1 = t1.aoshell[0];
    const SOTransformShell &s2 = t2.aoshell[0];
    const SOTransformShell &s3 = t3.aoshell[0];
    const SOTransformShell &s4 = t4.aoshell[0];

    int itr, itrfunc;
    int jtr, jtrfunc;
    int ktr, ktrfunc;
    int ltr, ltrfunc;

    for (itr=0, itrfunc=0; itr<n1; itr++, itrfunc++) {

        int ifunc = b1_->function(ish) + itr;
        int isym = b1_->irrep(ifunc);
        int irel = b1_->function_within_irrep(ifunc);
        int iabs = iirrepoff_[isym] + irel;
        int isooff = itr;

        for (jtr=0; jtr<n2; jtr++) {

            int jfunc = b2_->function(jsh) + jtr;
            int jsym = b2_->irrep(jfunc);
            int jrel = b2_->function_within_irrep(jfunc);
            int jabs = jirrepoff_[jsym] + jrel;
            int jsooff = isooff*nso2 + jtr;

            for (ktr=0; ktr<n3; ktr++) {

                int kfunc = b3_->function(ksh) + ktr;
                int ksym = b3_->irrep(kfunc);
                int krel = b3_->function_within_irrep(kfunc);
                int kabs = kirrepoff_[ksym] + krel;
                int ksooff = jsooff*nso3 + ktr;

                for (ltr=0; ltr<n4; ltr++) {

                    int lfunc = b4_->function(lsh) + ltr;
                    int lsym = b4_->irrep(lfunc);
                    int lrel = b4_->function_within_irrep(lfunc);
                    int labs = lirrepoff_[lsym] + lrel;
                    int lsooff = ksooff*nso4 + ltr;

                    int iiabs = iabs;
                    int jjabs = jabs;
                    int kkabs = kabs;
                    int llabs = labs;

                    int iiirrep = isym;
                    int jjirrep = jsym;
                    int kkirrep = ksym;
                    int llirrep = lsym;

                    int iirel = irel;
                    int jjrel = jrel;
                    int kkrel = krel;
                    int llrel = lrel;

                    if (ish == jsh) {
                        if (iabs < jabs)
                            continue;

                        if (ksh == lsh) {
                            if (kabs < labs)
                                continue;
                            if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                if (ish == ksh)   // IIII case
                                    continue;
                                else {            // IIJJ case
                                    SWAP_INDEX(ii, kk);
                                    SWAP_INDEX(jj, ll);
                                }
                            }
                        }
                        else{                     // IIJK case
                            if (labs > kabs) {
                                SWAP_INDEX(kk, ll);
                            }
                            if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                SWAP_INDEX(ii, kk);
                                SWAP_INDEX(jj, ll);
                            }
                        }
                    }
                    else {
                        if (ksh == lsh) {         // IJKK case
                            if (kabs < labs)
                                continue;
                            if (iabs < jabs) {
                                SWAP_INDEX(ii, jj);
                            }
                            if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                SWAP_INDEX(ii, kk);
                                SWAP_INDEX(jj, ll);
                            }
                        }
                        else {                   // IJIJ case
                            if (ish == ksh && jsh == lsh && INDEX2(iabs, jabs) < INDEX2(kabs, labs))
                                continue;
                            // IJKL case
                            if (iabs < jabs) {
                                SWAP_INDEX(ii, jj);
                            }
                            if (kabs < labs) {
                                SWAP_INDEX(kk, ll);
                            }
                            if (INDEX2(iabs, jabs) < INDEX2(kabs, labs)) {
                                SWAP_INDEX(ii, kk);
                                SWAP_INDEX(jj, ll);
                            }
                        }
                    }

#ifdef MINTS_TIMER
    timer_on("TwoBodySOInt::provide_IJKL functor");
#endif
                    for (int i=0; i<cdsalcs_->ncd(); ++i) {
                        if (fabs(deriv_[i][lsooff]) > 1.0e-14)
                            body(i, iiabs, jjabs, kkabs, llabs,
                                 iiirrep, iirel,
                                 jjirrep, jjrel,
                                 kkirrep, kkrel,
                                 llirrep, llrel,
                                 deriv_[i][lsooff]);
                    }
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL functor");
#endif
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("TwoBodySOInt::provide_IJKL overall");
#endif
}

typedef boost::shared_ptr<OneBodySOInt> SharedOneBodySOInt;
typedef boost::shared_ptr<TwoBodySOInt> SharedTwoBodySOInt;

}

#endif // _psi_src_lib_libmints_sointegral_h_
