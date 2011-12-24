#include "integraltransform.h"
#include "mospace.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libmints/matrix.h>
#include <libmints/molecule.h>
#include <libmints/wavefunction.h>
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace boost;
using namespace psi;

IntegralTransform::IntegralTransform(shared_ptr<Wavefunction> wfn,
                                     SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
            initialized_(false),
            psio_(_default_psio_lib_),
            wfn_(wfn),
            transformationType_(transformationType),
            uniqueSpaces_(spaces),
            moOrdering_(moOrdering),
            outputType_(outputType),
            frozenOrbitals_(frozenOrbitals),
            alreadyPresorted_(false),
            dpdIntFile_(PSIF_LIBTRANS_DPD),
            aHtIntFile_(PSIF_LIBTRANS_A_HT),
            bHtIntFile_(PSIF_LIBTRANS_B_HT),
            nTriSo_(0),
            nTriMo_(0),
            nfzc_(0),
            nfzv_(0),
            spaces_(0),
            labels_(0),
            tolerance_(1.0E-16),
            moIntFileAA_(0),
            moIntFileAB_(0),
            moIntFileBB_(0),
            myDPDNum_(1),
            print_(1),
            zeros_(0),
            sosym_(0),
            cacheFiles_(0),
            cacheList_(0),
            Ca_(wfn->Ca()),
            Cb_(wfn->Cb()),
            keepIwlSoInts_(false),
            keepIwlMoTpdm_(true),
            keepDpdSoInts_(false),
            keepDpdMoTpdm_(true),
            keepHtInts_(true),
            keepHtTpdm_(true),
            tpdmAlreadyPresorted_(false)
{
    // Implement set/get functions to customize any of this stuff.  Delayed initialization
    // is possible in case any of these variables need to be changed before setup.
    printTei_      = print_ > 5;
    useIWL_        = outputType_ == IWLAndDPD || outputType_ == IWLOnly;
    useDPD_        = outputType_ == IWLAndDPD || outputType_ == DPDOnly;
    iwlAAIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    iwlABIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    iwlBBIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    labels_  = Process::environment.molecule()->irrep_labels();
    nirreps_ = wfn->nirrep();
    nmo_     = wfn->nmo();
    nso_     = wfn->nso();
    sopi_    = wfn->nsopi();
    mopi_    = wfn->nmopi();
    clsdpi_  = wfn->doccpi();
    openpi_  = wfn->soccpi();
    frzcpi_  = wfn->frzcpi();
    frzvpi_  = wfn->frzvpi();

    common_moinfo_initialize();

    if(init) initialize();
}


IntegralTransform::IntegralTransform(SharedMatrix c,
                                     SharedMatrix i,
                                     SharedMatrix a,
                                     SharedMatrix v,
                                     SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
    initialized_(false),
    psio_(_default_psio_lib_),
    transformationType_(transformationType),
    uniqueSpaces_(spaces),
    moOrdering_(moOrdering),
    outputType_(outputType),
    frozenOrbitals_(frozenOrbitals),
    alreadyPresorted_(false),
    dpdIntFile_(PSIF_LIBTRANS_DPD),
    aHtIntFile_(PSIF_LIBTRANS_A_HT),
    bHtIntFile_(PSIF_LIBTRANS_B_HT),
    nirreps_(0),
    nmo_(0),
    nso_(0),
    nTriSo_(0),
    nTriMo_(0),
    nfzc_(0),
    nfzv_(0),
    spaces_(0),
    labels_(0),
    tolerance_(1.0E-16),
    memory_(250 * 1024 * 1024),
    moIntFileAA_(0),
    moIntFileAB_(0),
    moIntFileBB_(0),
    myDPDNum_(1),
    print_(1),
    zeros_(0),
    sopi_(0),
    sosym_(0),
    mopi_(0),
    clsdpi_(0),
    openpi_(0),
    frzcpi_(0),
    frzvpi_(0),
    cacheFiles_(0),
    cacheList_(0),
    aCorrToPitzer_(0),
    bCorrToPitzer_(0),
    keepIwlSoInts_(false),
    keepIwlMoTpdm_(true),
    keepDpdSoInts_(false),
    keepDpdMoTpdm_(true),
    keepHtInts_(true),
    keepHtTpdm_(true),
    tpdmAlreadyPresorted_(false)
{
//    wfn_           = NULL;
//    _Ca            = NULL;
//    _Cb            = NULL;
    printTei_      = print_ > 5;
    useIWL_        = outputType_ == IWLAndDPD || outputType_ == IWLOnly;
    useDPD_        = outputType_ == IWLAndDPD || outputType_ == DPDOnly;
    iwlAAIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    iwlABIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    iwlBBIntFile_  = transformationType_ == Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    nirreps_ = c->nirrep();
    nmo_     = c->ncol() + i->ncol() + a->ncol() + v->ncol();
    nso_     = i->nrow();
    sopi_    = i->rowspi();   // use i for this since there will always be occupied orbitals
    mopi_    = c->colspi() + i->colspi() + a->colspi() + v->colspi();
    clsdpi_  = i->colspi();
    openpi_  = Dimension(nirreps_); // This is the restricted constructor, there are no unpaired electrons
    frzcpi_  = c->colspi();
    frzvpi_  = v->colspi();

    // Need to smash together the C's only for them to be ripped apart elsewhere.
    std::vector<SharedMatrix > Cs;
    Cs.push_back(c); Cs.push_back(i); Cs.push_back(a); Cs.push_back(v);
    Ca_ = Matrix::horzcat(Cs);

    common_moinfo_initialize();

    if(init) initialize();
}


/**
 * Sets up the DPD buffers and performs semicanonicalization, if necessary.
 */
void
IntegralTransform::initialize()
{
    process_spaces();

    // Set up the DPD library
    // TODO implement caching of files
    int numSpaces = spacesUsed_.size();
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    cacheFiles_ = init_int_array(PSIO_MAXUNIT);
    cacheList_  = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    memory_ = Process::environment.get_memory();
    dpd_init(myDPDNum_, nirreps_, memory_, 0, cacheFiles_,
            cacheList_, NULL, numSpaces, spaceArray_);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if(transformationType_ == SemiCanonical){
        wfn_->semicanonicalize();
        Cb_ = wfn_->Cb();
    }
    process_eigenvectors();

    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);


    // Set up the correlated to Pitzer arrays.  These have to include the occupied core terms, because
    // the reference contributions are already folded into the TPDM.
    aCorrToPitzer_ = new int[nmo_];
    if(transformationType_ != Restricted){
        bCorrToPitzer_ = new int[nmo_];
    }else{
        bCorrToPitzer_ = aCorrToPitzer_;
    }
    size_t aCorrCount = 0;
    size_t bCorrCount = 0;
    size_t pitzerOffset = 0;

    // Frozen DOCC
    for(int h = 0; h < nirreps_; ++h){
        for(int n = 0; n < frzcpi_[h]; ++n){
            aCorrToPitzer_[aCorrCount++] = pitzerOffset + n;
            if(transformationType_ != Restricted)
                bCorrToPitzer_[bCorrCount++] = pitzerOffset + n;
        }
        pitzerOffset += mopi_[h];
    }
    // Active OCC
    pitzerOffset = 0;
    for(int h = 0; h < nirreps_; ++h){
        for(int n = frzcpi_[h]; n < clsdpi_[h] + openpi_[h]; ++n){
            aCorrToPitzer_[aCorrCount++] = pitzerOffset + n;
        }
        if(transformationType_ != Restricted)
            for(int n = frzcpi_[h]; n < clsdpi_[h]; ++n){
                bCorrToPitzer_[bCorrCount++] = pitzerOffset + n;
            }
        pitzerOffset += mopi_[h];
    }
    // Active VIR
    pitzerOffset = 0;
    for(int h = 0; h < nirreps_; ++h){
        for(int n = clsdpi_[h] + openpi_[h]; n < mopi_[h] - frzvpi_[h]; ++n){
            aCorrToPitzer_[aCorrCount++] = pitzerOffset + n;
        }
        if(transformationType_ != Restricted)
            for(int n = clsdpi_[h]; n < mopi_[h] - frzvpi_[h]; ++n){
                bCorrToPitzer_[bCorrCount++] = pitzerOffset + n;
            }
        pitzerOffset += mopi_[h];
    }


    initialized_ = true;
}

boost::shared_ptr<PSIO>
IntegralTransform::get_psio() const
{
    return psio_;
}

void
IntegralTransform::set_psio(boost::shared_ptr<PSIO> psio)
{
    psio_ = psio;
}

IntegralTransform::~IntegralTransform()
{
    if (initialized_) {
        dpd_close(myDPDNum_);
        free_int_matrix(cacheList_);
        free(cacheFiles_);
        free(zeros_);
    }
}

void IntegralTransform::check_initialized()
{
    if (initialized_ == false)
        throw PSIEXCEPTION("IntegralTransform::check_initialized: This instance is not initialized.");
}

