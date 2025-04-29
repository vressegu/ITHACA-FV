
#ifndef ITHACA_POD_PARAMETERS_H
#define ITHACA_POD_PARAMETERS_H

#include <map>

//#include "RedLUM.h"

#include <FieldField.H>
#include "fvCFD.H"
#include "ITHACAparameters.H"

#include "fvMesh.H"
#include "fvMeshSubset.H"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#pragma GCC diagnostic pop

#include <unsupported/Eigen/CXX11/Tensor>

namespace ITHACAPOD {


  /**
   *  @brief  Class that containing all parameters of the stochastic POD
   ***********************************************/
   #pragma message "CCR - Oublier temporairement le parent 'vide' Parameters..."
  class IthacaPODParameters {              // } : public Parameters {
  public:
    enum class PressureResolutionKind{
      FullOrder,
      ReducedOrder,
      Neglected,
      Undefined
    };
  private:
    /// Pointer to the IthacaPODParameters class for singleton pattern use
    static IthacaPODParameters* instance;// = nullptr ;

    autoPtr<argList> _args;
    autoPtr<Time> runTime0;

    Time* runTimeData;
    fvMesh* mesh;
    fvMesh* submesh;//! ref
    label nCells;

    /// Parameters from ITHACAparameters class from ITHACA-FV library
    ITHACAparameters* ithacaLibraryParameters;//! ref
    IOdictionary* ITHACAdict;
    wordList fieldlist;

    

    word eigensolver;
    int precision;
    bool exportPython;
    bool exportMatlab;
    bool exportTxt;
    std::_Ios_Fmtflags outytpe;

    bool useHypRedSto;
    /// boolean indicating if DEIM is used or not
    bool useDEIM;
    /// boolean indicating use DDES or not
    bool useDDES;
    /// boolean indicating use of DNS or not
    bool useDNS;

    /// boolean indicating use of data assimilation
    bool assimilation;

    /// boolean indicating on line reconstruction of velocity
    bool onLineReconstruct;

    /// string indicating what temporal scheme to used
    word ROMTemporalScheme;

    /// boolean indicating if the large-scale advection is corrected (LU style)
    bool advModifOrNot;
    label nBlocks;
    bool centeredOrNot;
    bool interpFieldCenteredOrNot;
    bool forcingOrNot;
    bool symDiff;
    bool forcingHasFfieldOrNot;
    bool stochasOrNot;
    bool useStratonovich;

    volVectorField* template_field_U;
    volScalarField* template_field_p;
    volScalarField* template_F_mask;
    volScalarField* template_field_nut;
    volScalarField* template_field_omega;
    volScalarField* template_field_k;

    // volScalarField* S2inflation;
    volScalarField* S2inflationMg;

    label nParticules;
    bool reducingNoisesMatrix;
    label nSimu;
    double intervalConfidenceRatio;

    HashTable<label,word> nModes;
    HashTable<word,word> hilbertSpacePOD;
    HashTable<double,word> resolvedVaryingEnergy;
    HashTable<double,word> varyingEnergy;
    HashTable<double,word> meanEnergy;
    double weightH1;
    double weightBC;
    word patchBC;
    double freqBC_0;
    List<word> field_name;
    List<word> field_type;

    fileName casenameData;

    label startTime;
    label endTime;
    label nSnapshots;

    label endTimeSimulation;
    label nSnapshotsSimulation;

    double saveTime;
    double FinalTime;
    double FinalTimeSimulation;
    label nSubsample;
    dimensionedScalar corTime;
    label interval_subsample;
    double time_interval_subsample;

    dimensionedScalar* nu;
    dimensionedScalar* freqRotActuDisk;
    dimensionedVector* centerActuDisk;
    dimensionedScalar* lengthActuDisk;
    dimensionedScalar* widthActuDisk;
    dimensionedScalar* initPhaseActuDisk;
    dimensionedScalar* timeBeforeRota;
    dimensionedScalar* magForcage;
    dimensionedScalar* simu;

    PressureResolutionKind pressureResolutionKind;

    Eigen::VectorXd eigenValues_U;
    Eigen::VectorXd eigenValues_p;
    Eigen::VectorXd lambda;

    /// Turbulent viscosity coefficient for Smagorinsky model
    /// Standard values used in OpenFOAM
    /// see : A. Yoshizawa, “Statistical theory for compressible turbulent shear flows, withthe application to subgrid modeling,” Phys. Fluids29, 2152 (1986).
    float Ck;
    float Ce;

    /// Number of magic points for build_DEIM_U
    label nMagicPoints;

    /// DEIM parameters
    PtrList<volTensorField> deformationTensorOfModesOnMagicPoints;
    PtrList<volTensorField> deformationTensorOfModesOnMagicNeighborhoods;
    PtrList<volVectorField> tracerGradOfModesOnMagicPoints;
    PtrList<volVectorField> tracerGradOfModesOnMagicNeighborhoods;
    List<label> magicPoints;
    List<label> localMagicPoints;

    /// DEIM matrix (high dimensioned)
    Eigen::MatrixXd K_DEIM_matrix;
    Eigen::MatrixXd projected_K_DEIM_matrix;
    Eigen::MatrixXd projected_MK_DEIM_matrix;

    /// Geometric Information
    volScalarField* volume;
    double totalVolume;
    volScalarField* delta;
    volScalarField* magicDelta;

    volVectorField* meanU;
    volVectorField* meanVectorDEIM;
    volScalarField* meanScalarDEIM;
    volVectorField* meanVectorDEIMMagic;
    volScalarField* meanScalarDEIMMagic;

    word DEIMInterpolatedField;
    word folder_DEIM;
    bool inflatNut;
    /// for DEIM on nut
    Eigen::Tensor<double, 3> xi_onMagicPts;
    Eigen::Tensor<double, 3 > M_xi_onMagicPts;

    Eigen::VectorXd* deltaWeight;
  private:
    
  public:
    IthacaPODParameters(int argc, char* argv[]);
    static IthacaPODParameters* getInstance();
    ~IthacaPODParameters();

    template<typename T>
    void read_snapshot(T& snapshot, const label& i_snap, word path = "default_path", word name = "default_name") const;

    /// \brief Getter to access all IthacaPODParameters members in the other classes
    ITHACAparameters* get_ITHACAparameters() const { return ithacaLibraryParameters; }

    word resolutionType;

      /// \brief Word indicating if we use SOTA or not and if yes which case of SOTA
      /// "0" means no SOTA
      /// "D" means deterministic SOTA
      /// "S" means stochastic SOTA
    word useSOTA;

    const wordList& get_fieldlist() const { return fieldlist; }
    const HashTable<label,word>& get_nModes() const { return nModes; }
    const HashTable<double,word>& get_resolvedVaryingEnergy() const { return resolvedVaryingEnergy; }
    const HashTable<double,word>& get_varyingEnergy() const { return varyingEnergy; }
    const HashTable<double,word>& get_meanEnergy() const { return meanEnergy; }
    const HashTable<word,word>& get_hilbertSpacePOD() const { return hilbertSpacePOD; } //Contains both the field name and Hilbert space
    word get_pathHilbertSpace(word fieldName);
    word get_pathHilbertSpace_fromHS(word hilbertSp);
    const double& get_weightH1() const { return weightH1; }
    const double& get_weightPOD() const { return weightBC; }
    const word& get_patchBC() const { return patchBC; }
    const double& get_freqBC() const { return freqBC_0; }
    const List<word>& get_field_name() const { return field_name; }
    const List<word>& get_field_type() const { return field_type; }

    // const fileName& get_pathData() const { return pathData; }
    const fileName& get_casenameData() const { return casenameData; }

    const dimensionedScalar& get_freqRotActuDisk() const { return *freqRotActuDisk; }
    const dimensionedVector& get_centerActuDisk() const { return *centerActuDisk; }
    const dimensionedScalar& get_lengthActuDisk() const { return *lengthActuDisk; }
    const dimensionedScalar& get_widthActuDisk() const { return *widthActuDisk; }
    const dimensionedScalar& get_initPhaseActuDisk() const { return *initPhaseActuDisk; }
    const dimensionedScalar& get_timeBeforeRota() const { return *timeBeforeRota; }
    const dimensionedScalar& get_magForcage() const { return *magForcage; }
    const dimensionedScalar& get_simu() const { return *simu; }

    // getting current F_0 force vector and F_mask
    const Foam::Vector<double>& get_F_0() const;
    void compute_F_mask();
    const volScalarField& get_F_mask() const { return *template_F_mask; }

    const PressureResolutionKind& get_pressureResolutionKind(){return pressureResolutionKind;}

    const label& get_startTime() const { return startTime; }
    const label& get_endTime() const { return endTime; }
    const label& get_nSnapshots() const { return nSnapshots; }
    const word& get_ROMTemporalScheme() const { return ROMTemporalScheme; }

    const double& get_initialTimeSimulation() const {
      return FinalTime; }
    const double& get_finalTime() const {
      return FinalTime; }
    const double& get_finalTimeSimulation() const {
      return FinalTimeSimulation; }

    const label& get_endTimeSimulation() const { return endTimeSimulation; }
    const label& get_nSnapshotsSimulation() const { return nSnapshotsSimulation; }

    const label& get_nBlocks() const { return nBlocks; }

    const bool& get_centeredOrNot() const { return centeredOrNot; }
    const bool& get_interpFieldCenteredOrNot() const { return interpFieldCenteredOrNot; }

    const bool& get_exportPython() const { return exportPython; }
    const bool& get_exportMatlab() const { return exportMatlab; }
    const bool& get_exportTxt() const { return exportTxt; }
    const int& get_precision() const { return precision; }
    const std::_Ios_Fmtflags& get_outytpe() const { return outytpe; }

    const word& get_eigensolver() const { return eigensolver; }
    const word& get_useSOTA() const {return useSOTA;}


    // const Time* get_runTime() const { return runTimeData; }
    // const Time& get_runTime() const { return *runTimeData; }
    const fvMesh& get_mesh() const { return *mesh; }
    const label& get_nCells() const { return nCells; }

    const volVectorField& get_template_field_U()   const { return *template_field_U; }
    const volScalarField& get_template_field_p()   const { return *template_field_p; }
    const volScalarField& get_template_field_nut() const { return *template_field_nut; }
    const volScalarField& get_template_field_omega()   const { return *template_field_omega; }
    const volScalarField& get_template_field_k()   const { return *template_field_k; }

    const double& get_saveTime() const { return saveTime; }

    const dimensionedScalar& get_nu() const { return *nu; }

    const label& get_nSubsample() const { return nSubsample; }
    const dimensionedScalar& get_corTime() const { return corTime; }
    const label& get_interval_subsample() const { return interval_subsample; }
    const double& get_time_interval_subsample() const { return time_interval_subsample; }

    const Eigen::VectorXd& get_eigenValues_U() const { return eigenValues_U; }
    const Eigen::VectorXd& get_eigenValues_p() const { return eigenValues_p; }
    const Eigen::VectorXd& get_lambda() const { return lambda; }

    const label& get_nParticules() const { return nParticules; }
    const bool& get_reducingNoisesMatrix() const { return reducingNoisesMatrix; }
    const label& get_nSimu() const { return nSimu; }
    const double& get_intervalConfidenceRatio() const { return intervalConfidenceRatio; }
    const bool& get_advModifOrNot() const {return advModifOrNot;}
    const bool& get_forcingOrNot() const {return forcingOrNot;}
    const bool& get_symDiff() const {return symDiff;}
    const bool& get_forcingHasFfieldOrNot() const {return forcingHasFfieldOrNot;}
    const bool& get_stochasOrNot() const {return stochasOrNot;}
    const bool& get_useStratonovich() const {return useStratonovich;}

    /// assimilation getter
    const bool& get_assimilation() const {return assimilation;}

    /// on line velocity Reconstruction getter
    const bool& get_onLineReconstruct() const {return onLineReconstruct;}

    /// DDES getter
    const bool& get_useDDES() const {return useDDES;}

    /// DNS getter
    const bool& get_useDNS() const {return useDNS;}

    /// Getters of DEIM parameters
    const bool& get_useHypRedSto() const {return useHypRedSto;}
    const bool& get_inflatNut() const {return inflatNut;}
    const bool& get_useDEIM() const {return useDEIM;}
    const float& get_Ck() const {return Ck;}
    const float& get_Ce() const {return Ce;}
    const label& get_nMagicPoints() const {return nMagicPoints;}
    const word& get_folder_DEIM() const {return folder_DEIM;}

    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicNeighborhoods() const {return tracerGradOfModesOnMagicNeighborhoods;}
    const PtrList<volVectorField>& get_tracerGradOfModesOnMagicPoints() const {return tracerGradOfModesOnMagicPoints;}
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicNeighborhoods() const {return deformationTensorOfModesOnMagicNeighborhoods;}
    const PtrList<volTensorField>& get_deformationTensorOfModesOnMagicPoints() const {return deformationTensorOfModesOnMagicPoints;}
    const List<label>& get_magicPoints() const {return magicPoints;}
    const List<label>& get_localMagicPoints() const {return localMagicPoints;}
    const Eigen::MatrixXd& get_projected_K_DEIM() const {return projected_K_DEIM_matrix;}
    const Eigen::MatrixXd& get_projected_MK_DEIM() const {return projected_MK_DEIM_matrix;}
    const volScalarField& get_delta() const {return *delta;}
    const volScalarField& get_magicDelta() const {return *magicDelta;}
    const fvMesh& get_submesh() const {return *submesh;}
    const volScalarField& get_volume() const {return *volume;}
    const double& get_totalVolume() const {return totalVolume;}
    const volVectorField& get_meanU() const {return *meanU;}
    const volVectorField& get_meanVectorDEIM() const { return *meanVectorDEIM; }
    const volScalarField& get_meanScalarDEIM() const { return *meanScalarDEIM; }
    void get_meanDEIM(volVectorField& meanVectorDEIM_in) {meanVectorDEIM_in = *meanVectorDEIM;}
    void get_meanDEIM(volScalarField& meanScalarDEIM_in) {meanScalarDEIM_in = *meanScalarDEIM;}
    const volVectorField& get_meanVectorDEIMMagic() const { return *meanVectorDEIMMagic; }
    const volScalarField& get_meanScalarDEIMMagic() const { return *meanScalarDEIMMagic; }
    // const volScalarField& get_S2inflation()   const { return *S2inflation; }
    const volScalarField& get_S2inflationMg()   const { return *S2inflationMg; }

    const word& get_DEIMInterpolatedField() const {return DEIMInterpolatedField;}

    const Eigen::Tensor<double, 3>& get_xi_onMagicPts() const {return xi_onMagicPts;}
    const Eigen::Tensor<double, 3>& get_M_xi_onMagicPts() const {return M_xi_onMagicPts;}

    Eigen::VectorXd get_deltaWeight() const {return *deltaWeight;}

    /// Getter and setter of DEIM matrix K
    const Eigen::MatrixXd& get_K_DEIM() const {return K_DEIM_matrix;}
    void set_K_DEIM(Eigen::MatrixXd& K_input) {K_DEIM_matrix=K_input;}


    /// Setters to the members of IthacaFVParameters that cannot initialize at the beginning of the program
    /// For instance, eigen values can be known only after POD was performed
    void set_nSnapshots(const label& input_nSnapshots) { nSnapshots = input_nSnapshots; }
    void set_nSnapshotsSimulation(const label& input_nSnapshotsSimulation) { nSnapshotsSimulation = input_nSnapshotsSimulation; }
    void set_eigenValues_U(const Eigen::VectorXd& input_eigenValues_U) { eigenValues_U = input_eigenValues_U; }
    void set_eigenValues_p(const Eigen::VectorXd& input_eigenValues_p) { eigenValues_p = input_eigenValues_p; }
    void set_lambda(const Eigen::VectorXd& input_lambda) { lambda = input_lambda; }
    void set_saveTime(const double& input_saveTime) { saveTime = input_saveTime; }
    void set_nSubsample(const label& input_nSubsample) { nSubsample = input_nSubsample; }
    void set_corTime(const dimensionedScalar& input_corTime) { corTime = input_corTime; }
    void set_interval_subsample(const label& input_interval_subsample) { interval_subsample = input_interval_subsample; }
    void set_time_interval_subsample(const double& input_time_interval_subsample) { time_interval_subsample = input_time_interval_subsample; }
    void set_xi_onMagicPts(const Eigen::Tensor<double, 3> xi_input) {xi_onMagicPts=xi_input;}
    void set_M_xi_onMagicPts(const Eigen::Tensor<double, 3 > M_xi_input) {M_xi_onMagicPts=M_xi_input;}
    void set_advModifOrNot(const bool& a) {advModifOrNot=a;}
    void set_useSOTA(word& input_useSOTA) {useSOTA=input_useSOTA;}

    /// set DDES
    void set_useDDES(const bool& a) {useDDES=a;}

    /// set DNS
    void set_useDNS(const bool & c) {useDNS= c;}

    /// Setters of POD and DEIM parameters
    void set_nModes(const word& field_name, const label& n)
    {
      nModes.insert(field_name,n);
    }
    void set_resolvedVaryingEnergy(const word& field_name, const double& n)
    {
      resolvedVaryingEnergy.set(field_name,n);
    }
    void set_varyingEnergy(const word& field_name, const double& n)
    {
      varyingEnergy.set(field_name,n);
    }
    void set_meanEnergy(const word& field_name, const double& n)
    {
      meanEnergy.insert(field_name,n);
    }
    void set_hilbertSpacePOD(const word& field_name, const word& hilbertSp)
    {
      hilbertSpacePOD.insert(field_name,hilbertSp);
    }
    void set_weightH1(const double & c) {weightH1= c;}
    void set_nMagicPoints(const label& n) {nMagicPoints=n;}
    void set_Ck(const float& c) {Ck=c;}
    void set_Ce(const float& c) {Ce=c;}
    void set_useDEIM(const bool& b) {useDEIM=b;}

    void set_useHypRedSto(const bool& x) {useHypRedSto=x;}  //hyperreduction
    void set_inflatNut(const bool& x) {inflatNut=x;}

    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volVectorField>& defT) {
      tracerGradOfModesOnMagicNeighborhoods.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        tracerGradOfModesOnMagicNeighborhoods.set(k, new volVectorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volVectorField>& defT) {
      tracerGradOfModesOnMagicPoints.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        tracerGradOfModesOnMagicPoints.set(k, new volVectorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicNeighborhoods(PtrList<volTensorField>& defT) {
      deformationTensorOfModesOnMagicNeighborhoods.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        deformationTensorOfModesOnMagicNeighborhoods.set(k, new volTensorField(defT[k]));
      }
    }

    void set_deformationTensorOfModesOnMagicPoints(PtrList<volTensorField>& defT) {
      deformationTensorOfModesOnMagicPoints.resize( nModes["U"]+1 );
      for (label k = 0 ; k < (nModes["U"]+1) ; k++)
      {
        deformationTensorOfModesOnMagicPoints.set(k, new volTensorField(defT[k]));
      }
    }

    void set_magicPoints(const List<label>& nn) {magicPoints = nn;}
    void set_localMagicPoints(const List<label>& nn) {localMagicPoints = nn;}
    void set_projected_K_DEIM(const Eigen::MatrixXd& K_input) {projected_K_DEIM_matrix=K_input;}
    void set_projected_MK_DEIM(const Eigen::MatrixXd& proj_MK_input) {projected_MK_DEIM_matrix= proj_MK_input;}
    void set_magicDelta(const volScalarField& mD) {magicDelta = new volScalarField(mD);}
    void set_submesh(fvMesh& s) {submesh = &s;}
    void set_meanU(const volVectorField& mm) {meanU = new volVectorField(mm);}
    void set_meanDEIM(const volVectorField& mm) {meanVectorDEIM = new volVectorField(mm);}
    void set_meanDEIM(const volScalarField& mm) {meanScalarDEIM = new volScalarField(mm);}
    void set_meanDEIMMagic(const volVectorField& mm)
    {meanVectorDEIMMagic = new volVectorField(mm);}
    void set_meanDEIMMagic(const volScalarField& mm)
    {meanScalarDEIMMagic = new volScalarField(mm);}
    // void set_S2inflation(const volScalarField& mm)
    //   {S2inflation = new volScalarField(mm);}
    void set_S2inflationMg(const volScalarField& mm)
    {S2inflationMg = new volScalarField(mm);}

    void set_deltaWeight(Eigen::VectorXd dw){deltaWeight = new Eigen::VectorXd(dw);}

  };

}

#endif
