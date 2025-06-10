/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Class
    ITHACAparameters
Description
    Some parameters for input output informations
SourceFiles
    ITHACAparameters.H
\*---------------------------------------------------------------------------*/
#ifndef Parameters_H
#define Parameters_H

#include <iostream>
#include "fvCFD.H"
#include "ITHACAassert.H"


//CHC : adding new includes for integrating RedLUM POD calculation parameters here
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

class Parameters    {public: virtual ~Parameters()    {}}; // only one derivatenamespace ITHACAPOD {


//namespace ITHACA {
/// Singleton Class for the definition of some general parameters, the parameters must be defined from the file ITHACAdict inside the
/// system folder.
#pragma message "CCR - adding temporary empty class  Parameters..."
class ITHACAparameters : public Parameters 
{
    private:

        
        //CHC : remove old constructor...  ITHACAparameters(fvMesh& mesh, Time& localTime, const word& myfield_name="");//ITHACAparameters(fvMesh& mesh, Time& localTime);

        //CHC remove tempo ITHACAparameters(int argc, char **argv, Foam::fvMesh&  mesh, Time& localTime, const word& myfield_name="");//CHC new POD specific constructor 
        ITHACAparameters(Foam::fvMesh&  mesh, Time& localTime, const word& myfield_name="");//CHC new POD specific constructor 
        
    public:
        /// type of eigensolver used in the eigenvalue decomposition can be either be eigen or spectra
        word eigensolver;

        /// precision of the output Market Matrix objects (i.e. reduced matrices, eigenvalues, ...)
        label precision;

        ///
        bool exportPython;
        bool exportMatlab;
        bool exportTxt;
        bool debug;
        bool warnings;

        /// type of output format can be fixed or scientific
        std::_Ios_Fmtflags outytpe;

        /// Dictionary for input objects from file
        IOdictionary* ITHACAdict;

        /// runTime defined locally
        Time& runTime;

        /// Mesh object defined locally
        fvMesh& mesh;

        /// Pointer to the ITHACAparameters singleton instance
        static ITHACAparameters* instance;
    // ----------------------------------------------------------------------
    // CHC : Adding attributes from RedLUM POD Calculation (named IthcaFVParameters, then IthacaPODParameters)
    //
    enum class PressureResolutionKind{
      FullOrder,
      ReducedOrder,
      Neglected,
      Undefined
    };


    word resolutionType;

      /// \brief Word indicating if we use SOTA or not and if yes which case of SOTA
      /// "0" means no SOTA
      /// "D" means deterministic SOTA
      /// "S" means stochastic SOTA
    word useSOTA;

  private:
    autoPtr<Time> runTime0;

    Time* runTimeData;
    fvMesh* submesh;//! ref
    label nCells;

    wordList fieldlist;

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

    word myfield_name;
    //
    // CHC : END OF Adding attributes from RedLUM POD Calculation (named IthcaFVParameters, then IthacaPODParameters)
    // ----------------------------------------------------------------------

    public:
        ///
        /// @brief      Gets an instance of ITHACAparameters, to be used if the instance is not existing.
        ///
        /// @param      mesh       The mesh
        /// @param      localTime  The runtime object
        ///
        /// @return     The instance.
        ///
        static ITHACAparameters* getInstance(fvMesh& mesh, Time& localTime, const word& myfield_name="");

        ///
        /// @brief      Gets an instance of ITHACAparameters, to be used if the instance is already existing.
        ///
        /// @return     The instance.
        ///
        static ITHACAparameters* getInstance();


        /// Delete empty constructor
        //~ITHACAparameters() = delete;
        #pragma message "recreate the destructor (after heriting Parameter empty class for RedLUM compatibility)"
        ~ITHACAparameters(){}

    public:
    // ----------------------------------------------------------------------
      template<typename T>
    void read_snapshot(T& snapshot, const label& i_snap, word path = "default_path", word name = "default_name") const;




  // CHC : Adding getters and setters for new POD calculation attributes (from RedLUM )
    //
    autoPtr<Time> get_runTime0(){return runTime0;}
    Foam::Time& get_ref_runTime0(){ return *runTime0;}

    wordList& get_fieldlist() { return fieldlist; }
    HashTable<label,word>& get_nModes() { return nModes; }
    HashTable<double,word>& get_resolvedVaryingEnergy() { return resolvedVaryingEnergy; }
    HashTable<double,word>& get_varyingEnergy() { return varyingEnergy; }
    HashTable<double,word>& get_meanEnergy() { return meanEnergy; }
    HashTable<word,word>& get_hilbertSpacePOD() { return hilbertSpacePOD; } //Contains both the field name and Hilbert space
/* //CHC ORIG 
    const wordList& get_fieldlist() const { return fieldlist; }
    const HashTable<label,word>& get_nModes() const { return nModes; }
    const HashTable<double,word>& get_resolvedVaryingEnergy() const { return resolvedVaryingEnergy; }
    const HashTable<double,word>& get_varyingEnergy() const { return varyingEnergy; }
    const HashTable<double,word>& get_meanEnergy() const { return meanEnergy; }
    const HashTable<word,word>& get_hilbertSpacePOD() const { return hilbertSpacePOD; } //Contains both the field name and Hilbert space
*/



    word get_pathHilbertSpace(word fieldName);
    word get_pathHilbertSpace_fromHS(word hilbertSp);
    double& get_weightH1() { return weightH1; }
    double& get_weightPOD() { return weightBC; }
    word& get_patchBC() { return patchBC; }
    double& get_freqBC() { return freqBC_0; }
    List<word>& get_field_name() { return field_name; }
    List<word>& get_field_type() { return field_type; }
/*  //CHC ORIG 
    const double& get_weightH1() const { return weightH1; }
    const double& get_weightPOD() const { return weightBC; }
    const word& get_patchBC() const { return patchBC; }
    const double& get_freqBC() const { return freqBC_0; }
    const List<word>& get_field_name() const { return field_name; }
    const List<word>& get_field_type() const { return field_type; }
*/



    // const fileName& get_pathData() const { return pathData; }
    fileName& get_casenameData() { return casenameData; }

    dimensionedScalar* get_pfreqRotActuDisk() { return freqRotActuDisk; }
    dimensionedVector* get_pcenterActuDisk() { return centerActuDisk; }
    dimensionedScalar* get_plengthActuDisk() { return lengthActuDisk; }
    dimensionedScalar* get_pwidthActuDisk() { return widthActuDisk; }
    dimensionedScalar* get_pinitPhaseActuDisk() { return initPhaseActuDisk; }
    dimensionedScalar* get_ptimeBeforeRota() { return timeBeforeRota; }
    dimensionedScalar* get_pmagForcage() { return magForcage; }
    dimensionedScalar* get_psimu() { return simu; }

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
    Foam::Vector<double>& get_F_0() ;
    void compute_F_mask();
    volScalarField& get_F_mask() { return *template_F_mask; }
    volScalarField* get_pF_mask() { return template_F_mask; }

    ITHACAparameters::PressureResolutionKind& get_pressureResolutionKind(){return pressureResolutionKind;}
    //PressureResolutionKind* get_ppressureResolutionKind(){return &pressureResolutionKind;}

    label& get_startTime() { return startTime; }
    label& get_endTime() { return endTime; }
    label& get_nSnapshots() { return nSnapshots; }
    word& get_ROMTemporalScheme() { return ROMTemporalScheme; }

    double& get_initialTimeSimulation() {
      return FinalTime; }
    double& get_finalTime() {
      return FinalTime; }
    double& get_finalTimeSimulation() {
      return FinalTimeSimulation; }

    label& get_endTimeSimulation() { return endTimeSimulation; }
    label& get_nSnapshotsSimulation() { return nSnapshotsSimulation; }

    label& get_nBlocks() { return nBlocks; }
/*  //CHC ORIG 
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
*/








const bool& get_centeredOrNot() const { return centeredOrNot; }
    const bool& get_interpFieldCenteredOrNot() const { return interpFieldCenteredOrNot; }

    const bool& get_exportPython() const { return exportPython; }
    const bool& get_exportMatlab() const { return exportMatlab; }
    const bool& get_exportTxt() const { return exportTxt; }
    const int& get_precision() const { return precision; }
    const std::_Ios_Fmtflags& get_outytpe() const { return outytpe; }

    const word& get_eigensolver() const { return eigensolver; }
    const word& get_useSOTA() const {return useSOTA;}


     Time* get_runTime() { return runTimeData; }    //CHC ORIG const Time* get_runTime() const { return runTimeData; }
     
    // const Time& get_runTime() const { return *runTimeData; }
    
    //CHC : mesh was a pointer in the RedLUM Parameters, but become a reference in Ithaca-FV    const fvMesh& get_mesh() const { return *mesh; }
    fvMesh& get_mesh() { return mesh; }//const fvMesh& get_mesh() const { return mesh; }
    fvMesh* get_pmesh() { return &mesh; }         //CHC ORIG const fvMesh* get_pmesh() const { return &mesh; }

    const label& get_nCells() const { return nCells; }

    const volVectorField& get_template_field_U()   const { return *template_field_U; }
    const volScalarField& get_template_field_p()   const { return *template_field_p; }
    const volScalarField& get_template_field_nut() const { return *template_field_nut; }
    const volScalarField& get_template_field_omega()   const { return *template_field_omega; }
    const volScalarField& get_template_field_k()   const { return *template_field_k; }
    volVectorField* get_ptemplate_field_U()       { return template_field_U; }
    volScalarField* get_ptemplate_field_p()       { return template_field_p; }
    volScalarField* get_ptemplate_field_nut()     { return template_field_nut; }
    volScalarField* get_ptemplate_field_omega()   { return template_field_omega; }
    volScalarField* get_ptemplate_field_k()       { return template_field_k; }

    const double& get_saveTime() const { return saveTime; }

    const dimensionedScalar& get_nu() const { return *nu; }
    dimensionedScalar* get_pnu() { return nu; }

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
    volScalarField* get_pdelta() {return delta;}
    const volScalarField& get_magicDelta() const {return *magicDelta;}
    volScalarField* get_pmagicDelta() {return magicDelta;}
    fvMesh& get_submesh() {return *submesh;}//CHC ORIG const fvMesh& get_submesh() const {return *submesh;}
    fvMesh* get_psubmesh() {return submesh;}//CHC ORIG const fvMesh& get_submesh() const {return *submesh;}
    const volScalarField& get_volume() const {return *volume;}
    volScalarField* get_pvolume() {return volume;}
    const double& get_totalVolume() const {return totalVolume;}
    const volVectorField& get_meanU() const {return *meanU;}
    volVectorField* get_pmeanU() {return meanU;}
    const volVectorField& get_meanVectorDEIM() const { return *meanVectorDEIM; }
    const volScalarField& get_meanScalarDEIM() const { return *meanScalarDEIM; }
    volVectorField* get_pmeanVectorDEIM() { return meanVectorDEIM; }
    volScalarField* get_pmeanScalarDEIM() { return meanScalarDEIM; }
    void get_meanDEIM(volVectorField& meanVectorDEIM_in) {meanVectorDEIM_in = *meanVectorDEIM;}
    void get_meanDEIM(volScalarField& meanScalarDEIM_in) {meanScalarDEIM_in = *meanScalarDEIM;}
    const volVectorField& get_meanVectorDEIMMagic() const { return *meanVectorDEIMMagic; }
    volVectorField* get_pmeanVectorDEIMMagic() { return meanVectorDEIMMagic; }
    const volScalarField& get_meanScalarDEIMMagic() const { return *meanScalarDEIMMagic; }
    volScalarField* get_pmeanScalarDEIMMagic() { return meanScalarDEIMMagic; }
    // const volScalarField& get_S2inflation()   const { return *S2inflation; }
    const volScalarField& get_S2inflationMg()   const { return *S2inflationMg; }
    volScalarField* get_pS2inflationMg()    { return S2inflationMg; }

    const word& get_DEIMInterpolatedField() const {return DEIMInterpolatedField;}

    const Eigen::Tensor<double, 3>& get_xi_onMagicPts() const {return xi_onMagicPts;}
    const Eigen::Tensor<double, 3>& get_M_xi_onMagicPts() const {return M_xi_onMagicPts;}

    Eigen::VectorXd get_deltaWeight() const {return *deltaWeight;}
    Eigen::VectorXd* get_pdeltaWeight() {return deltaWeight;}

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

    //
    // CHC : END OF Adding getters and setters for new POD calculation attributes (from RedLUM )
    // ----------------------------------------------------------------------

};
//}
#endif
