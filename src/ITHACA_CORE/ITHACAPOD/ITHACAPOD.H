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
    ITHACAPOD
Description
    Implementation of a POD decomposition of a general field
SourceFiles
    ITHACAPOD.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ITHACAPOD class.
/// \dir
/// Directory containing the header and source files for the ITHACAPOD class.

#ifndef ITHACAPOD_H
#define ITHACAPOD_H

#include "fvCFD.H"
#include "ITHACAutilities.H"
#include "ITHACAstream.H"
#include "ITHACAparameters.H"
#include "Foam2Eigen.H"
#include "EigenFunctions.H"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>
#pragma GCC diagnostic pop


#include "IthacaPODParameters.H"

/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// namespace for the computation of the POD, it exploits bot the method of snapshots and SVD.
namespace ITHACAPOD
{

/* => transfered to IthacaPODTemplate.H 
    /// \brief Structure to store the indexes of the added Triangular block to the cov matrix
    /// Triangular block
    struct indexTri
    {
        label index_start;
        label index_end;
    };


    /// \brief Structure to store the indexes of the added Square block to the cov matrix
    /// Square block
    struct indexSquare
    {
        label index1_start;
        label index1_end;
        label index2_start;
        label index2_end;
    };
*/



/* => transfered to IthacaPODTemplate.H 
template<typename T>
class ITHACAPODTemplate
{
    // CCR : TODO : think to a thread-safe method for this parameter class
//protected:
public:
    ITHACAparameters* ithacaFVParameters;
    IthacaPODParameters* ithacaPODParams;

    /// @brief Path where OpenFOAM simulation outputs are
    word field_name;
    /// @brief Type of data used in the POD
    const fileName casenameData;
    /// @brief Number of snapshot, also number of time step
    label l_nSnapshot;
    /// @brief number of blocks to load snapshots
    label l_nBlocks;
    /// @brief number of modes
    label l_nmodes;
    /// @brief Hilbert space used for the POD scalar product
    /// (and the associated norm to sort the modes)
    word l_hilbertSp;
    /// @brief Weight of the weighted H1 scalar product
    double weightH1;
    /// @brief Weight of the boundary conditions
    double weightBC;
    /// @brief Patch for the boundary conditions
    word patchBC;
    /// @brief start time of the POD decomposition
    label l_startTime;
    /// @brief end time of the POD decomposition
    label l_endTime;

    /// @brief Number of Simulation snapshot, also number of time step
    label l_nSnapshotSimulation;
    /// @brief end time of the Simulation
    label l_endTimeSimulation;

    /// @brief name used for save the covMatrix
    fileName name_covMatrix;
    /// @brief folder used to save the covMatrix
    word folder_covMatrix;
    /// @brief covMatrix existence boolean
    bool exist_covMatrix;

    // name and folder to save eigenValues et eigenVectors
    fileName name_eigenValues;
    fileName name_eigenValuesNormalized;
    fileName name_cumEigenValues;
    fileName name_eigenVector;
    /// @brief folder used to save the eigenDecomposition
    word folder_eigen;
    /// @brief eigenDecomposition existence boolean
    bool exist_eigenDecomposition;

    /// @brief folder of the spatialModes
    word folder_spatialModes;
    /// @brief spatialModes existence boolean
    bool exist_spatialModes;

    /// @brief folder of the temporal Modes
    word folder_temporalModes;
    /// @brief temporalModes existence boolean
    bool exist_temporalModes;

    /// @brief folder of the temporal Modes
    word folder_temporalModesSimulation;
    /// @brief temporalModesSimulation existence boolean
    bool exist_temporalModesSimulation;

    /// @brief template of field readed, useful for the ITHACAstream::read_fields mlethod
    T *f_field;

    
    // if true, the temporal mean of snapshot is computed and subtract to snapshot fields
    // else, temporal mean is not computed and snapshot fields are not modified
    bool b_centeredOrNot;
    T *f_meanField;
    word folder_mean;
    bool exist_noMean;

/// @brief Temporal modes energies
    Eigen::VectorXd eigenValueseigLam;
    Eigen::VectorXd lambda;

    word w_eigensolver;
    int i_precision;
    std::_Ios_Fmtflags ios_outytpe;

    Foam::Time runTime2;

    Eigen::MatrixXi inletIndex;
    bool lifting;
    PtrList<T> liftfield;

public:
    //IthacaFVPODTemplate(Parameters* myParameters, const word& myfield_name);
    ITHACAPODTemplate<T>(fvMesh& mesh, Time& localTime, const word& myfield_name);
    virtual ~ITHACAPODTemplate<T>();
 
    //--------------------------------------------------------------------------
    /// \brief Homogenize the snapshot matrix, it works with PtrList of volVectorField and volScalarField
    /// origin : redLUM 
    /// @param[in]      snapshots     The list of snapshots to be homogenized.
     void lift(PtrList<T>& snapshots);

    /// \brief Homogenize the snapshot matrix, it works with volVectorField and volScalarField
    /// origin : redLUM 
    ///
    /// @param[in]      snapshot     The single snapshot to be homogenized.
    void lift(T& snapshot);



    /// \brief computation of test temporal modes energies
    /// origin : redLUM 
    /// \param[inout] temporalModes
    void compute_lambda(Eigen::MatrixXd& temporalModes);

    /// \brief method to add a new square block to the covariance matrix
    /// used in buildCovMatrix method
    /// origin : redLUM 
    /// \param[inout] : covMatrix, covariance matrix
    /// \param[in] : snapshots1 and snapshots2, sample of the field
    /// \param[in] : indexes of the added block
    void addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots1,
        PtrList<T>& snapshots2,
        ITHACAPOD::indexSquare& indSquare);
    
    

    /// \brief method to add a new triangular block to the covariance matrix (near the diagonal)
    /// used in buildCovMatrix method
    /// origin : redLUM 
    /// \param[inout] covMatrix : the covariance matrix
    /// \param[in] snapshots : sample of the field
    /// \param[in] indTri : indexes of the added block
    void addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots,
        ITHACAPOD::indexTri& indTri);
    
    /// \brief temporary covMatrix filename
    ///
    /// filename format : folder_covMatrix + name_covMatrix + "_temp_<i>_<j>.npy"
    /// origin : redLUM 
    /// \param[in] i 
    /// \param[in] j
    word nameTempCovMatrix( int& i, int& j);

    /// \brief save the temporary covMatrix file
    /// origin : redLUM 
    /// \param[inout] covMatrix : the covariance matrix
    /// \param[in] i
    /// \param[in] j
    void saveTempCovMatrix(Eigen::MatrixXd& covMatrix, int i, int j);

    /// \brief delete the temporary covMatrix file
    /// origin : redLUM 
    /// \param[in] i
    /// \param[in] j
    void deleteTempCovMatrix(int i, int j);

    /// \brief delete previous [step-N] covMatrix temp file -> keep last N covMatrix temp file
    void deletePreviousTempCovMatrix_N(int* valI, int* valJ, int i, int j, int N);

    /// \brief find a covMatrix temporary file
    /// if a covMatrix temporary file is found in the covMat folder
    /// it is loaded and the process continues where it stopped
    /// \param[inout] covMat : the covariance matrix
    /// \param[in] index1
    /// \param[in] index2
    void findTempFile(Eigen::MatrixXd* covMat, int* index1, int* index2);

    
    void appendMeanfieldtoSpatialModes(PtrList<T>& spatialModes);
    
    
    PtrList<T> computeSpatialModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);

    Eigen::MatrixXd computeTemporalModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);

    void define_paths();

    void diagonalisation(Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);

    Eigen::MatrixXd buildCovMatrix();

    void computeMeanField();

    Eigen::MatrixXd computeSimulationTemporalModes(PtrList<T>& f_spatialModes);
    
    void getModes(PtrList<T>& spatialModes,
      Eigen::MatrixXd& temporalModes, Eigen::MatrixXd& temporalModesSimulation,
      Eigen::MatrixXd& covMatrix);
        
};//end of class ITHACAPODTemplate
*/

void computeLift(PtrList<volScalarField>& Lfield, PtrList<volScalarField>& liftfield, PtrList<volScalarField>& omfield, Eigen::MatrixXi inletIndex);
void computeLift(PtrList<volVectorField>& Lfield, PtrList<volVectorField>& liftfield, PtrList<volVectorField>& omfield, Eigen::MatrixXi inletIndex);
void computeLift(PtrList<volTensorField>& Lfield, PtrList<volTensorField>& liftfield, PtrList<volTensorField>& omfield, Eigen::MatrixXi inletIndex);

word get_pathHilbertSpace_fromHS(word hilbertSp);

//end of CCR ADD




// Public Members

//------------------------------------------------------------------------------
/// Computes the bases or reads them for a field
///
/// @param[in]  snapshots   List of snapshots.
/// @param[out] modes       A PtrList where modes are stored (it must be passed
///                         empty).
/// @param[in]  fieldName   The field name
/// @param[in]  podex       If 1, the functions read the stored mode. If 0, the
///                         function computes the modes and stores them.
/// @param[in]  supex       If 1, the functions read the stored supremizer mode.
///                         If 0, the function computes and stores them.
/// @param[in]  sup         If 1 it computes the supremizer modes.
/// @param[in]  nmodes      Number of modes to be stored. If 0, the maximum
///                         number of modes will computed.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex = 0, bool sup = 0,
    label nmodes = 0, bool correctBC = true);

//------------------------------------------------------------------------------
/// @brief      Gets the bases for a scalar field using SVD instead of the
///             method of snapshots
///
/// @param[in]  snapshots   List of snapshots.
/// @param[out] modes       A PtrList where modes are stored (it must be passed
///                         empty).
/// @param[in]  fieldName   The file name
/// @param[in]  podex       If 1, the functions read the stored mode. If 0, the
///                         function computes the modes and stores them.
/// @param[in]  supex       If 1, the functions read the stored supremizer mode.
///                         If 0, the function computes and stores them.
/// @param[in]  sup         If 1 it computes the supremizer modes.
/// @param[in]  nmodes      Number of modes to be stored. If 0, the maximum
///                         number of modes will computed.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void getModesSVD(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>&
    snapshots, PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex = 0, bool sup = 0,
    label nmodes = 0, bool correctBC = true);


//------------------------------------------------------------------------------
/// Nested-POD approach. Computes the nested snapshot matrix and weighted bases
/// for a vector field
///
/// @param[in]  snapshots   List of snapshots.
/// @param      modes       The modes
/// @param[in]  fieldName   The field name
/// @param[in]  Npar        Number of parameters
/// @param[in]  NnestedOut  Number of weighted modes to be stored for the
///                         assembly of the nested matrix
/// @param[out] UModesGlobal  List of snapshots where modes are stored (it must be
///                           passed empty).
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void getNestedSnapshotMatrix(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, label Npar,
    label NnestedOut);


//------------------------------------------------------------------------------
/// Computes the weighted bases (using the nested-pod approach) or read them for
/// a vector field
///
/// @param[in]  snapshots   List of snapshots.
/// @param[out] modes       A PtrList where modes are stored (it must be passed
///                         empty).
/// @param[in]  fieldName   The field name
/// @param[in]  podex       If 1, the functions read the stored mode. If 0, the
///                         function computes the modes and stores them.
/// @param[in]  supex       If 1, the functions read the stored supremizer mode.
///                         If 0, the function computes and stores them.
/// @param[in]  sup         If 1, it computes the supremizer modes.
/// @param[in]  nmodes      Number of modes to be stored. If 0, the maximum
///                         number of modes will computed.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void getWeightedModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex = 0, bool sup = 0,
    label nmodes = 0, bool correctBC = true);

//------------------------------------------------------------------------------
/// @brief      Performs GrammSchmidt orthonormalization on an Eigen Matrix
///
/// @param[out] Matrix  The matrix
///
void GrammSchmidt(Eigen::MatrixXd& Matrix);

//------------------------------------------------------------------------------
/// Computes the correlation matrix given a vector field snapshot Matrix using
/// different norms depending on the input snapshots
///
/// @param[in]  snapshots   List of snapshots.
///
/// @tparam     Field_type  Type of the input list you are passing, it can be:
///                         PtrList<volVectorField>, PtrList<volScalarField>,
///                         List<Eigen::SparseMatrix<double>>,
///                         List<Eigen::VectorXd>.
///
/// @return     the Eigen::MatrixXd correlation matrix.
///
template<class Field_type>
Eigen::MatrixXd corMatrix(Field_type& snapshots);


//------------------------------------------------------------------------------
/// Exports the basis for an OpenFOAM GeometricField into the ITHACAOutput/POD
/// or ITHACAOutput/supremizer
///
/// @param      snapshots   The snapshots
/// @param[in]  bases       a PtrList of GeometricField where the bases are
///                         stored.
/// @param[in]  fieldName   The field name
/// @param[in]  sup         a boolean variable 1 if you want to export the
///                         supremizer bases (in ITHACAOutput/supremizer) 0
///                         elsewhere (Default is 0).
/// @param[in]  s     a PtrList of GeometricField where the snapshots associated
///                   with the bases are stored.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void exportBases(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& bases,
    word fieldName,
    bool sup = 0);

//------------------------------------------------------------------------------
/// Exports the eigenvalues as a txt file
///
/// @param[in]  Eigenvalues  a scalarField of eigenvalues.
/// @param[in]  name         a fileName to indicate the name of the txt file.
/// @param[in]  sup          If 1, it exports the eigenValues in
///                          "ITHACAOutput/supremizer". If 0, no actions
///                          (Default is 0).
///
void exportEigenvalues(scalarField Eigenvalues, fileName name,
                       bool sup = 0);

//------------------------------------------------------------------------------
/// Export the cumulative eigenvalues as a txt file
///
/// @param[in]  cumEigenvalues  a scalarField of cumEigenvalues.
/// @param[in]  name            a fileName to indicate the name of the txt file.
/// @param[in]  sup             If 1, it exports the cumEigenValues in
///                             "ITHACAOutput/supremizer". If 0, no actions
///                             (Default is 0).
///
void exportcumEigenvalues(scalarField cumEigenvalues,  fileName name,
                          bool sup = 0);

//------------------------------------------------------------------------------
/// @brief      Get the DEIM modes for a generic non linear function
///
/// @param[in]  SnapShotsMatrix  The snapshots matrix
/// @param[in]  nmodes           The number of modes
/// @param[in]  FunctionName     The function name
///
/// @tparam     Type             templated object for the snapshots matrix type
///
/// @return     Type             the POD modes
///
template<class Type, template<class> class PatchField, class GeoMesh>
PtrList<GeometricField<Type, PatchField, GeoMesh>> DEIMmodes(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& SnapShotsMatrix,
            label nmodes, word FunctionName, word FieldName);

//------------------------------------------------------------------------------
/// @brief      Get the DEIM modes for a generic non-parametrized matrix coming
///             from a differential operator function
///
/// @param[in]  MatrixList   The matrix list as a list of pointers
/// @param[in]  nmodesA      The number of modes for A
/// @param[in]  nmodesB      The number of modes for B
/// @param[in]  MatrixName   The matrix name, used to save the matrix
///
/// @tparam     type_matrix  The   type of matrix, can be fvScalarMatrix or
///                          fvVectorMatrix
///
/// @return     It returns a tuple containing the list of modes for the matrix
///             and the list of modes for the source term
///
template<typename type_matrix>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(PtrList<type_matrix>& MatrixList, label nmodesA, label nmodesB,
                  word MatrixName = "Matrix");

//------------------------------------------------------------------------------
/// @brief      Get the DEIM modes for a generic non a parametrized matrix
///             coming from a differential operator function
///
/// @param      A           List of matrices for the lhs
/// @param      b           List of vectors for the rhs
/// @param[in]  nmodesA     The number of modes
/// @param[in]  nmodesB     The nmodes b
/// @param[in]  MatrixName  The matrix name, used to save the matrix
/// @param[in]  MatrixList  The matrix list as a list of pointers
///
/// @tparam     The   type of matrix, can be fvScalarMatrix or fvVectorMatrix
///
/// @return     It returns a tuple containing the list of modes for the matrix
///             and the list of modes for the source term
///
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(List<Eigen::SparseMatrix<double>>& A, List<Eigen::VectorXd>& b,
                  label nmodesA, label nmodesB, word MatrixName = "Matrix");


//------------------------------------------------------------------------------
/// @brief      Gets the modes.
///
/// @param[in]  snapshots   List of snapshots.
/// @param[out] modes       A PtrList where modes are stored (it must be passed
///                         empty).
/// @param      Volumes     The volumes
/// @param[in]  fieldName   The field name
/// @param[in]  podex       If 1, the functions read the stored mode. If 0, the
///                         function computes the modes and stores them.
/// @param[in]  supex       If 1, the functions read the stored supremizer mode.
///                         If 0, the function computes and stores them.
/// @param[in]  sup         If 1 it computes the supremizer modes.
/// @param[in]  nmodes      Number of modes to be stored. If 0, the maximum
///                         number of modes will computed.
///
/// @tparam     Type        vector or scalar.
/// @tparam     PatchField  fvPatchField or fvsPatchField.
/// @tparam     GeoMesh     volMesh or surfaceMesh.
///
template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
              PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
              PtrList<volScalarField>& Volumes, word fieldName, bool podex, bool supex = 0,
              bool sup = 0,
              label nmodes = 0);

//------------------------------------------------------------------------------
/// @brief      Gets Modes based on eigenvalue decomposition of another Field
///
/// @param      snapshots     The snapshots of which we want to create the modes
/// @param      modes         The modes of the snapshots
/// @param      fields2       The fields that we want to use for the eigenvalue
///                           decomposition
/// @param[in]  fieldName     The field name
/// @param[in]  podex         boolean variable 1 if the modes have been already
///                           computed and stored (in this case the function is
///                           reading them) 0 elsewhere.
/// @param[in]  supex         boolean variable 1 if the supremizer modes have
///                           been already computed and stored (in this case the
///                           function is reading them) 0 elsewhere.
/// @param[in]  sup           boolean variable 1 if you want to compute the
///                           supremizer modes 0 elsewhere.
/// @param[in]  nmodes        label variable to set the number of modes to be
///                           stored, if set to 0 the maximum number of modes
///                           will computed.
///
/// @tparam     Field_type    Type of field snapshots
/// @tparam     Field_type_2  Type of the other fields
///
template<class Field_type, class Field_type_2>
void getModes(PtrList<Field_type>& snapshots, PtrList<Field_type>& modes,
              PtrList<Field_type_2>& fields2, word fieldName, bool podex, bool supex = 0,
              bool sup = 0,
              label nmodes = 0);



}//End of namespace ITHACAPOD

#endif



