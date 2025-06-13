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
/// Header file of the ITHACAPODTemplate class

#ifndef __ITHACA_POD_TEMPLATE_H__
#define __ITHACA_POD_TEMPLATE_H__

#include "ITHACAPOD.H"
/* moved to "IthacaPODParameters.H"
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
*/

#include "IthacaPODParameters.H"

/*---------------------------------------------------------------------------*\
                        Class reductionProblem Declaration
\*---------------------------------------------------------------------------*/

/// namespace for the ITHACAPODTemplate class (ITHACAPOD::ITHACAPODTemplate)
namespace ITHACAPOD
{




template<typename T>
class ITHACAPODTemplate
{
public:
    ITHACAparameters* parameters;//ithacaFVParameters;
    //CHC remove, merged    IthacaPODParameters* ithacaPODParams;

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
    //CHC IthacaFVPODTemplate(Parameters* myParameters, const word& myfield_name);
    
    ITHACAPODTemplate<T>(ITHACAparameters *parameters, fvMesh& mesh, Time& localTime, const word& myfield_name);
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

    /*----------------------------------------------------------------------------*/
    /// \brief method to add a new square block to the covariance matrix
    /// used in buildCovMatrix method
    /// origin : redLUM 
    /// \param[inout] : covMatrix, covariance matrix
    /// \param[in] : snapshots1 and snapshots2, sample of the field
    /// \param[in] : indexes of the added block
    void addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots1,
        PtrList<T>& snapshots2,
        ITHACAutilities::indexSquare& indSquare);
    
    

    /// \brief method to add a new triangular block to the covariance matrix (near the diagonal)
    /// used in buildCovMatrix method
    /// origin : redLUM 
    /// \param[inout] covMatrix : the covariance matrix
    /// \param[in] snapshots : sample of the field
    /// \param[in] indTri : indexes of the added block
    void addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
        PtrList<T>& snapshots,
        ITHACAutilities::indexTri& indTri);
    
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

    const T& get_mean() const { return (*f_meanField); }

    virtual Eigen::MatrixXd buildCovMatrix();

    void computeMeanField();

    Eigen::MatrixXd computeSimulationTemporalModes(PtrList<T>& f_spatialModes);
    
    void getModes(PtrList<T>& spatialModes,
      Eigen::MatrixXd& temporalModes, Eigen::MatrixXd& temporalModesSimulation,
      Eigen::MatrixXd& covMatrix);
        
};//end of class ITHACAPODTemplate


} //end of namespace


#endif

