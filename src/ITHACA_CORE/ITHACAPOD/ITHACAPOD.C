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

\*---------------------------------------------------------------------------*/

/// \file
/// source file for the ITHACAPOD class

#include "ITHACAPOD.H"
#include "EigenFunctions.H"

namespace ITHACAPOD
{

template<class Type, template<class> class PatchField, class GeoMesh>
void getNestedSnapshotMatrix(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& ModesGlobal,
    word fieldName,
    label Npar, label NnestedOut)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<PtrList<GeometricField<Type, PatchField, GeoMesh>>>
    SnapMatrixNested;
    label Nt = snapshots.size() / Npar;
    SnapMatrixNested.setSize(Nt);

    for (label i = 0; i < Npar; i++)
    {
        SnapMatrixNested[i].resize(Nt);

        for (label j = 0; j < Nt; j++)
        {
            SnapMatrixNested[i].set(j, snapshots[j + Nt * i].clone());
        }
    }

    List<PtrList<GeometricField<Type, PatchField, GeoMesh>>> UModesNested;
    PtrList<GeometricField<Type, PatchField, GeoMesh>>  y;
    UModesNested.setSize(Nt);

    for (label i = 0; i < Npar; i++)
    {
        getWeightedModes(SnapMatrixNested[i], UModesNested[i], 0, 0, 0,
                         NnestedOut);
    }

    for (label i = 0; i < Npar; i++)
    {
        y = UModesNested[i];

        for (label j = 0; j < y.size(); j++)
        {
            ModesGlobal.append(y[j].clone());
        }
    }
}

template void getNestedSnapshotMatrix(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& ModesGlobal,
    word fieldName, label Npar, label NnestedOut);

template void getNestedSnapshotMatrix(
    PtrList<volVectorField>& snapshots, PtrList<volVectorField>& ModesGlobal,
    word fieldName, label Npar, label NnestedOut);


template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word PODkey = "POD_" + fieldName;
    word PODnorm = para->ITHACAdict->lookupOrDefault<word>(PODkey, "L2");
    M_Assert(PODnorm == "L2" ||
             PODnorm == "Frobenius", "The PODnorm can be only L2 or Frobenius");
    Info << "Performing POD for " << fieldName << " using the " << PODnorm <<
         " norm" << endl;

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        if (para->eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size() - 2;
            }

            M_Assert(nmodes <= snapshots.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }
        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size();
            }

            M_Assert(nmodes <= snapshots.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd _corMatrix;

        if (PODnorm == "L2")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots);
        }
        else if (PODnorm == "Frobenius")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots, 0, false);
        }

        if (Pstream::parRun())
        {
            reduce(_corMatrix, sumOp<Eigen::MatrixXd>());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             fieldName << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().array().reverse();
        }

        if (eigenValueseig.array().minCoeff() < 0)
        {
            eigenValueseig = eigenValueseig.array() + 2 * abs(
                                 eigenValueseig.array().minCoeff());
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        //Eigen::VectorXd eigenValueseigLam =
        //    eigenValueseig.real().array().abs().cwiseInverse().sqrt() ;
        //Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
        //                           eigenValueseigLam.head(nmodes).asDiagonal();
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig);
        // Computing Normalization factors of the POD Modes
        Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(snapshots[0]);
        Eigen::MatrixXd normFact(nmodes, 1);

        for (label i = 0; i < nmodes; i++)
        {
            if (PODnorm == "L2")
            {
                normFact(i, 0) = std::sqrt((modesEig.col(i).transpose() * V.asDiagonal() *
                                            modesEig.col(i))(0, 0));

                if (Pstream::parRun())
                {
                    normFact(i, 0) = (modesEig.col(i).transpose() * V.asDiagonal() *
                                      modesEig.col(i))(0, 0);
                }
            }
            else if (PODnorm == "Frobenius")
            {
                normFact(i, 0) = std::sqrt((modesEig.col(i).transpose() * modesEig.col(i))(0,
                                           0));

                if (Pstream::parRun())
                {
                    normFact(i, 0) = (modesEig.col(i).transpose() * modesEig.col(i))(0, 0);
                }
            }
        }

        if (Pstream::parRun())
        {
            reduce(normFact, sumOp<Eigen::MatrixXd>());
        }

        if (Pstream::parRun())
        {
            normFact = normFact.cwiseSqrt();
        }

        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig);
        }

        std::cout << normFact << std::endl;

        for (label i = 0; i < nmodes; i++)
        {
            modesEig.col(i) = modesEig.col(i).array() / normFact(i, 0);

            for (label j = 0; j < NBC; j++)
            {
                modesEigBC[j].col(i) = modesEigBC[j].col(i).array() / normFact(i, 0);
            }
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh>  tmp2(snapshots[0].name(),
                    snapshots[0]);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec, correctBC);

            for (label k = 0; k < NBC; k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModes(
    PtrList<volVectorField>& snapshots, PtrList<volVectorField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

template void getModes(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

template void getModes(
    PtrList<surfaceScalarField>& snapshots, PtrList<surfaceScalarField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

template<class Type, template<class> class PatchField, class GeoMesh>

void getWeightedModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if (nmodes == 0)
    {
        nmodes = snapshots.size() - 2;
    }

    M_Assert(nmodes <= snapshots.size() - 2,
             "The number of requested modes cannot be bigger than the number of Snapshots - 2");

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd _corMatrix = ITHACAutilities::getMassMatrix(snapshots);
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             snapshots[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                es(&op, nmodes, ncv);

        if (para->eigensolver == "spectra")
        {
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse();
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::VectorXd eigenValueseigWeigted = eigenValueseig.head(
                nmodes).real().array() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.head(nmodes).asDiagonal() *
                                   eigenValueseigWeigted.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal() * eigenValueseigWeigted.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp2(snapshots[0].name(),
                    snapshots[0]);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec, correctBC);

            for (label k = 0; k < tmp2.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        //exportBases(modes, snapshots, sup);
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getWeightedModes(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

template void getWeightedModes(
    PtrList<volVectorField>& snapshots, PtrList<volVectorField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);


template<class Type, template<class> class PatchField, class GeoMesh>
void getModesSVD(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        PtrList<volVectorField> Bases;
        modes.resize(nmodes);
        Info << "####### Performing POD using Singular Value Decomposition for " <<
             snapshots[0].name() << " #######" << endl;
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(snapshots[0]);
        Eigen::VectorXd V3dSqrt = V.array().sqrt();
        Eigen::VectorXd V3dInv = V3dSqrt.array().cwiseInverse();
        auto VMsqr = V3dSqrt.asDiagonal();
        auto VMsqrInv = V3dInv.asDiagonal();
        Eigen::MatrixXd SnapMatrix2 = VMsqr * SnapMatrix;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(SnapMatrix2,
                                              Eigen::ComputeThinU | Eigen::ComputeThinV);
        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        eigenValueseig = svd.singularValues().real();
        eigenVectoreig = svd.matrixU().real();
        Eigen::MatrixXd modesEig = VMsqrInv * eigenVectoreig;
        GeometricField<Type, PatchField, GeoMesh> tmb_bu(snapshots[0].name(),
                snapshots[0] * 0);

        for (label i = 0; i < nmodes; i++)
        {
            Eigen::VectorXd vec = modesEig.col(i);
            tmb_bu = Foam2Eigen::Eigen2field(tmb_bu, vec, correctBC);
            modes.set(i, tmb_bu.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        //exportBases(modes, snapshots, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup",
                                      "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModesSVD(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

template void getModesSVD(
    PtrList<volVectorField>& snapshots, PtrList<volVectorField>& modes,
    word fieldName, bool podex, bool supex, bool sup, label nmodes,
    bool correctBC);

/// Construct the Correlation Matrix for Scalar Field
template<>
Eigen::MatrixXd corMatrix(PtrList<volScalarField>& snapshots)
{
    Info << "########## Filling the correlation matrix for " << snapshots[0].name()
         << "##########" << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = fvc::domainIntegrate(snapshots[i] * snapshots[j]).value();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}


/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(PtrList<volVectorField>& snapshots)
{
    Info << "########## Filling the correlation matrix for " << snapshots[0].name()
         << "##########" << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = fvc::domainIntegrate(snapshots[i] & snapshots[j]).value();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}

/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(List<Eigen::SparseMatrix<double>>&
                          snapshots)
{
    Info << "########## Filling the correlation matrix for the matrix list ##########"
         << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            double res = 0;

            for (label k = 0; k < snapshots[i].cols(); k++)
            {
                res += snapshots[i].col(k).dot(snapshots[j].col(k));
            }

            matrix(i, j) = res;
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}

/// Construct the Correlation Matrix for Vector Field
template<>
Eigen::MatrixXd corMatrix(List<Eigen::VectorXd>& snapshots)
{
    Info << "########## Filling the correlation matrix for the matrix list ##########"
         << endl;
    Eigen::MatrixXd matrix( snapshots.size(), snapshots.size());

    for (label i = 0; i < snapshots.size(); i++)
    {
        for (label j = 0; j <= i; j++)
        {
            matrix(i, j) = (snapshots[i].transpose() * snapshots[j]).trace();
        }
    }

    for (label i = 1; i < snapshots.size(); i++)
    {
        for (label j = 0; j < i; j++)
        {
            matrix(j, i) = matrix(i, j);
        }
    }

    return matrix;
}



/// Export the Bases
template<class Type, template<class> class PatchField, class GeoMesh>
void exportBases(PtrList<GeometricField<Type, PatchField, GeoMesh>>& s,
                 PtrList<GeometricField<Type, PatchField, GeoMesh>>& bases,
                 word fieldName, bool sup)
{
    if (sup)
    {
        fileName fieldname;

        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/supremizer/" + name(i + 1));
            fieldname = "./ITHACAoutput/supremizer/" + name(i + 1) + "/" +
                        bases[i].name();
            OFstream os(fieldname);
            bases[i].writeHeader(os);
            os << s[i] << endl;
        }
    }
    else
    {
        fileName fieldname;

        for (label i = 0; i < s.size(); i++)
        {
            mkDir("./ITHACAoutput/POD/" + name(i + 1));
            fieldname = "./ITHACAoutput/POD/" + name(i + 1) + "/" + bases[i].name();
            OFstream os(fieldname);
            bases[i].writeHeader(os);
            os << s[i] << endl;
        }
    }
}
template void exportBases(PtrList<volVectorField>& s,
                          PtrList<volVectorField>& bases, word fieldName, bool sup);
template void exportBases(PtrList<volScalarField>& s,
                          PtrList<volScalarField>& bases, word fieldName, bool sup);

void exportEigenvalues(scalarField Eigenvalues, fileName name,
                       bool sup)
{
    if (sup)
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/supremizer/");
        fieldname = "./ITHACAoutput/supremizer/Eigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < Eigenvalues.size(); i++)
        {
            os << Eigenvalues[i] << endl;
        }
    }
    else
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/POD/");
        fieldname = "./ITHACAoutput/POD/Eigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < Eigenvalues.size(); i++)
        {
            os << Eigenvalues[i] << endl;
        }
    }
}

void exportcumEigenvalues(scalarField cumEigenvalues, fileName name,
                          bool sup)
{
    if (sup)
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/supremizer/");
        fieldname = "./ITHACAoutput/supremizer/cumEigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < cumEigenvalues.size(); i++)
        {
            os << cumEigenvalues[i] << endl;
        }
    }
    else
    {
        fileName fieldname;
        mkDir("./ITHACAoutput/POD/");
        fieldname = "./ITHACAoutput/POD/cumEigenvalues_" + name;
        OFstream os(fieldname);

        for (label i = 0; i < cumEigenvalues.size(); i++)
        {
            os << cumEigenvalues[i] << endl;
        }
    }
}


std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(List<Eigen::SparseMatrix<double>>& A,
                  List<Eigen::VectorXd>& b, label nmodesA, label nmodesB, word MatrixName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<Eigen::SparseMatrix<double>> ModesA(nmodesA);
    List<Eigen::VectorXd> ModesB(nmodesB);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + MatrixName))
    {
        if (nmodesA > A.size() - 2 || nmodesB > A.size() - 2 )
        {
            std::cout <<
                      "The number of requested modes cannot be bigger than the number of Snapshots - 2"
                      << std::endl;
            exit(0);
        }

        scalarField eigenValuesA(nmodesA);
        scalarField cumEigenValuesA(nmodesA);
        scalarField eigenValuesB(nmodesB);
        scalarField cumEigenValuesB(nmodesB);
        List<scalarField> eigenVectorA(nmodesA);
        List<scalarField> eigenVectorB(nmodesB);

        for (label i = 0; i < nmodesA; i++)
        {
            eigenVectorA[i].setSize(A.size());
        }

        for (label i = 0; i < nmodesB; i++)
        {
            eigenVectorB[i].setSize(A.size());
        }

        Eigen::MatrixXd corMatrixA = corMatrix(A);
        Eigen::MatrixXd corMatrixB = corMatrix(b);
        Info << "####### Performing the POD for the Matrix List #######" << endl;
        Spectra::DenseSymMatProd<double> opA(corMatrixA);
        Spectra::DenseSymMatProd<double> opB(corMatrixB);
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                esA(&opA, nmodesA, nmodesA + 10);
        Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                esB(&opB, nmodesB, nmodesB + 10);
        esA.init();
        esB.init();
        esA.compute();

        if (nmodesB != 1)
        {
            esB.compute();
        }

        Info << "####### End of the POD for the Matrix List #######" << endl;
        Eigen::VectorXd eigenValueseigA;
        Eigen::MatrixXd eigenVectorseigA;
        Eigen::VectorXd eigenValueseigB;
        Eigen::MatrixXd eigenVectorseigB;

        if (esA.info() == Spectra::SUCCESSFUL)
        {
            eigenValueseigA = esA.eigenvalues().real();
            eigenVectorseigA = esA.eigenvectors().real();

            if (esB.info() == Spectra::SUCCESSFUL && nmodesB != 1)
            {
                eigenValueseigB = esB.eigenvalues().real();
                eigenVectorseigB = esB.eigenvectors().real();
            }
            else if (nmodesB == 1)
            {
                eigenValueseigB.resize(1);
                eigenVectorseigB.resize(A.size(), nmodesB);
                eigenValueseigB(0) = 1;
                eigenVectorseigB = eigenVectorseigB * 0;
                eigenVectorseigB(0, 0) = 1;
            }
            else
            {
                Info << "The Eigenvalue solver in ITHACAPOD.H did not converge, exiting the code"
                     << endl;
                exit(0);
            }
        }
        else
        {
            Info << "The Eigenvalue solver in ITHACAPOD.H did not converge, exiting the code"
                 << endl;
            exit(0);
        }

        for (label i = 0; i < nmodesA; i++)
        {
            eigenValuesA[i] = eigenValueseigA(i) / eigenValueseigA.sum();
        }

        for (label i = 0; i < nmodesB; i++)
        {
            eigenValuesB[i] = eigenValueseigB(i) / eigenValueseigB.sum();
        }

        for (label i = 0; i < nmodesA; i++)
        {
            for (label k = 0; k < A.size(); k++)
            {
                eigenVectorA[i][k] = eigenVectorseigA(k, i);
            }
        }

        for (label i = 0; i < nmodesB; i++)
        {
            for (label k = 0; k < A.size(); k++)
            {
                eigenVectorB[i][k] = eigenVectorseigB(k, i);
            }
        }

        cumEigenValuesA[0] = eigenValuesA[0];
        cumEigenValuesB[0] = eigenValuesB[0];

        for (label i = 1; i < nmodesA; i++)
        {
            cumEigenValuesA[i] = cumEigenValuesA[i - 1] + eigenValuesA[i];
        }

        for (label i = 1; i < nmodesB; i++)
        {
            cumEigenValuesB[i] = cumEigenValuesB[i - 1] + eigenValuesB[i];
        }

        Eigen::SparseMatrix<double> tmp_A;
        Eigen::VectorXd tmp_B;

        for (label i = 0; i < nmodesA; i++)
        {
            tmp_A = eigenVectorA[i][0] * A[0];

            for (label k = 1; k < A.size(); k++)
            {
                tmp_A += eigenVectorA[i][k] * A[k];
            }

            ModesA[i] = tmp_A;
        }

        for (label i = 0; i < nmodesB; i++)
        {
            tmp_B = eigenVectorB[i][0] * b[0];

            for (label k = 1; k < A.size(); k++)
            {
                tmp_B += eigenVectorB[i][k] * b[k];
            }

            ModesB[i] = tmp_B;
        }

        ITHACAstream::exportList(eigenValuesA,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "eigenValuesA_" + MatrixName);
        ITHACAstream::exportList(cumEigenValuesA,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "cumEigenValuesA_" + MatrixName);
        ITHACAstream::exportList(eigenValuesB,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "eigenValuesB_" + MatrixName);
        ITHACAstream::exportList(cumEigenValuesB,
                                 "./ITHACAoutput/DEIM/" + MatrixName + "/", "cumEigenValuesB_" + MatrixName);

        for (label i = 0; i < ModesA.size(); i++)
        {
            ITHACAstream::SaveSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < ModesB.size(); i++)
        {
            ITHACAstream::SaveDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }
    else
    {
        for (label i = 0; i < nmodesA; i++)
        {
            ITHACAstream::ReadSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < nmodesB; i++)
        {
            ITHACAstream::ReadDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(ModesA, ModesB);
    return tupla;
}

void GrammSchmidt(Eigen::MatrixXd& Matrix)
{
    Eigen::MatrixXd Ortho = Matrix;
    Ortho = Matrix;

    for (label i = 0; i <  Matrix.cols(); i++)
    {
        for (label k = 0; k < i; k++)
        {
            double num = Ortho.col(k).transpose() * Matrix.col(i);
            double den = (Ortho.col(k).transpose() * Ortho.col(k));
            double fact = num / den;
            Ortho.col(i) -= fact * Ortho.col(k) ;
        }

        Ortho.col(i).normalize();
    }

    Matrix = Ortho;
}
template<class Type, template<class> class PatchField, class GeoMesh>
void getModes(
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    PtrList<volScalarField>& Volumes, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if (nmodes == 0 && para->eigensolver == "spectra")
    {
        nmodes = snapshots.size() - 2;
    }

    if (nmodes == 0 && para->eigensolver == "eigen")
    {
        nmodes = snapshots.size();
    }

    if (para->eigensolver == "spectra")
    {
        M_Assert(nmodes <= snapshots.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
    }

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd V = Foam2Eigen::PtrList2Eigen(Volumes);
        Eigen::MatrixXd _corMatrix(snapshots.size(), snapshots.size());
        Info << "Filling the correlation matrix for field " << snapshots[0].name() <<
             endl;

        for (label i = 0; i < snapshots.size(); i++)
        {
            for (label j = 0; j <= i; j++)
            {
                Eigen::VectorXd Mij = (V.col(i).array() * V.col(j).array());
                Mij = Mij.array().abs().sqrt();
                _corMatrix(i, j) = SnapMatrix.col(i).transpose() * Mij.asDiagonal() *
                                   SnapMatrix.col(j);
            }
        }

        std::cout << std::endl;

        for (label i = 1; i < snapshots.size(); i++)
        {
            for (label j = 0; j < i; j++)
            {
                _corMatrix(j, i) = _corMatrix(i, j);
            }
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition for " <<
             snapshots[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig) *
                            eigenValueseigLam.asDiagonal();
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh> tmp2(snapshots[0].name(),
                    snapshots[0] * 0);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec, correctBC);

            // Adjusting boundary conditions
            for (label k = 0; k < tmp2.boundaryField().size(); k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() << " #######"
             << endl;
        exportBases(modes, snapshots, fieldName, sup);
        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, "Usup", "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}
template void getModes(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& modes,
    PtrList<volScalarField>& Volumes, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC);

template void getModes(
    PtrList<volVectorField>& snapshots, PtrList<volVectorField>& modes,
    PtrList<volScalarField>& Volumes, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC);

template<typename type_matrix>
std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
        DEIMmodes(PtrList<type_matrix>& MatrixList, label nmodesA, label nmodesB,
                  word MatrixName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    List<Eigen::SparseMatrix<double>> ModesA(nmodesA);
    List<Eigen::VectorXd> ModesB(nmodesB);

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + MatrixName))
    {
        M_Assert(nmodesA <= MatrixList.size() - 2
                 && nmodesB <= MatrixList.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> snapshots =
                    Foam2Eigen::LFvMatrix2LSM(MatrixList);
        Eigen::MatrixXd corMatrixA = corMatrix(std::get<0>(snapshots));
        Eigen::MatrixXd corMatrixB = corMatrix(std::get<1>(snapshots));
        Eigen::VectorXd eigenValueseigA;
        Eigen::MatrixXd eigenVectorseigA;
        Eigen::VectorXd eigenValueseigB;
        Eigen::MatrixXd eigenVectorseigB;

        if (para->eigensolver == "spectra")
        {
            Info << "####### Performing the POD decomposition for the Matrix List using Spectra #######"
                 << endl;
            Spectra::DenseSymMatProd<double> opA(corMatrixA);
            Spectra::DenseSymMatProd<double> opB(corMatrixB);
            label ncvA = MatrixList.size();
            label ncvB = MatrixList.size();
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                    esA(&opA, nmodesA, ncvA);
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra:: DenseSymMatProd<double>>
                    esB(&opB, nmodesB, ncvB);
            esA.init();
            esB.init();
            esA.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(esA.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");

            if (nmodesB != 1)
            {
                esB.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
                M_Assert(esB.info() == Spectra::SUCCESSFUL,
                         "The Eigenvalue Decomposition did not succeed");
            }

            Info << "####### End of the POD decomposition for the Matrix List #######" <<
                 endl;
            eigenValueseigA = esA.eigenvalues().real();
            eigenVectorseigA = esA.eigenvectors().real();

            if (nmodesB != 1)
            {
                eigenValueseigB = esB.eigenvalues().real();
                eigenVectorseigB = esB.eigenvectors().real();
            }
            else
            {
                eigenValueseigB.resize(1);
                eigenVectorseigB.resize(MatrixList.size(), nmodesB);
                eigenValueseigB(0) = 1;
                eigenVectorseigB = eigenVectorseigB * 0;
                eigenVectorseigB(0, 0) = 1;
            }
        }
        else if (para->eigensolver == "eigen")
        {
            Info << "####### Performing the POD decomposition for the Matrix List using Eigen #######"
                 << endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEgA;
            esEgA.compute(corMatrixA);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEgB;
            esEgB.compute(corMatrixB);
            M_Assert(esEgA.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            M_Assert(esEgB.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectorseigA = esEgA.eigenvectors().real().rowwise().reverse().leftCols(
                                   nmodesA);
            eigenValueseigA = esEgA.eigenvalues().real().reverse().head(nmodesA);
            eigenVectorseigB = esEgB.eigenvectors().real().rowwise().reverse().leftCols(
                                   nmodesB);
            eigenValueseigB = esEgB.eigenvalues().real().reverse().head(nmodesB);
        }

        Eigen::SparseMatrix<double> tmp_A;
        Eigen::VectorXd tmp_B;

        for (label i = 0; i < nmodesA; i++)
        {
            tmp_A = eigenVectorseigA(0, i) * std::get<0>(snapshots)[0];

            for (label k = 1; k < MatrixList.size(); k++)
            {
                tmp_A += eigenVectorseigA(k, i) * std::get<0>(snapshots)[k];
            }

            ModesA[i] = tmp_A;
        }

        for (label i = 0; i < nmodesB; i++)
        {
            tmp_B = eigenVectorseigB(0, i) * std::get<1>(snapshots)[0];

            for (label k = 1; k < MatrixList.size(); k++)
            {
                tmp_B += eigenVectorseigB(k, i) * std::get<1>(snapshots)[k];
            }

            ModesB[i] = tmp_B;
        }

        eigenValueseigA = eigenValueseigA / eigenValueseigA.sum();
        eigenValueseigB = eigenValueseigB / eigenValueseigB.sum();
        Eigen::VectorXd cumEigenValuesA(eigenValueseigA);
        Eigen::VectorXd cumEigenValuesB(eigenValueseigB);

        for (label j = 1; j < cumEigenValuesA.size(); ++j)
        {
            cumEigenValuesA(j) += cumEigenValuesA(j - 1);
        }

        for (label j = 1; j < cumEigenValuesB.size(); ++j)
        {
            cumEigenValuesB(j) += cumEigenValuesB(j - 1);
        }

        for (label i = 0; i < ModesA.size(); i++)
        {
            ITHACAstream::SaveSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < ModesB.size(); i++)
        {
            ITHACAstream::SaveDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }

        Eigen::saveMarketVector(eigenValueseigA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesA", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(eigenValueseigB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/eigenValuesB", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValuesA,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesA", para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValuesB,
                                "./ITHACAoutput/DEIM/" + MatrixName + "/cumEigenValuesB", para->precision,
                                para->outytpe);
    }
    else
    {
        for (label i = 0; i < nmodesA; i++)
        {
            ITHACAstream::ReadSparseMatrix(ModesA[i],
                                           "./ITHACAoutput/DEIM/" + MatrixName + "/", "A_" + MatrixName + name(i));
        }

        for (label i = 0; i < nmodesB; i++)
        {
            ITHACAstream::ReadDenseMatrix(ModesB[i],
                                          "./ITHACAoutput/DEIM/" + MatrixName + "/", "B_" + MatrixName + name(i));
        }
    }

    std::tuple <List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>> tupla;
    tupla = std::make_tuple(ModesA, ModesB);
    return tupla;
}

template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
DEIMmodes(PtrList<fvScalarMatrix>& MatrixList, label nmodesA,
          label nmodesB,
          word MatrixName);

template std::tuple<List<Eigen::SparseMatrix<double>>, List<Eigen::VectorXd>>
DEIMmodes(PtrList<fvVectorMatrix>& MatrixList, label nmodesA,
          label nmodesB,
          word MatrixName);

template<class Type, template<class> class PatchField, class GeoMesh>
PtrList<GeometricField<Type, PatchField, GeoMesh>>DEIMmodes(
            PtrList<GeometricField<Type, PatchField, GeoMesh>>& snapshots, label nmodes,
            word FunctionName, word fieldName)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word PODkey = "POD_" + fieldName;
    word PODnorm = para->ITHACAdict->lookupOrDefault<word>(PODkey, "L2");
    M_Assert(PODnorm == "L2" ||
             PODnorm == "Frobenius", "The PODnorm can be only L2 or Frobenius");
    Info << "Performing POD for " << fieldName << " using the " << PODnorm <<
         " norm" << endl;
    PtrList<GeometricField<Type, fvPatchField, volMesh>> modes;
    bool correctBC = true;

    if (nmodes == 0 && para->eigensolver == "spectra")
    {
        nmodes = snapshots.size() - 2;
    }

    if (nmodes == 0 && para->eigensolver == "eigen")
    {
        nmodes = snapshots.size();
    }

    if (para->eigensolver == "spectra")
    {
        M_Assert(nmodes <= snapshots.size() - 2,
                 "The number of requested modes cannot be bigger than the number of Snapshots - 2");
    }

    if (!ITHACAutilities::check_folder("./ITHACAoutput/DEIM/" + FunctionName))
    {
        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(snapshots);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(snapshots);
        label NBC = snapshots[0].boundaryField().size();
        Eigen::MatrixXd _corMatrix;

        if (PODnorm == "L2")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots);
        }
        else if (PODnorm == "Frobenius")
        {
            _corMatrix = ITHACAutilities::getMassMatrix(snapshots, 0, false);
        }

        if (Pstream::parRun())
        {
            reduce(_corMatrix, sumOp<Eigen::MatrixXd>());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             fieldName << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().array().reverse();
        }

        if (eigenValueseig.array().minCoeff() < 0)
        {
            eigenValueseig = eigenValueseig.array() + 2 * abs(
                                 eigenValueseig.array().minCoeff());
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig);
        Eigen::VectorXd V = ITHACAutilities::getMassMatrixFV(snapshots[0]);
        Eigen::MatrixXd normFact(nmodes, 1);

        for (label i = 0; i < nmodes; i++)
        {
            if (PODnorm == "L2")
            {
                normFact(i, 0) = std::sqrt((modesEig.col(i).transpose() * V.asDiagonal() *
                                            modesEig.col(i))(0, 0));

                if (Pstream::parRun())
                {
                    normFact(i, 0) = (modesEig.col(i).transpose() * V.asDiagonal() *
                                      modesEig.col(i))(0, 0);
                }
            }
            else if (PODnorm == "Frobenius")
            {
                normFact(i, 0) = std::sqrt((modesEig.col(i).transpose() * modesEig.col(i))(0,
                                           0));

                if (Pstream::parRun())
                {
                    normFact(i, 0) = (modesEig.col(i).transpose() * modesEig.col(i))(0, 0);
                }
            }
        }

        if (Pstream::parRun())
        {
            reduce(normFact, sumOp<Eigen::MatrixXd>());
        }

        if (Pstream::parRun())
        {
            normFact = normFact.cwiseSqrt();
        }

        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < NBC; i++)
        {
            modesEigBC[i] = (SnapMatrixBC[i] * eigenVectoreig);
        }

        std::cout << normFact << std::endl;

        for (label i = 0; i < nmodes; i++)
        {
            modesEig.col(i) = modesEig.col(i).array() / normFact(i, 0);

            for (label j = 0; j < NBC; j++)
            {
                modesEigBC[j].col(i) = modesEigBC[j].col(i).array() / normFact(i, 0);
            }
        }

        for (label i = 0; i < modes.size(); i++)
        {
            GeometricField<Type, PatchField, GeoMesh>  tmp2(snapshots[0].name(),
                    snapshots[0]);
            Eigen::VectorXd vec = modesEig.col(i);
            tmp2 = Foam2Eigen::Eigen2field(tmp2, vec, correctBC);

            for (label k = 0; k < NBC; k++)
            {
                ITHACAutilities::assignBC(tmp2, k, modesEigBC[k].col(i));
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;
        ITHACAutilities::createSymLink("./ITHACAoutput/DEIM");

        for (label i = 0; i < modes.size(); i++)
        {
            ITHACAstream::exportSolution(modes[i], name(i + 1), "./ITHACAoutput/DEIM",
                                         fieldName);
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/DEIM/eigenValues_" + fieldName, para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/DEIM/cumEigenValues_" + fieldName, para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;
        ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/DEIM/");
    }

    return modes;
}

template<class Field_type, class Field_type_2>
void getModes(
    PtrList<Field_type>& snapshots, PtrList<Field_type>& modes,
    PtrList<Field_type_2>& fields2, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC)
{
    ITHACAparameters* para(ITHACAparameters::getInstance());

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        if (para->eigensolver == "spectra" )
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size() - 2;
            }

            M_Assert(nmodes <= snapshots.size() - 2,
                     "The number of requested modes cannot be bigger than the number of Snapshots - 2");
        }
        else
        {
            if (nmodes == 0)
            {
                nmodes = snapshots.size();
            }

            M_Assert(nmodes <= snapshots.size(),
                     "The number of requested modes cannot be bigger than the number of Snapshots");
        }

        Eigen::MatrixXd SnapMatrix = Foam2Eigen::PtrList2Eigen(fields2);
        List<Eigen::MatrixXd> SnapMatrixBC = Foam2Eigen::PtrList2EigenBC(fields2);
        label NBC = fields2[0].boundaryField().size();
        auto VM = ITHACAutilities::getMassMatrixFV(fields2[0]);
        Eigen::MatrixXd _corMatrix = SnapMatrix.transpose() * VM.asDiagonal() *
                                     SnapMatrix;

        if (Pstream::parRun())
        {
            reduce(_corMatrix, sumOp<Eigen::MatrixXd>());
        }

        Eigen::VectorXd eigenValueseig;
        Eigen::MatrixXd eigenVectoreig;
        modes.resize(nmodes);
        Info << "####### Performing the POD using EigenDecomposition " <<
             fields2[0].name() << " #######" << endl;
        label ncv = snapshots.size();
        Spectra::DenseSymMatProd<double> op(_corMatrix);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;

        if (para->eigensolver == "spectra")
        {
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>>
                    es(&op, nmodes, ncv);
            std::cout << "Using Spectra EigenSolver " << std::endl;
            es.init();
            es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            M_Assert(es.info() == Spectra::SUCCESSFUL,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = es.eigenvectors().real();
            eigenValueseig = es.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            std::cout << "Using Eigen EigenSolver " << std::endl;
            esEg.compute(_corMatrix);
            M_Assert(esEg.info() == Eigen::Success,
                     "The Eigenvalue Decomposition did not succeed");
            eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(
                                 nmodes);
            eigenValueseig = esEg.eigenvalues().real().reverse().head(nmodes);
        }

        Info << "####### End of the POD for " << snapshots[0].name() << " #######" <<
             endl;
        Eigen::VectorXd eigenValueseigLam =
            eigenValueseig.real().array().cwiseInverse().abs().sqrt() ;
        Eigen::MatrixXd modesEig = (SnapMatrix * eigenVectoreig) *
                                   eigenValueseigLam.asDiagonal();
        List<Eigen::MatrixXd> modesEigBC;
        modesEigBC.resize(NBC);

        for (label i = 0; i < modes.size(); i++)
        {
            Field_type tmp2(snapshots[0].name(), snapshots[0] * 0);

            for (label j = 0; j < snapshots.size(); j++)
            {
                tmp2 += eigenValueseigLam(i) * eigenVectoreig.col(i)(j) * snapshots[j];
            }

            modes.set(i, tmp2.clone());
        }

        eigenValueseig = eigenValueseig / eigenValueseig.sum();
        Eigen::VectorXd cumEigenValues(eigenValueseig);

        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j - 1);
        }

        Info << "####### Saving the POD bases for " << snapshots[0].name() <<
             " #######" << endl;

        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/",
                                       snapshots[0].name());
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", snapshots[0].name());
        }

        Eigen::saveMarketVector(eigenValueseig,
                                "./ITHACAoutput/POD/Eigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
                                "./ITHACAoutput/POD/CumEigenvalues_" + snapshots[0].name(), para->precision,
                                para->outytpe);
    }
    else
    {
        Info << "Reading the existing modes" << endl;

        if (sup == 1)
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}

template void getModes(
    PtrList<surfaceScalarField>& snapshots, PtrList<surfaceScalarField>& modes,
    PtrList<volVectorField>& fields2, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC);

template void getModes(
    PtrList<volScalarField>& snapshots, PtrList<volScalarField>& modes,
    PtrList<volVectorField>& fields2, word fieldName, bool podex, bool supex,
    bool sup, label nmodes, bool correctBC);

template PtrList<volScalarField>
DEIMmodes(
    PtrList<volScalarField>& SnapShotsMatrix,
    label nmodes,
    word FunctionName, word FieldName);

template PtrList<volVectorField>
DEIMmodes(
    PtrList<volVectorField>& SnapShotsMatrix,
    label nmodes,
    word FunctionName, word FieldName);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// CCR Code moved from RedLUM
/* => transfered to IthacaPODTemplate.H 
template<typename T>
ITHACAPODTemplate<T>::ITHACAPODTemplate(fvMesh& mesh, Time& localTime, const word& myfield_name)
{
    ithacaFVParameters = ITHACAparameters::getInstance(mesh, localTime);
    ithacaPODParams    = IthacaPODParameters::getInstance();

    field_name(myfield_name);
    casenameData(ithacaPODParams->get_casenameData());
    l_nSnapshot = ithacaPODParams->get_nSnapshots();//l_nSnapshot(ithacaPODParams->get_nSnapshots());
    l_nBlocks = ithacaPODParams->get_nBlocks();//l_nBlocks(ithacaPODParams->get_nBlocks());
    l_nmodes = ithacaPODParams->get_nModes()[field_name];//l_nmodes(ithacaPODParams->get_nModes()[field_name]);
    l_hilbertSp(ithacaPODParams->get_hilbertSpacePOD()[field_name]);
    weightH1 = ithacaPODParams->get_weightH1();//weightH1(ithacaPODParams->get_weightH1());
    weightBC = ithacaPODParams->get_weightPOD();//weightBC(ithacaPODParams->get_weightPOD());
    patchBC(ithacaPODParams->get_patchBC());
    l_startTime = ithacaPODParams->get_startTime();//l_startTime(ithacaPODParams->get_startTime());
    l_endTime = ithacaPODParams->get_endTime();//l_endTime(ithacaPODParams->get_endTime());
    l_nSnapshotSimulation = ithacaPODParams->get_nSnapshotsSimulation();//l_nSnapshotSimulation(ithacaPODParams->get_nSnapshotsSimulation());
    l_endTimeSimulation = ithacaPODParams->get_endTimeSimulation();//l_endTimeSimulation(ithacaPODParams->get_endTimeSimulation());
    b_centeredOrNot = ithacaPODParams->get_centeredOrNot();//b_centeredOrNot(ithacaPODParams->get_centeredOrNot());
    lambda(Eigen::VectorXd::Zero(l_nmodes));
    w_eigensolver(ithacaPODParams->get_eigensolver());
    i_precision = ithacaPODParams->get_precision();                         //i_precision(ithacaPODParams->get_precision());
    ios_outytpe = ithacaPODParams->get_outytpe();                           //ios_outytpe(ithacaPODParams->get_outytpe());
    #pragma message "CCR - TODO - s'occuper de l'attribut runTime2 ci-dessous"
    //runTime2(Foam::Time::controlDictName, ".", ithacaPODParams->get_casenameData());
}


template<typename T>
ITHACAPODTemplate<T>::~ITHACAPODTemplate()
{
  delete f_field;
  delete f_meanField;
}
*/



void computeLift(PtrList<volTensorField>& Lfield, PtrList<volTensorField>& liftfield, PtrList<volTensorField>& omfield, Eigen::MatrixXi inletIndex)
{
    scalar u_bc;
    scalar u_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        label l = inletIndex(k, 1);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        u_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]).component(l) / area;
        M_Assert(std::abs(u_lf) > 1e-5,
                 "The lift cannot be computed. Please, check your inletIndex definition");

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                u_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]).component(l) / area;
                volTensorField C(Lfield[0].name(), Lfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.append(C.clone());
            }
            else
            {
                u_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]).component(l) / area;
                volTensorField C(Lfield[0].name(), omfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

void computeLift(PtrList<volVectorField>& Lfield, PtrList<volVectorField>& liftfield, PtrList<volVectorField>& omfield, Eigen::MatrixXi inletIndex)
{
    scalar u_bc;
    scalar u_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        label l = inletIndex(k, 1);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        u_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]).component(l) / area;
        M_Assert(std::abs(u_lf) > 1e-5,
                 "The lift cannot be computed. Please, check your inletIndex definition");

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                u_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C(Lfield[0].name(), Lfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.append(C.clone());
            }
            else
            {
                u_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]).component(l) / area;
                volVectorField C(Lfield[0].name(), omfield[j] - liftfield[k]*u_bc / u_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

void computeLift(PtrList<volScalarField>& Lfield, PtrList<volScalarField>& liftfield, PtrList<volScalarField>& omfield, Eigen::MatrixXi inletIndex)
{
    scalar t_bc;
    scalar t_lf;
    scalar area;

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label p = inletIndex(k, 0);
        area = gSum(Lfield[0].mesh().magSf().boundaryField()[p]);
        t_lf = gSum(liftfield[k].mesh().magSf().boundaryField()[p] *
                    liftfield[k].boundaryField()[p]) / area;

        for (label j = 0; j < Lfield.size(); j++)
        {
            if (k == 0)
            {
                t_bc = gSum(Lfield[j].mesh().magSf().boundaryField()[p] *
                            Lfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), Lfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.append(C.clone());
            }
            else
            {
                t_bc = gSum(omfield[j].mesh().magSf().boundaryField()[p] *
                            omfield[j].boundaryField()[p]) / area;
                volScalarField C(Lfield[0].name(), omfield[j] - liftfield[k]*t_bc / t_lf);
                omfield.set(j, C.clone());
            }
        }
    }
}

/* => transfered to IthacaPODTemplate.H 

template<typename T>
void ITHACAPODTemplate<T>::lift(PtrList<T>& snapshots) // déclarer PtrList<T> liftfield, Eigen::MatrixXi inletIndex et T *f_meanField;
{
    if(!f_meanField)
    {
        f_meanField = new T(snapshots[0].name(), snapshots[0] * 0);
    }

    if (lifting)
    {
      PtrList<T> omfield;
      computeLift(snapshots, liftfield, omfield, inletIndex);
      snapshots = omfield;
    }
    if (b_centeredOrNot)
    {
      for (label k = 0; k < snapshots.size(); k++)
      {
        ITHACAutilities::subtractFields(snapshots[k], *(f_meanField));
      }
    }
}

template<typename T>
void ITHACAPODTemplate<T>::lift(T& snapshot)
{
    PtrList<T> snapshots(1);
    snapshots.set(0, new T(snapshot.name(), snapshot));
    lift(snapshots);
    snapshot = snapshots[0];
}
//specialisation
template void ITHACAPODTemplate<volTensorField>::lift(PtrList<volTensorField>& snapshots);
template void ITHACAPODTemplate<volTensorField>::lift(volTensorField& snapshot);
template void ITHACAPODTemplate<volVectorField>::lift(PtrList<volVectorField>& snapshots);
template void ITHACAPODTemplate<volVectorField>::lift(volVectorField& snapshot);
template void ITHACAPODTemplate<volScalarField>::lift(PtrList<volScalarField>& snapshots);
template void ITHACAPODTemplate<volScalarField>::lift(volScalarField& snapshot);






template<typename T>
void ITHACAPODTemplate<T>::compute_lambda(Eigen::MatrixXd& temporalModes)
{

  for (int p = 0; p < l_nmodes; p++)
  {
    for (int i = 0; i < l_nSnapshot; i++)
    {
        lambda(p) += temporalModes(i,p)*temporalModes(i,p);
    }
    lambda(p) /= l_nSnapshot;
  }
}
//specialisation
template void ITHACAPODTemplate<volScalarField>::compute_lambda(Eigen::MatrixXd& temporalModes);
template void ITHACAPODTemplate<volVectorField>::compute_lambda(Eigen::MatrixXd& temporalModes);
template void ITHACAPODTemplate<volTensorField>::compute_lambda(Eigen::MatrixXd& temporalModes);



template<typename T>
void ITHACAPODTemplate<T>::addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
  PtrList<T>& snapshots1,
  PtrList<T>& snapshots2,
  ITHACAPOD::indexSquare& indSquare)
  {

    Info << "Adding the square block [" << indSquare.index1_start << ":" << indSquare.index1_end -1 << "]x["
    << indSquare.index2_start << ":" << indSquare.index2_end -1 << "] to the covariance matrix" << endl;

    Eigen::MatrixXd covMatrixTemp(ITHACAutilities::dot_product_POD(snapshots1,snapshots2, l_hilbertSp, weightBC, patchBC));

    for (label i = 0; i < indSquare.index1_end - indSquare.index1_start; i++)
    {
      for (label j = 0; j < indSquare.index2_end - indSquare.index2_start; j++)
      {
        covMatrix(i + indSquare.index1_start, j + indSquare.index2_start) =
            covMatrixTemp(i,j);
      }
    }
    Info << endl;
 }


template<typename T>
void ITHACAPODTemplate<T>::addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix, PtrList<T>& snapshots, ITHACAPOD::indexTri& indTri)
{

    Info << "Adding the triangular block [" << indTri.index_start << ":" << indTri.index_end - 1 << "]x["
    << indTri.index_start << ":" << indTri.index_end - 1 << "] to the covariance matrix" << endl;

    Eigen::MatrixXd covMatrixTemp(ITHACAutilities::dot_product_POD(snapshots,snapshots, l_hilbertSp, weightBC, patchBC));

    for (label i = 0; i < indTri.index_end - indTri.index_start; i++)
    {
        for (label j = 0; j <= i; j++)
        {
            covMatrix(i + indTri.index_start, j + indTri.index_start) =
                covMatrixTemp(i,j);
        }
    }
    Info << endl;
}
//specialisation
template void ITHACAPODTemplate<volVectorField>::addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volVectorField>& snapshots,
    ITHACAPOD::indexTri& indTri);
template void ITHACAPODTemplate<volVectorField>::addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volVectorField>& snapshots1,
    PtrList<volVectorField>& snapshots2,
    ITHACAPOD::indexSquare& indSquare);
template void ITHACAPODTemplate<volScalarField>::addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volScalarField>& snapshots,
    ITHACAPOD::indexTri& indTri);
template void ITHACAPODTemplate<volScalarField>::addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volScalarField>& snapshots1,
    PtrList<volScalarField>& snapshots2,
    ITHACAPOD::indexSquare& indSquare);
template void ITHACAPODTemplate<volTensorField>::addCovMatrixTriCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volTensorField>& snapshots,
    ITHACAPOD::indexTri& indTri);
template void ITHACAPODTemplate<volTensorField>::addCovMatrixSquareCoeff(Eigen::MatrixXd& covMatrix,
    PtrList<volTensorField>& snapshots1,
    PtrList<volTensorField>& snapshots2,
    ITHACAPOD::indexSquare& indSquare);
 
            


template<typename T>
word ITHACAPODTemplate<T>::nameTempCovMatrix(int& i, int& j)
{
  char str1[10]; sprintf(str1, "%d", i);
  char str2[10]; sprintf(str2, "%d", j);
  word suffix = "_temp_" + static_cast<word>(str1) + "_" + static_cast<word>(str2);
  word filename = folder_covMatrix + name_covMatrix + suffix + ".npy";
  return filename;
}

template<typename T>
void ITHACAPODTemplate<T>::saveTempCovMatrix(Eigen::MatrixXd& covMatrix, int i, int j)
{
  word filename = this->nameTempCovMatrix(i,j);
  cnpy::save(covMatrix, filename);
}


template<typename T>
void ITHACAPODTemplate<T>::deleteTempCovMatrix(int i, int j)
{
  word filename = this->nameTempCovMatrix(i,j);
  std::ifstream fileStr(filename.c_str());
  if (fileStr) remove(filename.c_str());
}
//specialisation
template void ITHACAPODTemplate<volTensorField>::deleteTempCovMatrix(int i, int j);
 



template<typename T>
void ITHACAPODTemplate<T>::deletePreviousTempCovMatrix_N(int* valI, int* valJ, int i, int j, int N)
{
  // NUM = (sum_{n=1}^{n=i} n) = ( (m+(m+i))+(m+1+(m+i)-1)+...+((m+i)+m))/2 = (2*m*i+i+i*i)/2
  //  => NUM(m=0) = (i+i*i)/2
  //     NUM(m=0)+j+1 = (sum_{n=1}^{n=i} n) +j+1 = (1+i)*i/2 +j
  int NUM = (1+i)*i/2 +j;
  int dNUM = NUM-N;
  *valI=0;
  int S = 0;
  while ( dNUM >= S )
  {
    *valI = *valI+1;
    S = (1+*valI)**valI/2;
  }
  *valI = *valI-1;
  S = (1+*valI)**valI/2;
  *valJ = dNUM-S;
  if ( *valI != i )
  {
    *valJ = *valJ-1;
  }
  if ( *valJ == -1 )
  {
    *valI = *valI-1;
    *valJ = *valJ+1;
    }

  if ((*valI>=0) && (*valJ>=0))
  {
    Info << "Last matrix temp saved is ("<< i <<","<<j<<"); with " << N << " previous matrix hold => ";
    Info << "Deleting temp matrix (" << *valI << "," << *valJ << ")" << endl;
    this->deleteTempCovMatrix(*valI,*valJ);
  }
  *valJ = j;
}





template<typename T>
void ITHACAPODTemplate<T>::findTempFile(Eigen::MatrixXd* covMat, int* index1, int* index2)
{
    word pathTemp = name_covMatrix + "_temp_";
    DIR *dir;
    struct dirent *entry;
    dir = opendir(folder_covMatrix.c_str());
    word extTemp = ".npy";
    Info << "looking for " << pathTemp.c_str() << "*" << extTemp.c_str() << " in " << folder_covMatrix.c_str() << endl;
    // INDEX : index1 and index2 of oldest and most recent *_temp_*.npy file
    // INDEX[O][*] : index min (for example with *_temp_0_1.npy and *_temp_3_2.npy => INDEX[0][0]=0 and INDEX[0][1]=1 )
    // INDEX[1][*] : index max (for example with *_temp_0_1.npy and *_temp_3_2.npy => INDEX[1][0]=3 and INDEX[1][1]=2 )
    int INDEX[2][2], N_INDEX=0;
    for (int i=0;i<2;i++)
    { 
      for (int j=0;j<2;j++)
      {
        INDEX[i][j]=-1;
      } 
    }
    if (dir != NULL)
    {
        while ((entry = readdir(dir)) != NULL)
        {
          char ext_name[4]; strncpy(ext_name, "    ", 4);
          if ( strlen(entry->d_name) >= 4 )
          {
            for (int i=0; i<4; i++)
            {
              ext_name[i] = entry->d_name[strlen(entry->d_name)-(4-i)];
            }
          }
          if ((strncmp(entry->d_name, pathTemp.c_str(), pathTemp.size()) == 0) 
                && (strncmp(ext_name, extTemp.c_str(), 4) == 0))
            {
                int num1, num2;
                sscanf(entry->d_name, "%*[^0-9]%d_%d", &num1, &num2);
                *index1 = num1;
                *index2 = num2;
    
            if ( N_INDEX == 0 )
            {
              INDEX[0][0] = *index1; INDEX[0][1] = *index2;
              INDEX[1][0] = *index1; INDEX[1][1] = *index2;
            }
            else
            {
              // updating INDEX from oldest *_temp_*.npy file ( -> INDEX[0][*])
              if ( *index1 < INDEX[0][0] )
              {
                INDEX[0][0]=*index1; INDEX[0][1]=*index2;
              }
              else
              {
                if (( *index1 == INDEX[0][0] ) && ( *index2 < INDEX[0][1] ))
                {
                  INDEX[0][1]=*index2;
                }
              }
              // updating INDEX from most recent *_temp_*.npy file ( -> INDEX[1][*])
              if ( *index1 > INDEX[1][0] )
              {
                INDEX[1][0]=*index1; INDEX[1][1]=*index2;
              }
              else
              {
                if (( *index1 == INDEX[1][0] ) && ( *index2 > INDEX[1][1] ))
                {
                  INDEX[1][1]=*index2;
                }
              }
            }
            N_INDEX = N_INDEX+1;
            Info << "  " << entry->d_name << " FOUND => covMat  can be updated" << endl;
            }
        }
        closedir(dir);
    }
    // restart from oldest *_temp_*.npy file or most recent *_temp_*.npy file
    if ( N_INDEX > 0 )
    {
      // updating matrix from oldest *_temp_*.npy file ( -> INDEX[0][*])
      *index1 = INDEX[0][0]; *index2 = INDEX[0][1];
      Info << "    -> RESTART from INDEX : " << *index1 << " " << *index2 << " (oldest *_temp_*.npy file)" << endl;
      char str1[10]; sprintf(str1, "%d", *index1);
      char str2[10]; sprintf(str2, "%d", *index2);
      word suffix = static_cast<word>(str1) + "_" + static_cast<word>(str2);
      cnpy::load(*covMat, folder_covMatrix + pathTemp + suffix + extTemp);
      Info << "       with covMat size=(" << covMat->rows() << ", " << covMat->cols() << ")" << endl;
    }
}
//Specialisation
template void ITHACAPODTemplate<volTensorField>::findTempFile(Eigen::MatrixXd* covMat, int* index1, int* index2);

template<typename T>
void ITHACAPODTemplate<T>::define_paths()
{
  word pathCentered("");
  if (b_centeredOrNot && !(field_name=="U"))
  {
    pathCentered = "_centered";
  }
  if (lifting) 
  {
    pathCentered += "Lifted";
  }
  word pathHilbertSpace(get_pathHilbertSpace_fromHS(l_hilbertSp));

  // name and folder of the covariance
  // check if the matrix was already computed
  name_covMatrix = "covMatrix" + f_field->name();
  folder_covMatrix = "./ITHACAoutput/CovMatrices" + pathHilbertSpace
                                                  + pathCentered + "/";
  exist_covMatrix = ITHACAutilities::check_file(folder_covMatrix + name_covMatrix+".npy");

  // name and folder of the eigen decomposition
  name_eigenValues = "Eigenvalues_" + f_field->name();
  name_eigenValuesNormalized = "EigenvaluesNormalized_" + f_field->name();
  name_cumEigenValues = "/CumEigenvalues_" + f_field->name();
  name_eigenVector = "/Eigenvector_" + f_field->name() ;
  folder_eigen = "./ITHACAoutput/EigenValuesandVector"
                    + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
  exist_eigenDecomposition = ITHACAutilities::check_file(folder_eigen + name_eigenValues+".npy")
    && ITHACAutilities::check_file(folder_eigen + name_eigenVector+".npy");

  // folder of the spatial modes and check if every modes were already computed
  folder_spatialModes = "./ITHACAoutput/spatialModes"
                    + pathCentered +"_" + std::to_string(l_nmodes) + "modes/";
  exist_spatialModes = ITHACAutilities::check_file(folder_spatialModes + "1/" + f_field->name());

  // folder of the temporal modes and check if modes were already computed
  folder_temporalModes = "./ITHACAoutput/temporalModes"
                    + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
  exist_temporalModes = ITHACAutilities::check_file(folder_temporalModes + f_field->name() + ".npy");

  // folder of the temporal modes and check if modes were already computed
  folder_temporalModesSimulation = "./ITHACAoutput/temporalModesSimulation"
                    + pathCentered + "_" + std::to_string(l_nmodes) + "modes/";
  exist_temporalModesSimulation = ITHACAutilities::check_file(folder_temporalModesSimulation + f_field->name() + ".npy");

  // folder of mean field and check if mean was already computed
  folder_mean = "./ITHACAoutput/mean/";
  exist_noMean = !ITHACAutilities::check_file(folder_mean + "/" + std::to_string(1) + "/" + f_field->name());
}

template<typename T>
void ITHACAPODTemplate<T>::appendMeanfieldtoSpatialModes(PtrList<T>& spatialModes)
{
  if (b_centeredOrNot)
  {
      spatialModes.set(l_nmodes, f_meanField);
  }
}



template<typename T>
PtrList<T> ITHACAPODTemplate<T>::computeSpatialModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig)
{
  PtrList<T> spatialModes;

  // In case the modes were not computed
  if (!exist_spatialModes)
  {
    Info << "Computing the spatial modes of the " << f_field->name() << " field" << endl;

    spatialModes.resize(l_nmodes);
    for (int k = 0; k < l_nmodes; k++)
    {
      spatialModes.set(k, new T(f_field->name(), *f_field) );
      ITHACAutilities::setToZero(spatialModes[k]);
    }

    for (label j = 0; j < l_nSnapshot; j++)
    {
      T snapshotj = *f_field;
      ithacaPODParams->read_snapshot(snapshotj, l_startTime+j);
      if (ithacaPODParams->get_DEIMInterpolatedField() == "nut" && l_hilbertSp == "dL2"){
        ITHACAutilities::multField(snapshotj, ithacaPODParams->get_deltaWeight()); 
      }

      lift(snapshotj);

      for (label k = 0; k < l_nmodes; k++)
      {
        ITHACAutilities::addFields(spatialModes[k], snapshotj, eigenVectoreig(j,k));
      }
    }

    for (int k = 0; k < l_nmodes; k++)
    {
      ITHACAutilities::multField(spatialModes[k], 1/eigenValueseigLam(k));
    }

    mkDir(folder_spatialModes);
    ITHACAstream::exportFields(spatialModes, folder_spatialModes, f_field->name());
  }
  // in the case the spatial modes was already computed, it reads in the hard disk
  else
  {
    Info << "Reading the spatial modes of the " << f_field->name() << " field" << endl;
    ITHACAstream::read_fields(spatialModes, (*f_field), folder_spatialModes);
  }
  Info << endl;

  return spatialModes;
}


template<typename T>
Eigen::MatrixXd ITHACAPODTemplate<T>::computeTemporalModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig)
{
  Eigen::MatrixXd temporalModes(l_nSnapshot, l_nmodes);

  if (!exist_temporalModes)
  {
    Info << "Computing the temporal modes of the " << f_field->name() << " field" << endl;

    mkDir( folder_temporalModes );

    for(int i = 0; i < l_nmodes; i++)
    {
      temporalModes.col(i) = eigenVectoreig.col(i) * eigenValueseigLam(i);
    }

    cnpy::save(temporalModes, folder_temporalModes + f_field->name()+ ".npy");
  }
  else
  {
    Info << "Reading the temporal modes of the " << f_field->name() << " field" << endl;
    cnpy::load(temporalModes, folder_temporalModes + f_field->name()+ ".npy");
  }
  Info << endl;

  return temporalModes;
}



word get_pathHilbertSpace_fromHS(word hilbertSp)
{
  word pathHilbertSpace = "";
  if (hilbertSp == "L2" || hilbertSp == "dL2")
  {
    pathHilbertSpace = "";
  }
  else if (hilbertSp == "L2wBC")
  {
    pathHilbertSpace = "_L2wBC";
  }
  else if (hilbertSp == "H1")
  {
    pathHilbertSpace = "_H1";
  }
  else if (hilbertSp == "wH1")
  {
    pathHilbertSpace = "_wH1";
  }
  else
  {
    Info << "Error: hilbertSpacePOD type " <<
    hilbertSp
    << " is not valid." << endl;
    Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." << endl;
    abort();
  }
  return pathHilbertSpace;
}




template<typename T>
void ITHACAPODTemplate<T>::diagonalisation(Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig)
{
  if(!exist_eigenDecomposition)
  {
    mkDir(folder_eigen);

    Info << "Performing the eigen decomposition" << endl;
    if (w_eigensolver == "spectra")
    {
      Spectra::DenseSymMatProd<double> op(covMatrix);
      Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>> es(&op, l_nmodes, l_nSnapshot);
      std::cout << "Using Spectra EigenSolver " << std::endl;
      es.init();
      es.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
      M_Assert(es.info() == Spectra::SUCCESSFUL, "The Eigenvalue Decomposition did not succeed");
      eigenVectoreig = es.eigenvectors().real();
      eigenValueseig = es.eigenvalues().real();
    }
    else if (w_eigensolver == "eigen")
    {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esEg;
      std::cout << "Using Eigen EigenSolver " << std::endl;
      esEg.compute(covMatrix);
      M_Assert(esEg.info() == Eigen::Success, "The Eigenvalue Decomposition did not succeed");
      eigenVectoreig = esEg.eigenvectors().real().rowwise().reverse().leftCols(l_nmodes);
      eigenValueseig = esEg.eigenvalues().real().reverse().head(l_nmodes);
    }
    else
    {
    }
    // Compute the norm of each Modes
    eigenValueseigLam = eigenValueseig.real().array().abs().sqrt();

    // save the eigen values
    cnpy::save(eigenValueseig, folder_eigen + name_eigenValues + ".npy");

    // save the eigen vectors
    cnpy::save(eigenVectoreig, folder_eigen + "/Eigenvector_" + f_field->name()+ ".npy");

    // save the norm of each Modes
    cnpy::save(eigenValueseigLam, folder_eigen + "/EigenvectorLambda_" + f_field->name()+ ".npy");

    Eigen::VectorXd eigenValueseigNormalized = eigenValueseig / eigenValueseig.sum();
    Eigen::VectorXd cumEigenValues(eigenValueseigNormalized);

    for (int j = 1; j < cumEigenValues.size(); ++j)
    {
      cumEigenValues(j) += cumEigenValues(j - 1);
    }

    // save eigen values normalized
    cnpy::save(eigenValueseigNormalized, folder_eigen + name_eigenValuesNormalized + ".npy");
    // save the cumulated eigen values
    cnpy::save(cumEigenValues, folder_eigen + name_cumEigenValues + ".npy");

  }
  else
  {
    // in case the eigen decomposition was already performed
    // load eigen values
    cnpy::load(eigenValueseig, folder_eigen + name_eigenValues + ".npy");
    // load eigen vectors
    cnpy::load(eigenVectoreig, folder_eigen + "/Eigenvector_" + f_field->name()+ ".npy");
    // load the norm of each Modes
    cnpy::load(eigenValueseigLam, folder_eigen + "/EigenvectorLambda_" + f_field->name()+ ".npy");
  }
  double resolvedVaryingEnergy = eigenValueseig.sum()/l_nSnapshot;
  if (l_hilbertSp == "L2")
  {
    ithacaPODParams->set_resolvedVaryingEnergy( f_field->name(), resolvedVaryingEnergy);
  }
  else
  {
    ithacaPODParams->set_resolvedVaryingEnergy( f_field->name() + "_" + l_hilbertSp, resolvedVaryingEnergy);
  }
  Info << "% of varying " << l_hilbertSp << " energy captures by the modes for " << field_name 
       << " : " << 100.0*eigenValueseig.sum()/covMatrix.trace() << "%" << endl;
  Info << endl;
}

template<typename T>
Eigen::MatrixXd ITHACAPODTemplate<T>::buildCovMatrix()
{
  // the covariance matrix is initialized to a unreal value
  double CovUnrealValue = 999999.999999;
  Eigen::MatrixXd covMatrix;
  covMatrix.setConstant(l_nSnapshot,l_nSnapshot,CovUnrealValue);

  bool exist_covMatrix_bin = ITHACAutilities::check_file(folder_covMatrix + name_covMatrix);
  
  int valI,valJ;
  //initializing value of etapeI and etapeJ : no temp file found
  int etapeI, etapeJ;
  etapeI = -1; etapeJ = -1;

  // create folder
  mkDir(folder_covMatrix);

  // updating covMatrix from previous step (when a *_temp_* file exist)
  findTempFile(&covMatrix, &etapeI, &etapeJ);

  // in case the cov matrix was not computed
  if(!exist_covMatrix && !exist_covMatrix_bin)
  {

    Info << "Computing the covariance matrix of the " << f_field->name() << " field" << endl;
    // the size is cov matrix l_nSnapshot, it has been subdivised in l_nBlocks
    label q = l_nSnapshot/l_nBlocks;
    label r = l_nSnapshot%l_nBlocks;

    Info << "q = " << q << endl;
    Info << "r = " << r << endl;
    Info << "nSnapshot = " << l_nSnapshot << endl;
    Info << "nBlocks = " << l_nBlocks << endl;
    Info << endl;

    PtrList<T> snapshots;

    // If a temp file is found load it and go to the next step
    if(etapeI != -1)
    {
      if(etapeJ == etapeI-1){
        etapeJ = 0;
        etapeI +=1;
      }else{
        etapeJ +=1;
      }
    }else{
      etapeI = 0;
      etapeJ = 0;
    }

    // Number of previous tempory CovMatrix to keep in directory
    // (parameter used by [deletePreviousTempCovMatrix_N])
    int N_previous_temp_Mat = 3;
    
    for (label i = etapeI; i < l_nBlocks; i++)
    {
      ITHACAstream::read_fields(snapshots, (*f_field), casenameData, l_startTime -2 + i*q, q);
      lift(snapshots);

      ITHACAPOD::indexTri indTri;
      indTri.index_start = i*q;
      indTri.index_end = (i+1)*q;

      addCovMatrixTriCoeff(covMatrix, snapshots, indTri);

      PtrList<T> snapshots2;
      
      int valI=etapeI,valJ=etapeJ;

      for (label j = etapeJ; j < i; j++)
      {
        ITHACAstream::read_fields(snapshots2, (*f_field), casenameData, l_startTime -2 + j*q, q);
        lift(snapshots2);

        ITHACAPOD::indexSquare indSquare;
        indSquare.index1_start = i*q;
        indSquare.index1_end = (i+1)*q;
        indSquare.index2_start = j*q;
        indSquare.index2_end = (j+1)*q;

        addCovMatrixSquareCoeff(covMatrix, snapshots,snapshots2, indSquare);

        // Clear the pointer to snapshots2
        snapshots2.clear();

        // Save cove matrix temp file
        saveTempCovMatrix(covMatrix, i,j);

        // Delete previous covMatrix temp file after saving the current one
        deletePreviousTempCovMatrix_N(&valI,&valJ,i,j,N_previous_temp_Mat);
      }

      // Clear the pointer to snapshots
      snapshots.clear();
      valI = i;
      etapeJ = 0;
    }

    if (r != 0)
    {
      PtrList<T> snapshotsEnd;

      ITHACAstream::read_fields(snapshotsEnd, (*f_field), casenameData, l_startTime -2 + l_nBlocks*q, r);
      lift(snapshotsEnd);

      ITHACAPOD::indexTri indTri;
      indTri.index_start = l_nBlocks*q;
      indTri.index_end = l_nSnapshot;

      addCovMatrixTriCoeff(covMatrix, snapshotsEnd, indTri);

      PtrList<T> snapshotsEnd2;

      for (label j = 0; j < l_nBlocks; j++)
      {
        ITHACAstream::read_fields(snapshotsEnd2, (*f_field), casenameData, l_startTime -2 + j*q, q);
        lift(snapshotsEnd2);

        ITHACAPOD::indexSquare indSquare;
        indSquare.index1_start = l_nBlocks*q;
        indSquare.index1_end =  l_nSnapshot;
        indSquare.index2_start = j*q;
        indSquare.index2_end = (j+1)*q;

        addCovMatrixSquareCoeff(covMatrix, snapshotsEnd, snapshotsEnd2, indSquare);

        snapshotsEnd2.clear();
        
        // Save cove matrix temp file
        saveTempCovMatrix(covMatrix, l_nBlocks,j);
        
        // Delete previous covMatrix temp file after saving the current one
        deletePreviousTempCovMatrix_N(&valI,&valJ,l_nBlocks,j,N_previous_temp_Mat);
      }

      snapshotsEnd.clear();
    }
    // covMatrix is symetric, the lower part is used to build the upper part
    covMatrix = covMatrix.selfadjointView<Eigen::Lower>();
    cnpy::save(covMatrix, folder_covMatrix + name_covMatrix+ ".npy");
    
    // Delete previous covMatrix temp file after saving the current one
    if (r == 0) deletePreviousTempCovMatrix_N(&valI,&valJ,l_nBlocks-1,l_nBlocks-1,1);
    if (r != 0) deletePreviousTempCovMatrix_N(&valI,&valJ,l_nBlocks,l_nBlocks,1);
  }
  // in case the cov matrix was already computed, it reads in the hard disk
  else if (exist_covMatrix)
  {
    Info << "Reading the covariance matrix of the " << f_field->name() << " field" << endl;
    cnpy::load(covMatrix, folder_covMatrix + name_covMatrix+ ".npy");
  }
  // in case the cov matrix was already computed, it reads in the hard disk
  else
  {
    Info << "Reading (binary) the covariance matrix of the " << f_field->name() << " field" << endl;
    ITHACAstream::ReadDenseMatrix(covMatrix, folder_covMatrix, name_covMatrix);
    cnpy::save(covMatrix, folder_covMatrix + name_covMatrix+ ".npy");
  }

  // looking for CovUnrealValue
  int covMatrixOK = 1;
  int NbCovUnrealValue = 0;
  for (int i=0;i<l_nSnapshot;i++)
  {
    for (int j=i; j<l_nSnapshot; j++)
    {
      if (covMatrix(i,j) == CovUnrealValue)
      {
        covMatrixOK = 0;
        NbCovUnrealValue += 1;
      }
    }
  }
  if (covMatrixOK == 0) 
  {
    Info << "\n!!! OUPS !!! Unreal value [" << CovUnrealValue << "] found " << NbCovUnrealValue << " times in triangular up part of " << name_covMatrix << " !!!\n" << endl;
    abort();
  }
  
  double varyingEnergy = covMatrix.trace()/l_nSnapshot;
  if (l_hilbertSp == "L2" || l_hilbertSp == "dL2")
  {
    ithacaPODParams->set_varyingEnergy( f_field->name(), varyingEnergy);
  }
  Info << "Total varying " << l_hilbertSp << " energy for " << field_name << " : " << varyingEnergy << endl;
  Info << endl;

  return covMatrix;
}

template<typename T>
void ITHACAPODTemplate<T>::computeMeanField()
{
  if (lifting)
  {
    Info << "Loading lifting functions" << endl;
    liftfield.resize(inletIndex.rows());
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        T* f_lift = new T
        (
          IOobject
          (
            f_field->name() +"lift" + std::to_string(k),
            runTime2.path() + runTime2.times()[1].name(),
            ithacaPODParams->get_mesh(),
            IOobject::MUST_READ
          ),
          ithacaPODParams->get_mesh()
        );
        liftfield.set(k, f_lift );
    }
  }
  if (b_centeredOrNot)
  {

    if (exist_noMean)
    {
      ITHACAutilities::setToZero(*f_meanField);
      Info << "Computing the mean of " << f_field->name() << " field" << endl;
      T snapshotj = *f_field;
      b_centeredOrNot = false;
      for (label j = 0; j < l_nSnapshot; j++)
      {
        // Read the j-th field
        ithacaPODParams->read_snapshot(snapshotj, l_startTime+j);
        lift(snapshotj);
        // add j-th field to meanfield
        ITHACAutilities::addFields(*f_meanField,snapshotj);
      }
      b_centeredOrNot = true;

      ITHACAutilities::multField(*f_meanField, 1/double(l_nSnapshot));

      PtrList<T> meanExport(1);
      meanExport.set(0, new T(f_meanField->name(), *f_meanField));
      ITHACAstream::exportFields(meanExport, folder_mean, f_field->name());
    }
    else
    {
      PtrList<T> meanRead;
      Info << "Reading the mean of " << f_field->name() << " field" << endl;
      ITHACAstream::read_fields(meanRead, (*f_field), folder_mean);
      *f_meanField = meanRead[0];
    }
    double energyMean = ITHACAutilities::dot_product_L2(*f_meanField,*f_meanField);
    double energyHilbertMean = ITHACAutilities::dot_product_POD(*f_meanField,*f_meanField,l_hilbertSp);
    ithacaPODParams->set_meanEnergy( f_meanField->name(), energyMean);
    ithacaPODParams->set_meanEnergy( f_meanField->name() + "_" + l_hilbertSp, energyHilbertMean);
  }
  else
  {
    ITHACAutilities::setToZero(*f_meanField);
  }

  Info << endl;
}

template<typename T>
Eigen::MatrixXd ITHACAPODTemplate<T>::computeSimulationTemporalModes(PtrList<T>& f_spatialModes)
{
  Eigen::MatrixXd temporalModesSimulation(l_nSnapshotSimulation, l_nmodes);

  if (!exist_temporalModesSimulation)
  {
    Info << "Computing the Simulation temporal modes of the " << f_field->name() << " field" << endl;
    mkDir( folder_temporalModesSimulation );

    label l_startTimeSimulation(l_endTime);
    for (label j = 0; j < l_nSnapshotSimulation; j++)
    {

      T snapshotj = *f_field;
      ithacaPODParams->read_snapshot(snapshotj, l_startTimeSimulation + j);
      if (ithacaPODParams->get_DEIMInterpolatedField() == "nut" && l_hilbertSp == "dL2"){
        ITHACAutilities::multField(snapshotj, ithacaPODParams->get_deltaWeight()); 
      }
      lift(snapshotj);

      for (label i = 0; i < l_nmodes; i++)
      {
        temporalModesSimulation(j,i) = ITHACAutilities::dot_product_POD(snapshotj , f_spatialModes[i], l_hilbertSp, weightBC, patchBC,weightH1);
      }
    }
    cnpy::save(temporalModesSimulation, folder_temporalModesSimulation + f_field->name()+ ".npy");
  }
  else
  {
    Info << "Reading the Simulation temporal modes of the " << f_field->name() << " field" << endl;
    cnpy::load(temporalModesSimulation, folder_temporalModesSimulation + f_field->name()+ ".npy");
  }
  Info << endl;

  return temporalModesSimulation;
}


template<typename T>
void ITHACAPODTemplate<T>::getModes(PtrList<T>& spatialModes,
  Eigen::MatrixXd& temporalModes, Eigen::MatrixXd& temporalModesSimulation,
  Eigen::MatrixXd& covMatrix)
{
  Info << "-------------------------------------------------------------------------------------------" << endl;
  Info << "The POD is performing with the following parameters :" << endl;
  Info << "Field : " << f_field->name() << endl;
  Info << "Number of modes : " << l_nmodes << endl;
  Info << "POD Hilbert space : " << l_hilbertSp << endl;
  Info << "Weight for the boundary conditions (0 for L2 POD Hilbert space) : " << weightBC << endl;
  Info<< "Patch for the boundary conditions(inlet by default) used in L2wBC : "<<patchBC<<endl;

  Info << "start time : " << l_startTime << endl;
  Info << "End time : " << l_endTime << endl;
  Info << "Number of snapshots : " << l_nSnapshot << endl;
  Info << "Number of test snapshots : " << l_nSnapshotSimulation << endl;

  Info << "Number of blocks : " << l_nBlocks << endl;
  Info << "Centered datas or not : " << b_centeredOrNot << " (1 centered, 0 not centered)" << endl;

  Info << "Name of eigensolver used : " << w_eigensolver << endl;

  Info << "Results folder : " << "ITHACAoutput" << endl;

  computeMeanField();
  covMatrix = this->buildCovMatrix();
  // Reduction (parallelization)
  Eigen::VectorXd eigenValueseig = Eigen::VectorXd::Zero(l_nmodes);
  Eigen::MatrixXd eigenVectoreig = Eigen::MatrixXd::Zero(l_nSnapshot,l_nmodes);
  diagonalisation(covMatrix, eigenValueseig, eigenVectoreig);
  spatialModes = computeSpatialModes(eigenValueseig, eigenVectoreig);
  temporalModes = computeTemporalModes(eigenValueseig, eigenVectoreig);
  compute_lambda(temporalModes);
  temporalModesSimulation = computeSimulationTemporalModes(spatialModes);
  // Reduction (parallelization)

  if (field_name == "U")
  {
      ithacaPODParams->set_eigenValues_U(eigenValueseig);
      ithacaPODParams->set_lambda(lambda);
  }

  Info << "-------------------------------------------------------------------------------------------" << endl;
}



  // Specialisation
  template ITHACAPODTemplate<volTensorField>::ITHACAPODTemplate(fvMesh& mesh, Time& localTime, const word& myfield_name);//template ITHACAPODTemplate<volTensorField>::ITHACAPODTemplate(Parameters* myParameters, const word& myfield_name);
  template ITHACAPODTemplate<volTensorField>::~ITHACAPODTemplate();
  template void ITHACAPODTemplate<volTensorField>::define_paths();
  template void ITHACAPODTemplate<volTensorField>::computeMeanField();
  template void ITHACAPODTemplate<volTensorField>::appendMeanfieldtoSpatialModes(PtrList<volTensorField>& spatialModes);
  //template word ITHACAPODTemplate<volTensorField>::nameTempCovMatrix(int i, int j);
  template void ITHACAPODTemplate<volTensorField>::saveTempCovMatrix(Eigen::MatrixXd& covMatrix, int i, int j);
  template void ITHACAPODTemplate<volTensorField>::deletePreviousTempCovMatrix_N(int* valI, int* valJ, int i, int j, int N);
  template Eigen::MatrixXd ITHACAPODTemplate<volTensorField>::buildCovMatrix();
  template void ITHACAPODTemplate<volTensorField>::diagonalisation(Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template PtrList<volTensorField> ITHACAPODTemplate<volTensorField>::computeSpatialModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template Eigen::MatrixXd ITHACAPODTemplate<volTensorField>::computeTemporalModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template void ITHACAPODTemplate<volTensorField>::getModes(PtrList<volTensorField>& spatialModes, Eigen::MatrixXd& temporalModes,
    Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
  template Eigen::MatrixXd ITHACAPODTemplate<volTensorField>::computeSimulationTemporalModes(PtrList<volTensorField>& spatialModes);

  // Specialisation
  template ITHACAPODTemplate<volVectorField>::ITHACAPODTemplate(fvMesh& mesh, Time& localTime, const word& myfield_name);//template ITHACAPODTemplate<volVectorField>::ITHACAPODTemplate(Parameters* myParameters, const word& myfield_name);
  template ITHACAPODTemplate<volVectorField>::~ITHACAPODTemplate();
  template void ITHACAPODTemplate<volVectorField>::define_paths();
  template void ITHACAPODTemplate<volVectorField>::computeMeanField();
  template void ITHACAPODTemplate<volVectorField>::appendMeanfieldtoSpatialModes(PtrList<volVectorField>& spatialModes);
  template Eigen::MatrixXd ITHACAPODTemplate<volVectorField>::buildCovMatrix();
  template void ITHACAPODTemplate<volVectorField>::diagonalisation(Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template PtrList<volVectorField> ITHACAPODTemplate<volVectorField>::computeSpatialModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template Eigen::MatrixXd ITHACAPODTemplate<volVectorField>::computeTemporalModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template void ITHACAPODTemplate<volVectorField>::getModes(PtrList<volVectorField>& spatialModes, Eigen::MatrixXd& temporalModes,
    Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
  template Eigen::MatrixXd ITHACAPODTemplate<volVectorField>::computeSimulationTemporalModes(PtrList<volVectorField>& spatialModes);

  // Specialisation
  template ITHACAPODTemplate<volScalarField>::ITHACAPODTemplate(fvMesh& mesh, Time& localTime, const word& myfield_name);//template ITHACAPODTemplate<volScalarField>::ITHACAPODTemplate(Parameters* myParameters, const word& myfield_name);
  template ITHACAPODTemplate<volScalarField>::~ITHACAPODTemplate();
  template void ITHACAPODTemplate<volScalarField>::define_paths();
  template void ITHACAPODTemplate<volScalarField>::computeMeanField();
  template void ITHACAPODTemplate<volScalarField>::appendMeanfieldtoSpatialModes(PtrList<volScalarField>& spatialModes);
  template Eigen::MatrixXd ITHACAPODTemplate<volScalarField>::buildCovMatrix();
  template void ITHACAPODTemplate<volScalarField>::diagonalisation(Eigen::MatrixXd& covMatrix, Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template PtrList<volScalarField> ITHACAPODTemplate<volScalarField>::computeSpatialModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template Eigen::MatrixXd ITHACAPODTemplate<volScalarField>::computeTemporalModes(Eigen::VectorXd& eigenValueseig, Eigen::MatrixXd& eigenVectoreig);
  template void ITHACAPODTemplate<volScalarField>::getModes(PtrList<volScalarField>& spatialModes, Eigen::MatrixXd& temporalModes,
    Eigen::MatrixXd& temporalModesSimulation, Eigen::MatrixXd& covMatrix);
  template Eigen::MatrixXd ITHACAPODTemplate<volScalarField>::computeSimulationTemporalModes(PtrList<volScalarField>& spatialModes);
*/


// end of CCR Code moved from RedLUM
}
