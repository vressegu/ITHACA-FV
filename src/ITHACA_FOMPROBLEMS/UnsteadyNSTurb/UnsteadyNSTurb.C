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

#include "UnsteadyNSTurb.H"

/// \file
/// Source file of the unsteadyNS class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSTurb::UnsteadyNSTurb() {}

// Construct from zero
UnsteadyNSTurb::UnsteadyNSTurb(int argc, char* argv[])
{
    _args = autoPtr<argList>
            (
                new argList(argc, argv)
            );

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>
              (
                  new pimpleControl
                  (
                      mesh
                  )
              );
#include "createFields.H"
#include "createUfIfPresent.H"
#include "createFvOptions.H"
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    timeDerivativeSchemeOrder =
        ITHACAdict->lookupOrDefault<word>("timeDerivativeSchemeOrder", "second");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
             "The BC method must be set to lift or penalty in ITHACAdict");
    M_Assert(timeDerivativeSchemeOrder == "first"
             || timeDerivativeSchemeOrder == "second",
             "The time derivative approximation must be set to either first or second order scheme in ITHACAdict");
    para = ITHACAparameters::getInstance(mesh, runTime);
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
    /// Number of velocity modes to be calculated
    NUmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesUout", 15);
    /// Number of pressure modes to be calculated
    NPmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesPout", 15);
    /// Number of nut modes to be calculated
    NNutModesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesNutOut", 15);
    /// Number of velocity modes used for the projection
    NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
    /// Number of supremizers modes used for the projection
    NSUPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesSUPproj", 10);
    /// Number of pressure modes used for the projection
    NPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesPproj", 10);
    /// Number of nut modes used for the projection
    NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTurb::truthSolve(List<scalar> mu_now, std::string& offlinepath)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& nut = _nut();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    instantList Times = runTime.times();
    label& pRefCell = _pRefCell;
    scalar& pRefValue = _pRefValue;
    mesh.setFluxRequired(p.name());
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime + writeEvery;
    label nsnapshots = 0;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        ++runTime;
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            nsnapshots += 1;
            // Produces error when uncommented
            // volScalarField nut = turbulence->nut().ref();
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), offlinepath);
            ITHACAstream::exportSolution(p, name(counter), offlinepath);
            ITHACAstream::exportSolution(nut, name(counter), offlinepath);
            std::ofstream of(offlinepath + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            nutFields.append(tmp<volScalarField>(nut));
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
            // --- Fill in the mu_samples with parameters (time, mu) to be used for the PODI sample points
            mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size() + 1);
            mu_samples(mu_samples.rows() - 1, 0) = atof(runTime.timeName().c_str());

            for (label i = 0; i < mu_now.size(); i++)
            {
                mu_samples(mu_samples.rows() - 1, i + 1) = mu_now[i];
            }
        }
    }

    // Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == nsnapshots* mu.cols())
    {
        ITHACAstream::exportToFile(mu_samples, "mu_samples", "eigen",
                                   offlinepath);
    }
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceTensor1(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor;
    ct1Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                         nutModes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct1Tensor,
                                  "ct1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct1Tensor;
}


Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceAveTensor1(label NUmodes,
        label NSUPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1AveTensor;
    label samplesNumber = nutAve.size();
    ct1AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(
                                            nutAve[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct1AveTensor,
                                  "ct1Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return ct1AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPETensor1(label NUmodes,
        label NSUPmodes, label NPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1PPETensor;
    ct1PPETensor.resize(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct1PPETensor(i, j, k) = fvc::domainIntegrate(2 * Pmodes[i] * (fvc::laplacian(
                //                             L_U_SUPmodes[k]) & fvc::grad(nutModes[j]))).value();
                // ct1PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::div(
                //     fvc::laplacian(
                //         nutModes[j], L_U_SUPmodes[k])))).value();
                ct1PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (
                                            fvc::laplacian(
                                                nutModes[j], L_U_SUPmodes[k]))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct1PPETensor, 
                                  "ct1PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct1PPETensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPEAveTensor1(label NUmodes,
        label NSUPmodes, label NPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct1PPEAveTensor;
    label samplesNumber = nutAve.size();
    ct1PPEAveTensor.resize(NPmodes, samplesNumber, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(2 * Pmodes[i] * (fvc::laplacian(
                //                                L_U_SUPmodes[k]) & fvc::grad(nutAve[j]))).value();
                // ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::div(
                //  fvc::laplacian(
                //      nutAve[j], L_U_SUPmodes[k])))).value();
                ct1PPEAveTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & (
                                               fvc::laplacian(
                                                   nutAve[j], L_U_SUPmodes[k]))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct1PPEAveTensor, 
                                  "ct1PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct1PPEAveTensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceTensor2(label NUmodes,
        label NSUPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor;
    ct2Tensor.resize(cSize, nNutModes, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                         nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct2Tensor,
                                  "ct2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(nNutModes) + "_t");
    return ct2Tensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulenceAveTensor2(label NUmodes,
        label NSUPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2AveTensor;
    label samplesNumber = nutAve.size();
    ct2AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(
                                            nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct2AveTensor, 
                                  "ct2Ave_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_t");
    return ct2AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPETensor2(label NUmodes,
        label NSUPmodes, label NPmodes, label nNutModes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2PPETensor;
    ct2PPETensor.resize(NPmodes, nNutModes, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < nNutModes; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct2PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::grad(fvc::grad(
                //                             nutModes[j])) && (dev2((fvc::grad(L_U_SUPmodes[k]))() + fvc::grad(
                //                                         L_U_SUPmodes[k]))().T()))).value();
                // ct2PPETensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * ((fvc::div(fvc::div(
                //     nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))))).value();
                ct2PPETensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & ((fvc::div(
                                            nutModes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct2PPETensor, 
                                  "ct2PPE_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t");
    return ct2PPETensor;
}

Eigen::Tensor<double, 3> UnsteadyNSTurb::turbulencePPEAveTensor2(label NUmodes,
        label NSUPmodes, label NPmodes)
{
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::Tensor<double, 3> ct2PPEAveTensor;
    label samplesNumber = nutAve.size();
    ct2PPEAveTensor.resize(NPmodes, samplesNumber, cSize);

    for (label i = 0; i < NPmodes; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                // ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * (fvc::grad(
                //                                fvc::grad(
                //                                    nutAve[j])) && (dev2((fvc::grad(L_U_SUPmodes[k]))() + fvc::grad(
                //                                            L_U_SUPmodes[k]))().T()))).value();
                // ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(Pmodes[i] * ((fvc::div(fvc::div(
                //  nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T())))))).value();
                ct2PPEAveTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(Pmodes[i]) & ((
                                               fvc::div(
                                                   nutAve[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))).value();
            }
        }
    }

    // Export the tensor
    ITHACAstream::exportToFile(ct2PPEAveTensor, 
                                  "ct2PPEAve_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                                      NSUPmodes) + "_" + name(NPmodes) + "_t");
    return ct2PPEAveTensor;
}

Eigen::MatrixXd UnsteadyNSTurb::btTurbulence(label NUmodes, label NSUPmodes)
{
    label btSize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
    btMatrix = btMatrix * 0;

    // Project everything
    for (label i = 0; i < btSize; i++)
    {
        for (label j = 0; j < btSize; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(dev2((T(
                    fvc::grad(
                        L_U_SUPmodes[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportToFile(btMatrix,
                                  "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    return btMatrix;
}

void UnsteadyNSTurb::projectSUP(fileName folder, label NU, label NP, label NSUP,
                                label Nnut, bool rbfInterp)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    nNutModes = Nnut;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(liftfield[k]));
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(Umodes[k]));
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(supmodes[k]));
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word bStr = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes);

        if (!ITHACAstream::importNpy(B_matrix, bStr))
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (!ITHACAstream::importNpy(btMatrix, btStr))
        {
            btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        word kStr = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes) + "_" + name(NPmodes);

        if (!ITHACAstream::importNpy(K_matrix, kStr))
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word pStr = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes) + "_" + name(NPmodes);

        if (!ITHACAstream::importNpy(P_matrix, pStr))
        {
            P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        }

        word mStr = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                        NSUPmodes);

        if (!ITHACAstream::importNpy(M_matrix, mStr))
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (!ITHACAstream::importNpy(C_tensor, C_str))
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct1Tensor, ct1Str))
        {
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct2Tensor, ct2Str))
        {
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct1AveTensor, ct1AveStr))
            {
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            }

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct2AveTensor, ct2AveStr))
            {
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
            }
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "python");
        ITHACAstream::exportToFile(K_matrix, "K", "python");
        ITHACAstream::exportToFile(P_matrix, "P", "python");
        ITHACAstream::exportToFile(M_matrix, "M", "python");
        ITHACAstream::exportToFile(btMatrix, "bt", "python");
        ITHACAstream::exportToFile(C_tensor, "C", "python");
        ITHACAstream::exportToFile(ct1Tensor, "ct1", "python");
        ITHACAstream::exportToFile(ct2Tensor, "ct2", "python");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave", "python");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave", "python");
        }
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "matlab");
        ITHACAstream::exportToFile(K_matrix, "K", "matlab");
        ITHACAstream::exportToFile(P_matrix, "P", "matlab");
        ITHACAstream::exportToFile(M_matrix, "M", "matlab");
        ITHACAstream::exportToFile(btMatrix, "bt", "matlab");
        ITHACAstream::exportToFile(C_tensor, "C", "matlab");
        ITHACAstream::exportToFile(ct1Tensor, "ct1", "matlab");
        ITHACAstream::exportToFile(ct2Tensor, "ct2", "matlab");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave", "matlab");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave", "matlab");
        }
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "eigen");
        ITHACAstream::exportToFile(K_matrix, "K", "eigen");
        ITHACAstream::exportToFile(P_matrix, "P", "eigen");
        ITHACAstream::exportToFile(M_matrix, "M", "eigen");
        ITHACAstream::exportToFile(btMatrix, "bt", "eigen");
        ITHACAstream::exportToFile(C_tensor, "C", "eigen", "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportToFile(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportToFile(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2Ave");
        }
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    // Define coeffL2
    coeffL2 = ITHACAutilities::getCoeffs(nutFields,
                                         nutModes, nNutModes);
    ITHACAstream::exportToFile(coeffL2, "coeffL2", "python");
    ITHACAstream::exportToFile(coeffL2, "coeffL2", "matlab");

    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;
    }

    if (rbfInterp == true && (!Pstream::parRun()))
    {
        if (ITHACAutilities::check_file("./radii.txt"))
        {
            radii = ITHACAstream::readMatrix("./radii.txt");
            M_Assert(radii.size() == nNutModes,
                     "Thes size of the shape parameters vector must be equal to the number of eddy viscosity modes nNutModes");
        }
        else
        {
            radii = Eigen::MatrixXd::Ones(nNutModes,
                                          1) * e;
        }

        samples.resize(nNutModes);
        rbfSplines.resize(nNutModes);
        Eigen::MatrixXd weights;

        for (label i = 0; i < nNutModes; i++)
        {
            word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                              + name(NUmodes) + "_" + name(NSUPmodes) ;

            if (ITHACAstream::importNpy(weights, weightName, "./ITHACAoutput/weightsSUP/"))
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i],
                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, radii(i));
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
            else
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i],
                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, radii(i));
                ITHACAstream::exportToFile(rbfSplines[i]->weights, weightName, "cnpy",
                                              "./ITHACAoutput/weightsSUP/");
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
        }
    }
}



void UnsteadyNSTurb::projectPPE(fileName folder, label NU, label NP, label NSUP,
                                label Nnut, bool rbfInterp)
{
    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = 0;
    nNutModes = Nnut;
    L_U_SUPmodes.resize(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(liftfield[k]));
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            L_U_SUPmodes.append(tmp<volVectorField>(Umodes[k]));
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (!ITHACAstream::importNpy(B_matrix,B_str))
        {
            B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        }

        word btStr = "bt_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (!ITHACAstream::importNpy(btMatrix,btStr))
        {
            btMatrix = btTurbulence(NUmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes);

        if (!ITHACAstream::importNpy(K_matrix,K_str))
        {
            K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        }

        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes);

        if (!ITHACAstream::importNpy(M_matrix,M_str))
        {
            M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        }

        word D_str = "D_" + name(NPmodes);

        if (!ITHACAstream::importNpy(D_matrix,D_str))
        {
            D_matrix = laplacian_pressure(NPmodes);
        }

        word bc1_str = "BC1_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (!ITHACAstream::importNpy(BC1_matrix,bc1_str))
        {
            BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        }

        word bc2_str = "BC2_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(
                           NSUPmodes) + "_" + name(NPmodes) + "_t";

        if (!ITHACAstream::importNpy(bc2Tensor, bc2_str))
        {
            bc2Tensor = pressureBC2(NUmodes, NPmodes);
        }

        word bc3_str = "BC3_" + name(liftfield.size()) + "_" + name(
                           NUmodes) + "_" + name(NPmodes);

        if (!ITHACAstream::importNpy(BC3_matrix,bc3_str))
        {
            BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_t";

        if (!ITHACAstream::importNpy(C_tensor, C_str))
        {
            C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        }

        word ct1Str = "ct1_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct1Tensor, ct1Str))
        {
            ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        }

        word ct2Str = "ct2_" + name(liftfield.size()) + "_" + name(
                          NUmodes) + "_" + name(
                          NSUPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct2Tensor, ct2Str))
        {
            ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        }

        word G_str = "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(
                         NSUPmodes) + "_" + name(NPmodes) + "_t";

        if (!ITHACAstream::importNpy(gTensor, G_str))
        {
            gTensor = divMomentum(NUmodes, NPmodes);
        }

        if (nutAve.size() != 0)
        {
            word ct1AveStr = "ct1Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct1AveTensor, ct1AveStr))
            {
                ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            }

            word ct2AveStr = "ct2Ave_" + name(liftfield.size()) + "_" + name(
                                 NUmodes) + "_" + name(
                                 NSUPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct2AveTensor, ct2AveStr))
            {
                ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
            }
        }

        word ct1PPEStr = "ct1PPE_" + name(liftfield.size()) + "_" + name(
                             NUmodes) + "_" + name(
                             NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct1PPETensor, ct1PPEStr))
        {
            ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);
        }

        word ct2PPEStr = "ct2PPE_" + name(liftfield.size()) + "_" + name(
                             NUmodes) + "_" + name(
                             NSUPmodes) + "_" + name(NPmodes) + "_" + name(nNutModes) + "_t";

        if (!ITHACAstream::importNpy(ct2PPETensor, ct2PPEStr))
        {
            ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);
        }

        if (nutAve.size() != 0)
        {
            word ct1PPEAveStr = "ct1PPEAve_" + name(liftfield.size()) + "_" + name(
                                    NUmodes) + "_" + name(
                                    NSUPmodes) + "_" + name(NPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct1PPEAveTensor, ct1PPEAveStr))
            {
                ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);
            }

            word ct2PPEAveStr = "ct2PPEAve_" + name(liftfield.size()) + "_" + name(
                                    NUmodes) + "_" + name(
                                    NSUPmodes) + "_" + name(NPmodes) + "_t";

            if (!ITHACAstream::importNpy(ct2PPEAveTensor, ct2PPEAveStr))
            {
                ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
            }
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }
    else
    {
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        D_matrix = laplacian_pressure(NPmodes);
        gTensor = divMomentum(NUmodes, NPmodes);
        BC1_matrix = pressure_BC1(NUmodes, NPmodes);
        bc2Tensor = pressureBC2(NUmodes, NPmodes);
        BC3_matrix = pressure_BC3(NUmodes, NPmodes);
        btMatrix = btTurbulence(NUmodes, NSUPmodes);
        ct1Tensor = turbulenceTensor1(NUmodes, NSUPmodes, nNutModes);
        ct2Tensor = turbulenceTensor2(NUmodes, NSUPmodes, nNutModes);
        ct1PPETensor = turbulencePPETensor1(NUmodes, NSUPmodes, NPmodes, nNutModes);
        ct2PPETensor = turbulencePPETensor2(NUmodes, NSUPmodes, NPmodes, nNutModes);

        if (nutAve.size() != 0)
        {
            ct1AveTensor = turbulenceAveTensor1(NUmodes, NSUPmodes);
            ct2AveTensor = turbulenceAveTensor2(NUmodes, NSUPmodes);
            ct1PPEAveTensor = turbulencePPEAveTensor1(NUmodes, NSUPmodes, NPmodes);
            ct2PPEAveTensor = turbulencePPEAveTensor2(NUmodes, NSUPmodes, NPmodes);
        }

        if (bcMethod == "penalty")
        {
            bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
            bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        }
    }

    // Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "python");
        ITHACAstream::exportToFile(K_matrix, "K", "python");
        ITHACAstream::exportToFile(D_matrix, "D", "python");
        ITHACAstream::exportToFile(M_matrix, "M", "python");
        ITHACAstream::exportToFile(BC1_matrix, "BC1", "python");
        ITHACAstream::exportToFile(BC3_matrix, "BC3", "python");
        ITHACAstream::exportToFile(C_tensor, "C", "python");
        ITHACAstream::exportToFile(gTensor, "G", "python");
        ITHACAstream::exportToFile(bc2Tensor, "BC2", "python");
        ITHACAstream::exportToFile(ct1Tensor, "ct1", "python");
        ITHACAstream::exportToFile(ct2Tensor, "ct2", "python");
        ITHACAstream::exportToFile(ct1PPETensor, "ct1PPE", "python");
        ITHACAstream::exportToFile(ct2PPETensor, "ct2PPE", "python");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave", "python");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave", "python");
            ITHACAstream::exportToFile(ct1PPEAveTensor, "ct1PPEAve", "python");
            ITHACAstream::exportToFile(ct2PPEAveTensor, "ct2PPEAve", "python");
        }
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "matlab");
        ITHACAstream::exportToFile(K_matrix, "K", "matlab");
        ITHACAstream::exportToFile(D_matrix, "D", "matlab");
        ITHACAstream::exportToFile(M_matrix, "M", "matlab");
        ITHACAstream::exportToFile(BC1_matrix, "BC1", "matlab");
        ITHACAstream::exportToFile(BC3_matrix, "BC3", "matlab");
        ITHACAstream::exportToFile(C_tensor, "C", "matlab");
        ITHACAstream::exportToFile(gTensor, "G", "matlab");
        ITHACAstream::exportToFile(bc2Tensor, "BC2", "matlab");
        ITHACAstream::exportToFile(ct1Tensor, "ct1", "matlab");
        ITHACAstream::exportToFile(ct2Tensor, "ct2", "matlab");
        ITHACAstream::exportToFile(ct1PPETensor, "ct1PPE", "matlab");
        ITHACAstream::exportToFile(ct2PPETensor, "ct2PPE", "matlab");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave", "matlab");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave", "matlab");
            ITHACAstream::exportToFile(ct1PPEAveTensor, "ct1PPEAve", "matlab");
            ITHACAstream::exportToFile(ct2PPEAveTensor, "ct2PPEAve", "matlab");
        }
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportToFile(B_matrix, "B", "eigen");
        ITHACAstream::exportToFile(K_matrix, "K", "eigen");
        ITHACAstream::exportToFile(D_matrix, "D", "eigen");
        ITHACAstream::exportToFile(M_matrix, "M", "eigen");
        ITHACAstream::exportToFile(BC1_matrix, "BC1", "eigen");
        ITHACAstream::exportToFile(BC3_matrix, "BC3", "eigen");
        ITHACAstream::exportToFile(C_tensor, "C", "eigen",
                                   "./ITHACAoutput/Matrices/C");
        ITHACAstream::exportToFile(gTensor, "G", "eigen",
                                   "./ITHACAoutput/Matrices/G");
        ITHACAstream::exportToFile(bc2Tensor, "BC2_", "eigen",
                                   "./ITHACAoutput/Matrices/BC2");
        ITHACAstream::exportToFile(btMatrix, "bt", "eigen");
        ITHACAstream::exportToFile(ct1Tensor, "ct1_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1");
        ITHACAstream::exportToFile(ct2Tensor, "ct2_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2");
        ITHACAstream::exportToFile(ct1PPETensor, "ct1PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct1PPE");
        ITHACAstream::exportToFile(ct2PPETensor, "ct2PPE_", "eigen",
                                   "./ITHACAoutput/Matrices/ct2PPE");

        if (nutAve.size() != 0)
        {
            ITHACAstream::exportToFile(ct1AveTensor, "ct1Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1Ave");
            ITHACAstream::exportToFile(ct2AveTensor, "ct2Ave_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2Ave");
            ITHACAstream::exportToFile(ct1PPEAveTensor, "ct1PPEAve_", "eigen",
                                       "./ITHACAoutput/Matrices/ct1PPEAve");
            ITHACAstream::exportToFile(ct2PPEAveTensor, "ct2PPEAve_", "eigen",
                                       "./ITHACAoutput/Matrices/ct2PPEAve");
        }
    }

    bTotalMatrix = B_matrix + btMatrix;
    label cSize = NUmodes + NSUPmodes + liftfield.size();
    cTotalTensor.resize(cSize, nNutModes, cSize);
    cTotalTensor = ct1Tensor + ct2Tensor;
    cTotalPPETensor.resize(NPmodes, nNutModes, cSize);
    cTotalPPETensor = ct1PPETensor + ct2PPETensor;

    if (nutAve.size() != 0)
    {
        cTotalAveTensor.resize(cSize, nutAve.size(), cSize);
        cTotalAveTensor = ct1AveTensor + ct2AveTensor;
        cTotalPPEAveTensor.resize(NPmodes, nutAve.size(), cSize);
        cTotalPPEAveTensor = ct1PPEAveTensor + ct2PPEAveTensor;
    }

    if (rbfInterp == true && (!Pstream::parRun()))
    {
        if (ITHACAutilities::check_file("./radii.txt"))
        {
            radii = ITHACAstream::readMatrix("./radii.txt");
            M_Assert(radii.size() ==  nNutModes,
                     "Thes size of the shape parameters vector must be equal to the number of eddy viscosity modes nNutModes");
        }
        else
        {
            radii = Eigen::MatrixXd::Ones(nNutModes,
                                          1) * e;
        }

        samples.resize(nNutModes);
        rbfSplines.resize(nNutModes);
        Eigen::MatrixXd weights;

        for (label i = 0; i < nNutModes; i++)
        {
            word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_"
                              + name(NUmodes) + "_" + name(NSUPmodes) ;

            if (ITHACAstream::importNpy(weights, weightName, "./ITHACAoutput/weightsPPE/"))
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i],
                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, radii(i));
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
            else
            {
                samples[i] = new SPLINTER::DataTable(1, 1);

                for (label j = 0; j < coeffL2.cols(); j++)
                {
                    samples[i]->addSample(velRBF.row(j), coeffL2(i, j));
                }

                rbfSplines[i] = new SPLINTER::RBFSpline( * samples[i],
                    SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, radii(i));
                ITHACAstream::exportToFile(rbfSplines[i]->weights, weightName, "cnpy"
                                              "./ITHACAoutput/weightsPPE/");
                std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
            }
        }
    }
}

List < Eigen::MatrixXd > UnsteadyNSTurb::velDerivativeCoeff(Eigen::MatrixXd A,
        Eigen::MatrixXd G,
        Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    Eigen::MatrixXd pars;
    label parsSamplesNum = initSnapInd.size();
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum;
    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum;
    newCoeffs[0].resize(newRowsNum, newColsNum);
    newCoeffs[1].resize(newRowsNum, G.cols());

    for (label j = 0; j < parsSamplesNum; j++)
    {
        Eigen::MatrixXd b0 = A.middleRows(j* timeSnapshotsPerSample,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd b2 = A.middleRows(j* timeSnapshotsPerSample + 1,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b2.cols());
        bNew << b2, (b2 - b0) / (timeSnap(j, 0));
        newCoeffs[0].block(j* timeSnapshotsPerSample - j, 0,
                           timeSnapshotsPerSample - 1, newColsNum) = bNew;
        newCoeffs[1].middleRows(j* timeSnapshotsPerSample - j,
                                timeSnapshotsPerSample - 1) = G.middleRows(j* timeSnapshotsPerSample + 1,
                                    timeSnapshotsPerSample - 1);
    }

    interChoice = 3;
    return newCoeffs;
}

List < Eigen::MatrixXd > UnsteadyNSTurb::velParCoeff(Eigen::MatrixXd A,
        Eigen::MatrixXd G)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    Eigen::MatrixXd pars;
    pars = z.leftCols(z.cols() - 1);
    newCoeffs[0].resize(A.rows(), A.cols() + z.cols() - 1);
    newCoeffs[1].resize(G.rows(), G.cols());
    newCoeffs[0] << pars, A;
    newCoeffs[1] = G;
    interChoice = 2;
    return newCoeffs;
}

List < Eigen::MatrixXd > UnsteadyNSTurb::velParDerivativeCoeff(
    Eigen::MatrixXd A, Eigen::MatrixXd G,
    Eigen::VectorXd initSnapInd, Eigen::VectorXd timeSnap)
{
    List < Eigen::MatrixXd > newCoeffs;
    newCoeffs.setSize(2);
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    Eigen::MatrixXd pars;
    pars = z.leftCols(z.cols() - 1);
    label parsSamplesNum = initSnapInd.size();
    label timeSnapshotsPerSample = snapshotsNum / parsSamplesNum;
    label newColsNum = 2 * velCoeffsNum;
    label newRowsNum = snapshotsNum - parsSamplesNum;
    newCoeffs[0].resize(newRowsNum, newColsNum + z.cols() - 1);
    newCoeffs[1].resize(newRowsNum, G.cols());

    for (label j = 0; j < parsSamplesNum; j++)
    {
        Eigen::MatrixXd b0 = A.middleRows(j* timeSnapshotsPerSample,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd b2 = A.middleRows(j* timeSnapshotsPerSample + 1,
                                          timeSnapshotsPerSample - 1);
        Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b2.cols());
        bNew << b2, (b2 - b0) / (timeSnap(j, 0));
        newCoeffs[0].block(j* timeSnapshotsPerSample - j, 0,
                           timeSnapshotsPerSample - 1, z.cols() - 1) = pars.middleRows(
                                   j* timeSnapshotsPerSample + 1,
                                   timeSnapshotsPerSample - 1);
        newCoeffs[0].block(j* timeSnapshotsPerSample - j, z.cols() - 1,
                           timeSnapshotsPerSample - 1, newColsNum) = bNew;
        newCoeffs[1].middleRows(j* timeSnapshotsPerSample - j,
                                timeSnapshotsPerSample - 1) = G.middleRows(j* timeSnapshotsPerSample + 1,
                                    timeSnapshotsPerSample - 1);
    }

    interChoice = 4;
    return newCoeffs;
}

Eigen::MatrixXd UnsteadyNSTurb::velParDerivativeCoeff(Eigen::MatrixXd A,
        Eigen::VectorXd par, double timeSnap)
{
    Eigen::MatrixXd newCoeffs;
    label velCoeffsNum = A.cols();
    label snapshotsNum = A.rows();
    label parsSamplesNum = par.size();
    label newColsNum = 2 * velCoeffsNum + parsSamplesNum;
    label newRowsNum = snapshotsNum - 1;
    newCoeffs.resize(newRowsNum, newColsNum);
    Eigen::MatrixXd b0 = A.topRows(A.rows() - 1);
    Eigen::MatrixXd b1 = A.bottomRows(A.rows() - 1);
    Eigen::MatrixXd bNew(b0.rows(), b0.cols() + b1.cols());
    bNew << b1, ((b1 - b0) / (timeSnap));
    newCoeffs.leftCols(parsSamplesNum) = Eigen::MatrixXd::Ones(newRowsNum,
                                         parsSamplesNum) * par;
    newCoeffs.rightCols(newColsNum - parsSamplesNum) = bNew;
    return newCoeffs;
}
