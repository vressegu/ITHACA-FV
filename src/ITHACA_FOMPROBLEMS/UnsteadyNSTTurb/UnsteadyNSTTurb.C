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

#include "UnsteadyNSTTurb.H"
#include "viscosityModel.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"


/// \file
/// Source file of the UnsteadyNSTTurb class.

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Construct Null
UnsteadyNSTTurb::UnsteadyNSTTurb() {};

// Construct from zero
UnsteadyNSTTurb::UnsteadyNSTTurb(int argc, char* argv[])
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
    /// Number of temperature modes to be calculated
    NTmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesTout", 15);
    /// Number of pressure modes to be calculated
    NPmodesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesPout", 15);
    /// Number of nut modes to be calculated
    NNutModesOut = para->ITHACAdict->lookupOrDefault<label>("NmodesNutOut", 15);
    /// Number of velocity modes used for the projection
    NUmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
    /// Number of temperature modes used for the projection
    NTmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesTproj", 10);
    /// Number of supremizers modes used for the projection
    NSUPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesSUPproj", 10);
    /// Number of pressure modes used for the projection
    NPmodes = para->ITHACAdict->lookupOrDefault<label>("NmodesPproj", 10);
    /// Number of nut modes used for the projection
    NNutModes = para->ITHACAdict->lookupOrDefault<label>("NmodesNutProj", 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void UnsteadyNSTTurb::truthSolve(List<scalar> mu_now, std::string& offlinepath)
{
    Time& runTime = _runTime();
    surfaceScalarField& phi = _phi();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& T = _T();
    volScalarField& nut = _nut();
    dimensionedScalar nu = _nu();
    dimensionedScalar Pr = _Pr();
    dimensionedScalar Prt = _Prt();
    volScalarField alphat = _alphat();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    IOMRFZoneList& MRF = _MRF();
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
            #include "TEqn.H"
            //TEqn.solve();
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
            ITHACAstream::exportSolution(T, name(counter), offlinepath);
            ITHACAstream::exportSolution(alphat, name(counter), offlinepath);
            std::ofstream of(offlinepath + name(counter) + "/" +
                             runTime.timeName());
            Ufield.append(tmp<volVectorField>(U));
            Pfield.append(tmp<volScalarField>(p));
            nutFields.append(tmp<volScalarField>(nut));
            Tfield.append(tmp<volScalarField>(T));
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

    if (mu_samples.rows() == nsnapshots * mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
                                   offlinepath);
    }
}

bool UnsteadyNSTTurb::checkWrite(Time& timeObject)
{
    scalar diffnow = mag(nextWrite - atof(timeObject.timeName().c_str()));
    scalar diffnext = mag(nextWrite - atof(timeObject.timeName().c_str()) -
                          timeObject.deltaTValue());

    if ( diffnow < diffnext)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void UnsteadyNSTTurb::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField& T = _T();
        simpleControl simple(mesh);
        //IOMRFZoneList& MRF = _MRF();
        //singlePhaseTransportModel& laminarTransport = _laminarTransport();
        //volScalarField& nut = _nut();
        volScalarField& alphat = _alphat();
        //dimensionedScalar& nu = _nu();
        dimensionedScalar& Pr = _Pr();
        dimensionedScalar& Prt = _Prt();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        scalar t1 = 1;
        scalar t0 = 0;

        for (label j = 0; j < T.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Tlift, j, t1);
            }
            else if (T.boundaryField()[BCind].type() == "zeroGradient")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t1);
            }
            else if (T.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t0);
            }
            else
            {
            }
        }

        while (simple.correctNonOrthogonal())
        {
            alphat = turbulence->nut() / Prt;
            alphat.correctBoundaryConditions();
            volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);
            fvScalarMatrix TEqn
            (
                fvm::ddt(Tlift)
                - fvm::laplacian(alphaEff, Tlift)
            );
            TEqn.solve();
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        scalar totalTime = mesh.time().value();
        scalar dt = mesh.time().deltaTValue();
        forAll(Tlift, i)
        {
            Tlift[i] = (totalTime * Tlift[i] + dt * Tlift[i] ) / (totalTime + dt);
        }
        Tlift.write();
        liftfieldT.append(Tlift.clone());
    }
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::turbulenceTerm1(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> CT1_matrix;
    CT1_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT1_matrix[j].resize(Nnutmodes, Csize);
        CT1_matrix[j] = CT1_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT1_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT1_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::laplacian(
                                          nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "eigen",
                               "./ITHACAoutput/Matrices/CT1");
    return CT1_matrix;
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::turbulenceTerm2(label NUmodes,
        label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List <Eigen::MatrixXd> CT2_matrix;
    CT2_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT2_matrix[j].resize(Nnutmodes, Csize);
        CT2_matrix[j] = CT2_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT2_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT2_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & (fvc::div(
                                          nuTmodes[j] * dev((fvc::grad(Together[k]))().T())))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "eigen",
                               "./ITHACAoutput/Matrices/CT2");
    return CT2_matrix;
}

Eigen::MatrixXd UnsteadyNSTTurb::BTturbulence(label NUmodes, label NSUPmodes)
{
    label BTsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd BT_matrix(BTsize, BTsize);
    BT_matrix = BT_matrix * 0;
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }

    // Project everything
    for (label i = 0; i < BTsize; i++)
    {
        for (label j = 0; j < BTsize; j++)
        {
            BT_matrix(i, j) = fvc::domainIntegrate(Together[i] & (fvc::div(dev((T(fvc::grad(
                    Together[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "eigen",
                               "./ITHACAoutput/Matrices/");
    return BT_matrix;
}

List <Eigen::MatrixXd> UnsteadyNSTTurb::temperatureTurbulenceTerm(
    label NTmodes, label Nnutmodes)
{
    label Stsize = NTmodes + liftfieldT.size();
    List <Eigen::MatrixXd> S_matrix;
    S_matrix.setSize(Stsize);

    for (label j = 0; j < Stsize; j++)
    {
        S_matrix[j].resize(Nnutmodes, Stsize);
    }

    PtrList<volScalarField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            Together.append(liftfieldT[k].clone());
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            Together.append(Tmodes[k].clone());
        }
    }

    for (label i = 0; i < Stsize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix S_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Stsize; k++)
            {
                S_matrix[i](j, k) = fvc::domainIntegrate(Together[i] * fvc::laplacian(
                                        nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "python",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "matlab",
                               "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "eigen",
                               "./ITHACAoutput/Matrices/S");
    return S_matrix;
}

void UnsteadyNSTTurb::projectSUP(fileName folder, label NU, label NP,
                                 label NSUP, label Nnut, label NT)
{
    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.npy");
        C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C", "C");
        K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.npy");
        P_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/P_mat.npy");
        M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.npy");
        BT_matrix =
            ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BT_matrix_mat.npy");
        CT1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT1",
                                              "CT1_matrix");
        CT2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT2",
                                              "CT2_matrix");
        Q_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/Q", "Q");
        Y_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/Y_mat.npy");
        S_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/S_mat.npy");
        MT_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/MT_mat.npy");
    }
    else
    {
        NUmodes = NU;
        NPmodes = NP;
        NSUPmodes = NSUP;
        Nnutmodes = Nnut;
        NTmodes = NT;
        B_matrix = diffusive_term(NUmodes, NPmodes, NSUPmodes);
        C_matrix = convective_term(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPmodes, NSUPmodes);
        M_matrix = mass_term(NUmodes, NPmodes, NSUPmodes);
        BT_matrix = BTturbulence(NUmodes, NSUPmodes);
        CT1_matrix = turbulenceTerm1(NUmodes, NSUPmodes, Nnutmodes);
        CT2_matrix = turbulenceTerm2(NUmodes, NSUPmodes, Nnutmodes);
        S_matrix = temperatureTurbulenceTerm(NTmodes, Nnutmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        MT_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
    }

    B_total_matrix = B_matrix + BT_matrix;
    C_total_matrix.setSize(C_matrix.size());

    for (label i = 0; i < C_matrix.size(); i++)
    {
        C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
    }

    NUmodes = NU;
    NPmodes = NP;
    NSUPmodes = NSUP;
    Nnutmodes = Nnut;
    // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    Eigen::MatrixXd Ncoeff = ITHACAutilities::getCoeffs(nutFields, nuTmodes);
    ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python",
                               "./ITHACAoutput/Matrices/");
    samples.resize(Nnutmodes);
    rbfsplines.resize(Nnutmodes);

    for (label i = 0; i < Nnutmodes; i++)
    {
        samples[i] = new SPLINTER::DataTable(1, 1);

        for (label j = 0; j < Ncoeff.cols(); j++)
        {
            samples[i]->addSample(mu.row(j), Ncoeff(i, j));
        }

        rbfsplines[i] = new SPLINTER::RBFSpline(*samples[i],
                                                SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    }
}
