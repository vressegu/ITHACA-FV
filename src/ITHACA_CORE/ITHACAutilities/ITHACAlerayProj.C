#include "ITHACAlerayProj.H"
#include "ITHACAstream.H"
#include "IFstream.H"

namespace ITHACAutilities
{
ITHACAlerayProj::ITHACAlerayProj(const fvMesh& mesh)
{
  // IFstream ithacaDictPath("./system/ITHACAdict");
  // Foam::dictionary ithacaDict(ithacaDictPath);
  // fileName rootPath = ithacaDict.lookupOrDefault<Foam::fileName>("rootPath", "./");
  // fileName casename = ithacaDict.lookupOrDefault<Foam::fileName>("casename", "./");

  // const Time& runTime = mesh.time();
  Time& runTime = const_cast<Time&>(mesh.time());

  ITHACAparameters* para = ITHACAparameters::getInstance(mesh, runTime);

  // Initialize the solution field using the persistent mesh
    solution.reset(new volScalarField(
        IOobject(
            "p",
            runTime.times()[1].name(), 
            mesh,
            IOobject::MUST_READ
        ),
        mesh
    ));


    // Make a clean backup copy. This prevents the solver from using 
    // the previous timestep's projection as its initial guess.
    initialGuess.reset(new volScalarField(*solution));

}

const volVectorField ITHACAlerayProj::freeDivergenceProjectionOrthogonal(const volVectorField& function) 
{
  dimensionedScalar nu_fake(1);

  *solution = *initialGuess;
  

  fvScalarMatrix LapEqn
  (
      fvm::laplacian(nu_fake,*solution)
    ==
      fvc::div(function)
  );

  LapEqn.solve();

  return fvc::grad(*solution);

}

volVectorField ITHACAlerayProj::freeDivergenceProjection(volVectorField& function)
{
  Info << "" << endl;
  Info << "Performing projection onto the space of divergence-free functions" << endl;

  return (function - freeDivergenceProjectionOrthogonal(function));
}
}