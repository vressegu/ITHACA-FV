#ifndef ITHACA_FV_LERAY_PROJECTION
#define ITHACA_FV_LERAY_PROJECTION

#include "ITHACAparameters.H"
#include "StoredParameters.H"


  // ----------------------------------------------
  // Methods to perform free-divergence projection
  // and free-divergence projection orthogonal
  // by solving Poisson's equation in full-order space



namespace ITHACAutilities
{

class ITHACAlerayProj
{
  private :
    ITHACAparameters* para;
    autoPtr<volScalarField> solution;
    std::unique_ptr<volScalarField> initialGuess;
  public :
    bool export_debug;
    dimensionedScalar nu_fake;
    ITHACAlerayProj(const fvMesh& mesh);
    const volVectorField freeDivergenceProjectionOrthogonal(const volVectorField& function);
    volVectorField freeDivergenceProjection(volVectorField& function);
};

}

#endif
