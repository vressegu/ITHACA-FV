#include "ITHACAparameters.H"
#include "ITHACAutilities.H"
#include "Foam2Eigen.H"

template<class Enum>
struct StrLookup: public std::map<std::string,Enum>{
  using Base=std::map<std::string,Enum>;
  using Base::Base;
  Enum lookup(std::string str,Enum defaultValue){
    if (this->count(str)){
      return Base::at(str);
    }else{
      return defaultValue;
    }
  }
};



ITHACAparameters* ITHACAparameters::instance = nullptr;

///CHC new POD specific constructor 
ITHACAparameters::ITHACAparameters(Foam::fvMesh& mesh, Time& localTime, const word& myfield_name) :
        runTime(localTime),
        mesh(mesh),
        corTime
            (
                "corTime",
                dimensionSet(0,0,1,0,0), 
                scalar(0)
            ),
        myfield_name(myfield_name)
{
    Info << "=============>>>>>     ITHACAParameters constructor, step 0 : create dictionary and get init infos" << endl;

    //Add constructor code of "standard" ithacaParams
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            true
        )
    );
    precision = ITHACAdict->lookupOrDefault<label>("OutPrecision", 10);
    word typeout = ITHACAdict->lookupOrDefault<word>("OutType", "fixed");

    if (typeout == "fixed")
    {
        outytpe = std::ios_base::fixed;
    }

    if (typeout == "scientific")
    {
        outytpe = std::ios_base::scientific;
    }

    eigensolver = ITHACAdict->lookupOrDefault<word>("EigenSolver", "spectra");
    exportPython = ITHACAdict->lookupOrDefault<bool>("exportPython", 0);
    exportMatlab = ITHACAdict->lookupOrDefault<bool>("exportMatlab", 0);
    exportTxt = ITHACAdict->lookupOrDefault<bool>("exportTxt", 0);
    debug = ITHACAdict->lookupOrDefault<bool>("debug", 0);
    warnings = ITHACAdict->lookupOrDefault<bool>("warnings", 0);

  //Add to constructor code of "RedLUM" IthacaFVParams
    Info << "=============>>>>>     ITHACAParameters constructor, step 1 : get runTime0" << endl;
    runTime0 = autoPtr<Foam::Time>(&runTime) ;

    Info << "=============>>>>>     ITHACAParameters constructor, step 2 : get nb cells in mesh" << endl;
    nCells = mesh.cells().size();
    Info << "=============>>>>>     mesh.cells().size() = "<< nCells << endl;  

    Info << "=============>>>>>     ITHACAParameters constructor, step 2.1 : get casenameData in dict" << endl;        
    casenameData = ITHACAdict->lookupOrDefault<fileName>("casename", "./");
    Info << "=============>>>>>     ITHACAParameters casenameData=" << casenameData << endl;

    Info << "=============>>>>>     ITHACAParameters constructor, step 2.2 : find 'fields' in dictionary" << endl;

//#define _REDLUM_INIT_FOR_PARAM
#ifdef _REDLUM_INIT_FOR_PARAM

    fieldlist = static_cast<List<word>>(ITHACAdict->lookup("fields"));


    Info << "=============>>>>>     ITHACAParameters constructor, step 2.3 !" << endl;     
    pressureResolutionKind = StrLookup<PressureResolutionKind>(
          {
            {"FullOrder",PressureResolutionKind::FullOrder},
            {"ReducedOrder",PressureResolutionKind::ReducedOrder},
            {"Neglected",PressureResolutionKind::Neglected},
          }).lookup(ITHACAdict->lookupOrDefault<word>("pressureResolutionKind", ""),PressureResolutionKind::Undefined);
    nBlocks = ITHACAdict->lookupOrDefault<label>("nBlocks", 1);
    advModifOrNot = ITHACAdict->lookupOrDefault<bool>("advModifOrNot", 1);
    centeredOrNot = ITHACAdict->lookupOrDefault<bool>("centeredOrNot", 1);
    interpFieldCenteredOrNot = ITHACAdict->lookupOrDefault<bool>("interpFieldCenteredOrNot", 0);
    // nMagicPoints = ITHACAdict->lookupOrDefault<label>("nMagicPoints", 1);
    DEIMInterpolatedField = ITHACAdict->lookupOrDefault<word>("DEIMInterpolatedField", "fullStressFunction");
    assimilation = ITHACAdict->lookupOrDefault<bool>("assimilation", 0);
    onLineReconstruct = ITHACAdict->lookupOrDefault<bool>("onLineReconstruct", 0);
    forcingOrNot = ITHACAdict->lookupOrDefault<bool>("forcingOrNot", 0);
    symDiff = ITHACAdict->lookupOrDefault<bool>("symDiff", 0);
    forcingHasFfieldOrNot = ITHACAdict->lookupOrDefault<bool>("forcingHasFfieldOrNot", 0);
    stochasOrNot = ITHACAdict->lookupOrDefault<bool>("stochasOrNot", 1);
    useStratonovich = ITHACAdict->lookupOrDefault<bool>("useStratonovich", 0);
    ROMTemporalScheme = ITHACAdict->lookupOrDefault<word>("ROMTemporalScheme", "euler");
    useHypRedSto = ITHACAdict->lookupOrDefault<bool>("HypRedSto", 0);
    inflatNut = ITHACAdict->lookupOrDefault<bool>("inflatNut", 0);
           // SOTA can be 0 (no SOTA), D (deterministic version), S (stochastic version)
    useSOTA = ITHACAdict->lookupOrDefault<word>("useSOTA", "None");
    set_useSOTA(useSOTA);
    Info << "=============>>>>>     ITHACAParameters constructor, step 3 !" << endl;

  if (useSOTA != "None")
  {
    Info << "===============================================================" << endl;
    Info << "  RedLUM is launched in " << useSOTA << "-SOTA mode" << endl;
    Info << "  Ithacadict file will be overidden by parameters specified in IthacaFVParameters.C" << endl;
    Info << "===============================================================" << endl;

    if (useSOTA == "D")
    {
    advModifOrNot = 0;
    centeredOrNot = 1;
    DEIMInterpolatedField =  "fullStressFunction";
    assimilation = 0;
    stochasOrNot = 0;
    useHypRedSto = 0;
    inflatNut = 0;
    useStratonovich=0;
    }
    
  }
   Info << "=============>>>>>     ITHACAParameters constructor, step 4 !" << endl;
   
    if ( (!(ROMTemporalScheme == "adams-bashforth"))
         && (!(ROMTemporalScheme == "euler"))
         && (!(ROMTemporalScheme == "euler–maruyama")) )
    {
      Info << "This temporal scheme is not implemented." << endl;
      abort();
    }
    if ( (ROMTemporalScheme == "adams-bashforth") && stochasOrNot )
    {
      Info << "This temporal scheme is deterministic." << endl;
      abort();
    }
    if ( useStratonovich && !stochasOrNot )
    {
      Info << "Stratonovich correction needs stochastic model to be performed." << endl;
      abort();
    }
    
    Info << "=============>>>>>     ITHACAParameters constructor, step 5 !" << endl;

    // Object Time to read OpenFOAM data in the correct folder
    runTimeData = new Foam::Time(Foam::Time::controlDictName, ".", casenameData);

    template_field_U = new volVectorField
        (
          IOobject
          (
            "U",
            runTimeData->path() + runTimeData->times()[1].name(),
            mesh,
            IOobject::MUST_READ
            ),
          mesh
          );

    template_field_p = new volScalarField
        (
          IOobject
          (
            "p",
            runTimeData->path() + runTimeData->times()[1].name(),
            mesh,
            IOobject::MUST_READ
            ),
          mesh
          );
Info << "=============>>>>>     ITHACAParameters constructor, step 6 !" << endl;

    nParticules = ITHACAdict->lookupOrDefault<label>("nParticules", 1);
     
    reducingNoisesMatrix = ITHACAdict->lookupOrDefault("reducingNoisesMatrix", true);
    nSimu = ITHACAdict->lookupOrDefault("nSimu", 100);
    intervalConfidenceRatio = ITHACAdict->lookupOrDefault("intervalConfidenceRatio", 0.95);


    // Initialize field_name, field_type and nModes
    field_name.resize(fieldlist.size());
    field_type.resize(fieldlist.size());
    Info << "=============>>>>>     ITHACAParameters constructor, step 7 !" << endl;

    for (label k = 0; k < fieldlist.size(); k++)
    {
      dictionary& subDict = ITHACAdict->subDict(fieldlist[k]);
      field_name[k] = static_cast<word>(subDict.lookup("field_name"));
      field_type[k] = static_cast<word>(subDict.lookup("field_type"));
      nModes.insert(field_name[k],subDict.lookupOrDefault<label>("nmodes",1));
      hilbertSpacePOD.insert(field_name[k],
                             subDict.lookupOrDefault<word>("hilbertSpacePOD","L2"));
      varyingEnergy.insert(field_name[k], 0);
      resolvedVaryingEnergy.insert(field_name[k], 0);
    }

    weightH1 = 1.0;
    weightBC = ITHACAdict->lookupOrDefault<double>("weightBC", 0.);
    patchBC = ITHACAdict->lookupOrDefault<word>("patchBC", "inlet");
    freqBC_0 = ITHACAdict->lookupOrDefault<double>("freqBC_0", 0);

    // Initialize startTime, endTime and nSnapshots

    // Get times list from the case folder
    instantList Times = runTimeData->times();
Info << "=============>>>>>     ITHACAParameters constructor, step 8 !" << endl;

    // Read Initial and last time from the POD dictionary
    const entry* existnsnap = ITHACAdict->findEntry("Nsnapshots");
    const entry* existLT = ITHACAdict->findEntry("FinalTime");

    scalar InitialTime(0);
    // Initiate variable from PODSolverDict
    if ((existnsnap) && (existLT))
    {
      Info << "Error you cannot define LatestTime and NSnapShots together" << endl;
      abort();
    }
    else if (existnsnap)
    {
      InitialTime = ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0);
      FinalTime = ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000);
      nSnapshots = readScalar(ITHACAdict->lookup("Nsnapshots"));
      startTime = Time::findClosestTimeIndex(runTimeData->times(), InitialTime);
      nSnapshots = min(nSnapshots , Times.size() - startTime);
      endTime = startTime + nSnapshots - 1;
      FinalTime = std::stof(runTimeData->times()[endTime].name());
    }
    else
    {
      InitialTime = ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0);
      FinalTime = ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000);
      endTime = Time::findClosestTimeIndex(runTimeData->times(), FinalTime);
      startTime = Time::findClosestTimeIndex(runTimeData->times(), InitialTime);
      nSnapshots = endTime - startTime + 1;
      if (InitialTime > FinalTime)
      {
        Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAdict file\n" << endl;
        abort();
      }
      FinalTime = std::stof(runTimeData->times()[endTime].name());
    }
  
    // Read Initial and last time from the POD dictionary
    const entry* existnsnapSimulation = ITHACAdict->findEntry("NsnapshotsSimulation");
    const entry* existLTSimulation = ITHACAdict->findEntry("FinalTimeSimulation");

    scalar InitialTimeSimulation(FinalTime);
    label startTimeSimulation(Time::findClosestTimeIndex(runTimeData->times(), InitialTimeSimulation));
Info << "=============>>>>>     ITHACAParameters constructor, step 9 !" << endl;

    if ((existnsnapSimulation) && (existLTSimulation))
    {
      Info << "Error you cannot define LatestTimeSimulation and NSnapShotsSimulation together" << endl;
      abort();
    }
    else if (existnsnapSimulation)
    {
      nSnapshotsSimulation = readScalar(ITHACAdict->lookup("NsnapshotsSimulation"));
      nSnapshotsSimulation = min(nSnapshotsSimulation , Times.size() - startTimeSimulation);
      endTimeSimulation = startTimeSimulation + nSnapshotsSimulation - 1;
      FinalTimeSimulation = std::stof(runTimeData->times()[endTimeSimulation].name());
    }
    else
    {
      FinalTimeSimulation = ITHACAdict->lookupOrDefault<scalar>("FinalTimeSimulation", 100000000000000);
      endTimeSimulation = Time::findClosestTimeIndex(runTimeData->times(), FinalTimeSimulation);
      nSnapshotsSimulation = endTimeSimulation - startTimeSimulation + 1;
      if (InitialTimeSimulation > FinalTimeSimulation)
      {
        Info << "FinalTimeSimulation cannot be smaller than the InitialTimeSimulation check your ITHACAdict file\n" << endl;
        abort();
      }
      FinalTimeSimulation = std::stof(runTimeData->times()[endTimeSimulation].name());
    }
Info << "=============>>>>>     ITHACAParameters constructor, step 10 !" << endl;

        // Initialize saveTime
    IOdictionary controlDict
        (
          IOobject
          (
            "controlDict",
            runTimeData->caseName() + runTimeData->system(),
            get_mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
          );
    saveTime = controlDict.lookupOrDefault<double>("writeInterval", 1);
    corTime.value() = saveTime;
Info << "=============>>>>>     ITHACAParameters constructor, step 11 !" << endl;

    IOdictionary transportProperties
        (
          IOobject
          (
            "transportProperties",
            runTimeData->caseName() + runTimeData->constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
          );

    nu = new dimensionedScalar
        (
          "nu",
          dimViscosity,
          transportProperties
          );
Info << "=============>>>>>     ITHACAParameters constructor, step 12 !" << endl;

    // if TAG forcing term taken into account, set all TAG related values to their
    //value
    if(forcingOrNot)
    {
      if(forcingHasFfieldOrNot)
      {
        Info<<"Entered forcing-specific value assignment"<<endl;
        freqRotActuDisk = new dimensionedScalar
            (
              "freqRotActuDisk",
              dimless,
              transportProperties
              );

        centerActuDisk = new dimensionedVector
            (
              "centerActuDisk",
              dimless,
              transportProperties
              );

        lengthActuDisk = new dimensionedScalar
            (
              "lengthActuDisk",
              dimless,
              transportProperties
              );

        widthActuDisk = new dimensionedScalar
            (
              "widthActuDisk",
              dimless,
              transportProperties
              );

        initPhaseActuDisk = new dimensionedScalar
            (
              "initPhaseActuDisk",
              dimless,
              transportProperties
              );

        timeBeforeRota = new dimensionedScalar
            (
              "timeBeforeRota",
              dimless,
              transportProperties
              );

        magForcage = new dimensionedScalar
            (
              "magForcage",
              dimless,
              transportProperties
              );

        simu = new dimensionedScalar
            (
              "simu",
              dimless,
              transportProperties
              );
        
        compute_F_mask();
      }
      // if TAG forcing term taken into account, set all TAG related values to 0
      else
      {
        freqRotActuDisk = new dimensionedScalar();

        centerActuDisk = new dimensionedVector();

        lengthActuDisk = new dimensionedScalar();

        widthActuDisk = new dimensionedScalar();

        initPhaseActuDisk = new dimensionedScalar();

        timeBeforeRota = new dimensionedScalar();

        magForcage = new dimensionedScalar();

        simu = new dimensionedScalar();
      }
    }
Info << "=============>>>>>     ITHACAParameters constructor, step 13 !" << endl;
        volume = new volScalarField(
          IOobject(
            "volume",
            runTimeData->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
          mesh,
          dimensionedScalar("volume", dimensionSet(0, 3, 0, 0, 0, 0, 0), 1.0)
          );
    volume->primitiveFieldRef() = mesh.V().field();
    Eigen::VectorXd volVect = Foam2Eigen::field2Eigen(mesh.V().field());
    totalVolume = volVect.sum();

    delta = new volScalarField(pow(*volume,1.0/3.0));

    //TO DO : rewrite the following method to search for the object turbulenceProperties and its attibutes (simulationType and LESModel)
    string simulationType;
    string LESModel;

    string DDESModel;

Info << "=============>>>>>     ITHACAParameters constructor, step 14 !" << endl;

    string filename = "constant/turbulenceProperties";
    std::ifstream strm( filename );

    string line;
    getline( strm, line );

    while (line.find( "simulationType" ) == string::npos || line.find( "//" ) != string::npos)
    {
      getline( strm, line );
    }

    int N = line.size();
    int i =0;
    while ( i<N )
    {
      if (line[i]==' ')
      {
        line.erase(i,1);
        N=N-1;
      }
      else
      {
        i++;
      }
    }
Info << "=============>>>>>     ITHACAParameters constructor, step 15 !" << endl;

    line.erase(line.size()-1,1);
    line.erase(0,14);
    simulationType=line;



    Info << "---------------------------------------------------------" << endl;
    Info << "Simulation type:" << simulationType << endl;


    if (simulationType=="laminar")
    {
      Info<< "DNS simulation will be performed" <<endl;
      set_useDNS(true);
      useHypRedSto = false;
    }


    else if(simulationType =="LES")
    {
      template_field_nut = new volScalarField
          (
            IOobject
            (
              "nut",
              runTimeData->path() + runTimeData->times()[1].name(),
              mesh,
              IOobject::MUST_READ
              ),
            mesh
            );

      getline( strm, line );
      while (line.find( "LESModel" ) == string::npos || line.find( "//" ) != string::npos )
      {
        getline( strm, line );
      }

      N = line.size();
      i =0;
      while ( i<N )
      {
        if (line[i]==' ')
        {
          line.erase(i,1);
          N=N-1;
        }
        else
        {
          i++;
        }
      }
      line.erase(N-1,1);
      line.erase(0,8);
      LESModel=line;

      Info << "LESmodel:" << LESModel << endl;

      if (LESModel=="Smagorinsky")
      {
        Info << "Simulation type LES Smagorinsky, DEIM will be performed on non linear term" << endl;
        Info << "DEIM will be perfomed on " << DEIMInterpolatedField << " terms." << endl;
        set_useDEIM(true);
        set_Ck(ITHACAdict->lookupOrDefault<double>("Ck", 0.094));
        set_Ce(ITHACAdict->lookupOrDefault<double>("Ce", 1.048));
        set_useDDES(false);
        if (useHypRedSto)
        { Info <<" Stochastic Hyperreduction  will be perfomed" << endl;
        } else{Info <<"Stochastic Hyperreduction will not be perfomed" << endl;}


      }
      else if (LESModel=="kOmegaSSTDDES")
      { Info<< "New Simulationtype of DDES model will be performed"<<endl;
        set_useDDES(true);
        set_useDEIM(false);
        template_field_omega = new volScalarField
            (
              IOobject
              (
                "omega",
                runTimeData->path() + runTimeData->times()[1].name(),
                mesh,
                IOobject::MUST_READ
                ),
              mesh
              );
        template_field_k = new volScalarField
            (
              IOobject
              (
                "k",
                runTimeData->path() + runTimeData->times()[1].name(),
                mesh,
                IOobject::MUST_READ
                ),
              mesh
              );
      }
      else
      {
        Info << "Only turbulence model Smagorinsky is supported with simulation type LES" << endl;
        Info << "DDES NOT USED**0**" << endl;


      }


    }

    else{
      Info <<"Simulation type not LES, no DEIM will be performed" << endl;
      set_useDEIM(false);
      Info<<"\n"<<"No DDES was  performed"<<endl;
      set_useDDES(false);
      Info <<"\n"<<"No DNS simulation was perfomed"<<endl;
      set_useDNS(false);
      abort();
    }
    if (useDEIM)
    {
      if(get_hilbertSpacePOD()["nut"] == "dL2"){
        set_deltaWeight(ITHACAutilities::getMassMatrixFV(*template_field_nut).array().pow(-2.0/3.0));
      }
      nMagicPoints = get_nModes()[DEIMInterpolatedField];
      if (interpFieldCenteredOrNot)
      {
        folder_DEIM = "./ITHACAoutput/DEIM_centered/";
        if ((DEIMInterpolatedField == "nut") && useHypRedSto && !useStratonovich)
        {
          Info << "Error: The stochastic hyperreduction requires a Stratonovich integration scheme." << endl;
          abort();
        }
      }
      else
      {
        folder_DEIM = "./ITHACAoutput/DEIM/";
        if (useHypRedSto)
        {
          Info << "Error: The stochastic hyperreduction requires a centered DEIM." << endl;
          abort();
        }
      }
      folder_DEIM += DEIMInterpolatedField + "/";
      folder_DEIM += std::to_string(nMagicPoints) + "magicPoints/";
    }
    else
    {
      folder_DEIM = "./ITHACAoutput/";
    }

    Info << "---------------------------------------------------------" << endl;
    

#endif
}



ITHACAparameters* ITHACAparameters::getInstance(fvMesh& mesh, Time& localTime, const word& myfield_name)
{
    if (instance == nullptr)
    {
      Info << "=============>>>>>     create new ITHACAParameters !" << endl;
        instance = new ITHACAparameters(mesh, localTime, myfield_name);
    }
Info << "=============>>>>>     return existing ITHACAParameters instance!" << endl;
    return instance;
}


ITHACAparameters* ITHACAparameters::getInstance()
{
    M_Assert(instance != nullptr,
             "ITHACAparameters needs to be initialized, call ITHACAparameters::getInstance(mesh, runTime) first");
    return instance;
}



void ITHACAparameters::compute_F_mask()
{
    volScalarField& result_F_mask = *(new volScalarField
                                      (
                                        IOobject
                                        (
                                          "p",
                                          runTimeData->path() + runTimeData->times()[1].name(),
                                          mesh,
                                          IOobject::MUST_READ
                                          ),
                                        mesh
                                        ));

    if(forcingHasFfieldOrNot)
    {
      Info << "Computing F mask from characteristics values" << endl;

      forAll ( result_F_mask, i )
      {
        result_F_mask[i] = 0.;
        Foam::Vector<double> normalActuDisk = get_F_0();

        // Getting coordinates of current grid cell.

        double x = result_F_mask.mesh().C()[i][0];
        double y = result_F_mask.mesh().C()[i][1];
        double z = result_F_mask.mesh().C()[i][2];

        // Calculating the vector going from the center of the actuator disk to current point. Used to check if current point
        // is part of the actuator disk field F.

        Foam::Vector<double> centerToCurrentPoint
            (
              (centerActuDisk->value()[0]-x),
              (centerActuDisk->value()[1]-y),
              (centerActuDisk->value()[2]-z)
              );

        // Calculating the coordinates of the previous vector in rotating vector coordinates.

        double yRota = centerToCurrentPoint[0]*normalActuDisk[0]
            + centerToCurrentPoint[1]*normalActuDisk[1]
            + centerToCurrentPoint[2]*normalActuDisk[2];

        double rRota = mag(centerToCurrentPoint-yRota*normalActuDisk);


        if ( std::abs(yRota) < widthActuDisk->value() && rRota < lengthActuDisk->value() )
        {
          result_F_mask[i]=1;
        }
      }
    }
    else
    {
      Info << "Computing F_mask from geometry" << endl;

      forAll ( result_F_mask, i )
      {
        result_F_mask[i] = 0;
      }

      string filename = "constant/polyMesh/sets/actuationDisk1";
      std::ifstream strm( filename );

      string line;
      string lineTemp;
      getline( strm, line );
      lineTemp = line;

      while (line.find( "(" ) == string::npos || line.find( "//" ) != string::npos)
      {
        lineTemp = line;
        getline( strm, line );
      }

      label size = stoi(lineTemp);
      for(int i = 0; i < size; i++)
      {
        getline( strm, line );
        label current_cell = stoi(line);

        result_F_mask[current_cell] = 1;
      }
    }

    template_F_mask = &result_F_mask;
  }

  const Foam::Vector<double>& ITHACAparameters::get_F_0() const
  {
    Foam::Vector<double>& result_F_0 = *(new Foam::Vector<double>
                                         (
                                           std::cos(initPhaseActuDisk->value()),
                                           std::sin(initPhaseActuDisk->value()),
                                           0.
                                           ));

    return(result_F_0);
  }



  //}
