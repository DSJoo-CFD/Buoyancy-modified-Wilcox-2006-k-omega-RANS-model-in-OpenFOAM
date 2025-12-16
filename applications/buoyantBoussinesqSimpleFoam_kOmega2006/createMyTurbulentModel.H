/*---------------------------------------------------------------------------*\
Models for k, omega, and nut are implemented similarly to those in 
"LaunderSharmaKE.C" in OpenFOAM. The algebraic turbulent heat flux model and 
the temperature variance equation are implemented as described in the paper by 
Kenjereš, S. and Hanjalić, K., 1995, Int. J. Heat Fluid Flow, 16(5).

My new modeling approach only algebraically modifies the buoyancy-related term 
in the omega equation by incorporating a global effect of buoyancy, 
represented by available potential energy.
\*---------------------------------------------------------------------------*/

// Model transport equation for k (turbulent kinetic energy)
Info<< "Reading field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Model transport equation for omega
Info<< "Reading field omega\n" << endl;
volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Turbulent viscosity
volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    k/omega
);


volScalarField sigmad
(
    IOobject
    (
        "sigmad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimensionedScalar("zero",dimless,0.125)
);


volScalarField fbeta
(
    IOobject
    (
        "fbeta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimless
);

volScalarField absChi
(
    IOobject
    (
        "absChi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimless
);


volScalarField Cbetas_lowRe
(
    IOobject
    (
        "Cbetas_lowRe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimensionedScalar("zero",dimless,1.0)
);


volScalarField Calpha_lowRe
(
    IOobject
    (
        "Calpha_lowRe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimensionedScalar("zero",dimless,1.0)
);

volScalarField Calphas_lowRe
(
    IOobject
    (
        "Calphas_lowRe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh, dimensionedScalar("zero",dimless,1.0)
);




// volScalarField Cbetas_lowRe("Cbetas_lowRe",dimless,scalar(1.0));

scalar kMin(laminarTransport.getOrDefault<scalar>("kMin",1.e-12));
scalar kMax(laminarTransport.getOrDefault<scalar>("kMax",1.e+5));
scalar omegaMin(laminarTransport.getOrDefault<scalar>("omegaMin",1.e-12));
scalar omegaMax(laminarTransport.getOrDefault<scalar>("omegaMax",1.e+10));
dimensionedScalar kMin_(dimensionSet(0,2,-2,0,0,0,0),kMin);
dimensionedScalar omegaMin_(dimensionSet(0,0,-1,0,0,0,0),omegaMin);
dimensionedScalar kMax_(dimensionSet(0,2,-2,0,0,0,0),kMax);
dimensionedScalar omegaMax_(dimensionSet(0,0,-1,0,0,0,0),omegaMax);
bound(k, kMin_);
bound(omega, omegaMin_);


// Turbulent heat flux
volVectorField thf
(
    IOobject
    (
        "thf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U*T*scalar(0)
);

volVectorField Nu
(
    IOobject
    (
        "Nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U*T*scalar(0)
);

Info<< "Transport properties:" << endl;
Info<< "nu = " << nu.value() << endl;
Info<< "beta = " << beta.value() << endl;
Info<< "Q = " << Q << endl;
Info<< "Pr = " << Pr.value() << endl;
dimensionedScalar Ra(mag(g)*beta/nu/nu*Pr);
Info<< "Ra = " << Ra.value() << endl;
dimensionedScalar Rai(mag(g)*beta*Q/nu/nu/nu*Pr*Pr);
Info<< "Rai = " << Rai.value() << endl;
Info<< " " << endl;


// see Wilcox (2006) Turbulent modeling for CFD 3rd, Chapter 4.3.1 The k-w model, page 126-127
dimensionedScalar Clim(laminarTransport.getOrDefault<scalar>("Clim", 7.0/8.0));
dimensionedScalar Calpha(laminarTransport.getOrDefault<scalar>("Calpha", 13.0/25.0));
dimensionedScalar Cbeta0(laminarTransport.getOrDefault<scalar>("Cbeta0", 0.0708));
dimensionedScalar Cbetas(laminarTransport.getOrDefault<scalar>("Cbetas", 0.09));
dimensionedScalar SigmaK(laminarTransport.getOrDefault<scalar>("SigmaK", 0.5));
dimensionedScalar SigmaW(laminarTransport.getOrDefault<scalar>("SigmaW", 0.6));
dimensionedScalar Csigmad0(laminarTransport.getOrDefault<scalar>("Csigmad0", 0.125));
dimensionedScalar Cg(laminarTransport.getOrDefault<scalar>("Cg", -2.0));
dimensionedScalar Prt(laminarTransport.getOrDefault<scalar>("Prt", 0.89));

Switch lowReVersion(laminarTransport.getOrDefault<Switch>("lowReVersion",false));


Info<< "Wilcox (2006) k-omega Model" << endl;
Info<< "Wilcox, D. C., Turbulence Modeling for CFD, 3rd edition, 2006" << endl;
Info<< " " << endl;

Info<< "Modeling coefficients:" << endl;
Info<< "C_lim = " << Clim.value() << endl;
Info<< "alpha = " << Calpha.value() << endl;
Info<< "beta_0 = " << Cbeta0.value() << endl;
Info<< "beta_str = " << Cbetas.value() << endl;
Info<< "sigma = " << SigmaK.value() << endl;
Info<< "sigma_str = " << SigmaW.value() << endl;
Info<< "sigma_d0 = " << Csigmad0.value() << endl;

Info<< " " << endl;
Info<< "C_g = " << Cg.value() << endl; // buoyant production in omega equation
Info<< "Pr_t = " << Prt.value() << endl; // turbulent Prandtl number

Info<< " " << endl;
Info<< "low Reynolds Version activate: " << lowReVersion << endl;

Info<< " " << endl;
Info<< "kMin = " << kMin << endl;
Info<< "kMax = " << kMax << endl;
Info<< "omegaMin = " << omegaMin << endl;
Info<< "omegaMax = " << omegaMax << endl;