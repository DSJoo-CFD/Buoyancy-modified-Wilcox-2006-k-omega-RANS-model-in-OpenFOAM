// Info<< "myTurbulentCorrect.H\n" << endl;

volScalarField k_(max(k,kMin_));
volScalarField omega_(max(omega,omegaMin_));


volScalarField divU(fvc::div(fvc::absolute(phi, U)));
volTensorField gradU(fvc::grad(U));
volVectorField gradT(fvc::grad(T));


volSymmTensorField Sij(symm(gradU));
// volTensorField Omij(-skew(gradU));
volScalarField StressLim(Clim*sqrt(2.0/Cbetas)*mag(Sij));

volSymmTensorField tauij(2.0*nut*Sij-((2.0/3.0)*I)*k);
volVectorField gradK(fvc::grad(k));
volVectorField gradOmega(fvc::grad(omega));

volScalarField G(tauij && gradU);
// volScalarField G(min(tauij && gradU, 20.0*Cbetas*omega*k)); // SST flux limiter: used for stabilization in the early stage of the simulation
volScalarField Gb(-beta*(g&thf)); // buoyant production





volScalarField sigmadCheck(gradK&gradOmega);
forAll(sigmad,celli)
{
   if (sigmadCheck[celli] <= 0.001)
   {
      sigmad[celli]=scalar(0.0);
   }else
   {
      sigmad[celli]=scalar(0.125);
   }
}

if (lowReVersion)
{
   // Info<< "low Re coefficient tuning on\n" << endl;

   Cbetas_lowRe = (100.0*Cbeta0/27.0+pow(k/omega/nu/8.0,4))/(1+pow(k/omega/nu/8.0,4));
   Calphas_lowRe = (Cbeta0/3.0+k/omega/nu/6.0)/(1.0+k/omega/nu/6.0);
   Calpha_lowRe = (13.0/25.0)*(1.0/9.0+k/omega/nu/2.61)/(1.0+k/omega/nu/2.61)/Calphas_lowRe;

   nut = Calphas_lowRe*k/max(omega,Clim*sqrt(2.0/Cbetas*Calphas_lowRe*magSqr(Sij)));

} else
{
   // Info<< "low Re coefficient tuning off\n" << endl;

   Cbetas_lowRe = scalar(1.0);
   Calphas_lowRe = scalar(1.0);
   Calpha_lowRe = scalar(1.0);

   nut = k/max(omega,Clim*sqrt(2.0/Cbetas*magSqr(Sij)));
}


tmp<fvScalarMatrix> omegaEqn
(
    fvm::ddt(omega)
  + fvm::div(phi, omega)
  - fvm::laplacian(nu+Calphas_lowRe*SigmaW*k_/omega_,omega)
 ==
    Calpha*omega/k_*(G + Cg*max(Gb,scalar(0)*kMin_*omegaMin_) -2.0*min(Gb,scalar(0)*kMin_*omegaMin_) )*Calpha_lowRe
    // Calpha*omega/k_*(G + Cg*Gb)*Calpha_lowRe
  - fvm::Sp(Cbeta0*omega,omega)
  // - fvm::Sp(Cbeta0*fbeta*omega,omega) // fbeta = 1 for 2D
  + sigmad/omega*(gradK&gradOmega)
);
omegaEqn.ref().relax();
omegaEqn.ref().boundaryManipulate(omega.boundaryFieldRef());
solve(omegaEqn);
bound(omega, omegaMin_);


tmp<fvScalarMatrix> kEqn
(
    fvm::ddt(k)
  + fvm::div(phi, k)
  - fvm::laplacian(nu+Calphas_lowRe*SigmaK*k/omega,k)
 ==
    G + Gb
 - fvm::Sp(Cbetas*Cbetas_lowRe*omega, k)
);
kEqn.ref().relax();
solve(kEqn);
bound(k, kMin_);


k=min(k,kMax_);
omega=min(omega,omegaMax_);


if (lowReVersion)
{
// Info<< "low Reynolds Version activate: " << lowReVersion << endl;
   nut = Calphas_lowRe*k/max(omega,Clim*sqrt(2.0/Cbetas*Calphas_lowRe*magSqr(Sij)));
} else
{
   nut = k/max(omega,Clim*sqrt(2.0/Cbetas*magSqr(Sij)));
}



thf = -nut/Prt*gradT;

dimensionedScalar fG_B_avg1("fG_B_avg1",dimless,gSum(mesh.V()*G.primitiveField())/gSum(mesh.V()));
dimensionedScalar fG_B_avg2("fG_B_avg2",dimless,gSum(mesh.V()*Gb.primitiveField())/gSum(mesh.V()));
dimensionedScalar fG_B_avg4("fG_B_avg4",dimless,gSum(mesh.V()*Cbetas*omega.primitiveField()*k.primitiveField())/gSum(mesh.V()));

Info<< "k production shear : " << fG_B_avg1.value()/fG_B_avg4.value() << endl;
Info<< "k production buoyancy : " << fG_B_avg2.value()/fG_B_avg4.value() << endl;




