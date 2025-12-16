// Info<< "TEqn.H\n" << endl;

// volScalarField alphaEff("alphaEff", nu/Pr + nut/Prt_);
volScalarField alphaEff("alphaEff", nu/Pr + alphaT);


fvScalarMatrix TEqn
(
    fvm::ddt(T)
  + fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
 ==
    radiation->ST(rhoCpRef, T)
  + fvOptions(T)
  + Q_
);

TEqn.relax();

fvOptions.constrain(TEqn);

TEqn.solve();

radiation->correct();

fvOptions.correct(T);

rhok = 1.0 - beta*(T - TRef);
