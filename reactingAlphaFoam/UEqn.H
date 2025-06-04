// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

/*
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
 ==
    fvModels.source(rho, U)
);
*/

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(alphaEl*rho, U)                 // storage term
  + fvm::div(phiEl,       U)                 // convection uses liquid flux
  + MRF.DDt(alphaEl*rho,  U)                 // if you keep MRF
  //+ turbulence->divDevTau(U)*alphaEl         // viscosity only in liquid
  + alphaEl * turbulence->divDevTau(U)

 ==
    fvModels.source(alphaEl*rho, U)          // any sources scaled
);


fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
