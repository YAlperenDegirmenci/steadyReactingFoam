tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
		mesh.divScheme("div(phiEl,Yi_h)")
    )
);

reaction->correct();

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

/*
        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
          + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport->divj(Yi)
         ==
            reaction->R(Yi)
          + fvModels.source(rho, Yi)
        );
*/

		fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi)
			+ mvConvection->fvmDiv(phiEl, Yi)
			+ thermophysicalTransport->divj(Yi)

         ==
            alphaEl*reaction->R(Yi)
		+ fvModels.source(alphaEl*rho, Yi)
        );
		

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
