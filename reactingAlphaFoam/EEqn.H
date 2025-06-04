{
    volScalarField& he = thermo.he();

/*
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        reaction->Qdot()
      + fvModels.source(rho, he)
    );
*/

	fvScalarMatrix EEqn
	(
		fvm::ddt(alphaEl*rho, he)
	  + mvConvection->fvmDiv(phiEl, he)
	  + fvc::ddt(alphaEl*rho, K)
	  + fvc::div(phiEl, K)
	  + (
			he.name() == "e"
		  ? mvConvection->fvcDiv(fvc::absolute(phiEl, rho, U), p/rho)
		  : -alphaEl*dpdt
		)
	  + alphaEl*thermophysicalTransport->divq(he)
	 ==
		alphaEl*reaction->Qdot()
	  + fvModels.source(alphaEl*rho, he)
	);



    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
