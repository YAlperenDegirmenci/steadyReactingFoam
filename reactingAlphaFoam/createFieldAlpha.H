// ───────── Alpha field (prescribed liquid volume-fraction) ──────────────
Info<< "Reading field alphaEl\n" << endl;

volScalarField alphaEl
(
    IOobject
    (
        "alphaEl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE        // you may switch to NO_WRITE
    ),
    mesh
);

// face-interpolated liquid fraction
surfaceScalarField alphaElF
(
    IOobject("alphaElF", runTime.timeName(), mesh),
    fvc::interpolate(alphaEl)
);

// liquid-phase flux   φ_el = α_elF * φ
surfaceScalarField phiEl
(
    IOobject("phiEl", runTime.timeName(), mesh),
    alphaElF * phi
);
