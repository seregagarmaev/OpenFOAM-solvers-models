    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + fvc::div( ((2.0/3.0)*turbulence->k()*I + aij), "div(devRhoReff)")
      + fvc::laplacian(turbulence->nut(), U, "laplacian(nuEff,U)")
      - fvm::laplacian(nu+turbulence->nut(), U)
      - fvc::div(nu*dev2(T(fvc::grad(U))))
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
