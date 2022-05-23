    // Momentum predictor

    MRF.correctBoundaryVelocity(U);
    
    dimensioned<scalar> tend = runTime.endTime().value();
    dimensioned<scalar> tnow = atoi(runTime.timeName().c_str());

    dimensionedScalar xi = ((0.8 - 0.3) / (tend - 20000)) * (tnow - 20000) + 0.3;
    
    if (tnow > tend - 30000)
    {
    	xi = 0.8;
    }
    else
    {
    	xi = xi;
    }
    
    Info << "Momentum equation xi = " << xi << nl << endl;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + (1 - xi) * MRF.DDt(U)
      + (1 - xi) * turbulence->divDevReff(U)
      + xi * fvc::div(R)
      + xi * fvc::laplacian(turbulence->nut(), U, "laplacian(nuEff, U)")
      - xi * fvm::laplacian(nu+turbulence->nut(), U)
      - xi * fvc::div(nu*dev2(T(fvc::grad(U))))
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
