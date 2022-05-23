    // Momentum predictor

    MRF.correctBoundaryVelocity(U);
    
    dimensioned<scalar> tend = runTime.endTime().value();
    dimensioned<scalar> tnow = atoi(runTime.timeName().c_str());

    dimensionedScalar xi = ((0.8 - 0) / (tend - 20000 - 10000)) * (tnow - 20000) + 0;
    
    if (tnow > tend - 20000)
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
      + xi * fvc::div( ((2.0/3.0)*turbulence->k()*I + aij), "div(devRhoReff)")
      + xi * fvc::laplacian(turbulence->nut(), U, "laplacian(nuEff,U)")
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
