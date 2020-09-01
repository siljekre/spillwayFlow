/*              COMPUTE THE FREE SURFACE DETECTION FUNCTION                   */

label n = 4;

// Recompute the gradient of alpha1
gradAlpha1 = fvc::grad(alpha1);

// Compute the critical value of grad(alpha1)
dimensionedScalar gradAlpha1Cr("gradAlpha1Cr", 1/(n*averageCellSize));

Info<< "Critical grad(alpha1) value used in free surface detection = " 
    << gradAlpha1Cr.value() << endl; 

// Distance from the focus point of the parabola to its top point
dimensionedScalar aParabola
(
    "aParabola",
     0.25*sqr(gradAlpha1Cr - deltaFunctionCutoff)
);

deltaFreeSurface = 
    posPart(-0.25/aParabola*sqr(mag(gradAlpha1) - gradAlpha1Cr) + 1);

// List of cells which are definitely close to the free surface
DynamicList<label> interfaceCells(1000);

// Set the function to 1 if above gradAlphaCr
forAll(deltaFreeSurface, cellI)
{
    if (mag(gradAlpha1[cellI]) > gradAlpha1Cr.value())
    {
        deltaFreeSurface[cellI] = 1;
    }
    if (deltaFreeSurface[cellI] > 0.9)
    {
        interfaceCells.append(cellI);
    }
}

interfaceCells.shrink();
Info << interfaceCells.size() << nl;
const volVectorField & C = mesh.C();

forAll(deltaFreeSurface, i)
{

    if ((deltaFreeSurface[i] < SMALL) || (deltaFreeSurface[i] > 0.9))
    {
        continue;
    }
    else
    {
        scalar dist = 0;
        bool flag = 1;

        forAll (interfaceCells, cellI)
        {

            const vector vdist = C[interfaceCells[cellI]] - C[i];
            dist = Foam::sqrt(vdist & vdist);

            if (dist < interfaceThickness.value())
            {
                flag = 0;
                break;
            }
        }
        // No cell marked as interface found in 4*average(cellSize) proximity
        // Set deltaFunction to 0 for the cell
        if (flag)
        {
            deltaFreeSurface[i] = 0;
        }
    }
}

/*           COMPUTE THE SOURCE TERM IN THE ALPHAG EQUATION                   */

// Interface-normal unit vectors
volVectorField nHat(gradAlpha1/(mag(gradAlpha1) + mixture.deltaN()));

// Interface-normal velocity derivative, keep positive values only
UnHat = U & nHat;
dUdn = Foam::posPart(fvc::grad(UnHat) & nHat);

// Compute the surface roughness
volScalarField a = 0.5*Foam::pow(0.09, 0.75)
                   *Foam::pow(turbulence->k(), 1.5)
                   /turbulence->epsilon();

// Compute the source term for the alphag equation
alphagSource = a / phiEnt * dUdn * deltaFreeSurface;

// Bound the source term in space using appropriate criteria
alphagSource *= pos(1.0 - alphag);
alphagSource *= pos(alpha1 - 0.1);
alphagSource *= pos(UnHat);

if (kCritical.value() > 0)
{
    alphagSource *= pos(turbulence->k() - kCritical);
}
else
{
    Pt = turbulence->k()*rho;
    Pd = rho*mag(g & nHat)*a + sigma/a;  
    alphagSource *= pos(Pt - Pd);
}

Info<< "Air entrainment source: "
    << "  Min(alphagSource) = " << min(alphagSource).value()
    << "  Max(alphagSource) = " << max(alphagSource).value()
    << endl;


/*              ASSEMBLE AND SOLVE THE ALPHAG EQUATION                        */

scalarField & penaltySourceField = penaltySource.primitiveFieldRef();
penaltySourceField = -neg(alpha1 - 0.1)/runTime.deltaTValue();

fvScalarMatrix alphagEqn
(
    fvm::ddt(alphag)
  + fvm::div(phi, alphag)
  - diffusionFlag*Ct*fvm::laplacian(mixture.nuf(), alphag)
        ==
    alphagSource
  + bba*fvm::Sp(penaltySource, alphag)
);

// Only introduce air below h90
alphagEqn.solve();

alphag.correctBoundaryConditions();

Info<< "Phase-g volume fraction = "
    << alphag.weightedAverage(mesh.Vsc()).value()
    << "  Min(alphag) = " << min(alphag).value()
    << "  Max(alphag) = " << max(alphag).value()
    << endl;
