Info<< "Reading alphag" << endl;
volScalarField alphag
(
    IOobject
    (
        "alphag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alphagSource
(
    IOobject
    (
        "alphagSource",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphagSource", dimless/dimTime, 0)
);

volScalarField penaltySource
(
    IOobject
    (
        "penaltySource",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("penaltySource", dimless/dimTime, 0)
);

volScalarField UnHat
(
    IOobject
    (
        "UnHat",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("UnHat", dimVelocity, 0)
);

volScalarField a
(
    IOobject
    (
        "a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("a", dimLength, 0)
);

volVectorField gradAlpha1
(
    IOobject
    (
        "gradAlpha1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradAlpha1", dimless/dimLength, pTraits<vector>::zero)
);

volScalarField dUdn
(
IOobject
   (
        "dUdn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
   ),
   mesh,
   dimensionedScalar("dUdn",dimVelocity/dimLength, 0)
);

volScalarField deltaFreeSurface
(
    IOobject
    (
        "deltaFreeSurface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("deltaFreeSurface", dimless, 0)
);

volScalarField Pd
(
    IOobject
    (
        "Pd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,    
    dimensionedScalar("Pd", dimEnergy/dimVolume, 0)
);

volScalarField Pt
(
    IOobject
    (
        "Pt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,    
    dimensionedScalar("Pt", dimEnergy/dimVolume, 0)
);

Info<< "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar sigma
(
    "sigma",
    dimensionSet(1, 0, -2, 0, 0),
    transportProperties
);

volVectorField Ug
(
    IOobject
    (
        "Ug",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    U
);
Ug.write();

volVectorField Ur
(
    IOobject
    (
        "Ur",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (Ug - U)
);
Ur.write();

// *** Variables stored for ddt and div

surfaceScalarField phiUg
(
    IOobject
    (
        "phiUg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(Ug) & mesh.Sf()
);

volScalarField ws
(
IOobject
   (
        "ws",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
   ),
   mesh,
   dimensionedScalar("ws",dimensionSet(0,0,0,0,0,0,0),1.0)
);


// Estimate of the mesh length scale
DimensionedField<scalar, volMesh> cellSize 
(
    IOobject
    (
        "cellSize",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimLength,
    mesh.V()
);

// Get the dimensionless field
scalarField & cellSizeValues =  cellSize.field();

// Handle 2d cases specially since the volume of the cell is arbitrary in 2d.
// Loop through all mesh directions and find the one invalid.
// Divide the cellSize by the span in that direction.
if (mesh.nGeometricD() == 3)
{
    cellSizeValues = cbrt(cellSizeValues);
}
else
{
    const vector meshSpan = mesh.bounds().span();
    const vector geometricD = mesh.geometricD();
    forAll (geometricD, i)
    {
        if (geometricD[i] == -1)
        {
            cellSizeValues /= meshSpan[i];
        }
    }
    cellSizeValues = sqrt(cellSizeValues);
}

dimensionedScalar averageCellSize("averageCellSize", cellSize.average());




