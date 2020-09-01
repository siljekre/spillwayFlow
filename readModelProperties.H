Info<< "Reading airModelProperties\n" << endl;

IOdictionary airModelProperties
(
    IOobject
    (
        "airModelProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar deltaFunctionCutoff
(
    "deltaFunctionCutoff",
     dimless/dimLength,
     airModelProperties
);

dimensionedScalar phiEnt("phiEnt", dimLength, airModelProperties);

dimensionedScalar Ct
(
    "diffusionConst",
    dimless,
    airModelProperties.lookupOrDefault("diffusionConst", 0)
);

dimensionedScalar interfaceThickness
(
    "interfaceThickness",
    dimLength,
    airModelProperties
);

dimensionedScalar kCritical("kCritical", dimEnergy/dimMass, airModelProperties);

dimensionedScalar diffusionFlag
(
    "diffusionFlag",
    dimless,
    airModelProperties.lookupOrDefault("diffusiveTerm", false)
);

dimensionedScalar bba
(
    "BBA",
    dimless,
    airModelProperties.lookupOrDefault("BBA", false)
);

Info<< "Air entrainment model parameters" << endl;
Info<< "    phiEnt\t\t\t" << phiEnt.value() << endl;
Info<< "    Diffusion term active?\t" << diffusionFlag.value() << endl;
Info<< "    Diffusion constant value:\t" << Ct.value() << endl;
Info<< "    BBA active?\t\t\t" << bba.value() << endl;
Info<< "    interfaceThickness\t\t" << interfaceThickness.value() << endl;
if (kCritical.value() <= 0)
{
    Info<< "    Critical k value not positive, will use the criterion based on"
        << " energy balance." << endl;
}
else
{
    Info<< "    kCritical\t\t\t" << kCritical.value() << endl;
}
