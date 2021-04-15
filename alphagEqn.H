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

// Want to find the index opsition where alpha1=0.9, starting from the wall


label mypatch(1000); 

// Info << "mesh.boundaryMesh().patchID()=" << mesh.boundaryMesh().patchID()[2] << endl;
// Info << "alpha1.boundaryField().size()=" << alpha1.boundaryField().size() << endl;
// Info << "alpha1.boundaryField()[0].size()=" << alpha1.boundaryField()[0].size() << endl;
// Info << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

forAll(mesh.boundary(), patchi)
{
    // Info << "mesh.boundary()[patchi].type()=" << mesh.boundary()[patchi].type() << endl;
    // Info << "mesh.boundary()[patchi].name()=" << mesh.boundary()[patchi].name() << endl;
    if(mesh.boundary()[patchi].type()=="wall"){
	mypatch = patchi;
    	Info << "Found our patch: " << mesh.boundary()[patchi].type() << ", " <<
	    mesh.boundary()[patchi].name() << ", index=" << mypatch << endl;
    }
}

const fvPatch & patch(mesh.boundary()[mypatch]);


const labelList & owner = mesh.faceOwner();
const labelList & neighbour = mesh.faceNeighbour();
const vectorField & faceCentres = mesh.Cf().primitiveField();
const List<cell> & cells = mesh.cells();

const vectorField & patchFaceCentres = patch.Cf();
const tmp<vectorField> tfaceNormals = patch.nf();
const vectorField faceNormals = tfaceNormals();
const tmp<vectorField> tcellCentres = patch.Cn();
const vectorField cellCentres = tcellCentres();
const volVectorField & C = mesh.C();
const UList<label> & faceCells = patch.faceCells();
const List<face> & faces = mesh.faces();
const label patchStart = patch.start();
const polyBoundaryMesh & boundaryMesh = mesh.boundaryMesh();

//labelList indexList(patch.size());
DynamicList<label> indexList;

//indexList[0]=1;

// Ok, then I have one single patch which is the wall patch, which I want to go
// through - goes through all these patches to find the h_90 surface on each of
// them

// Info << "patch.size()=" << patch.size() << " " << std::round(patch.size()) << endl;


// Info << "alpha1.boundaryField()[mypatch][0] " << alpha1.boundaryField()[mypatch][0] << endl; 

label nLayers = 200;


// Info << "faceCells[0]=" << faceCells[0] << " " << faces[4]  << endl;
// Info << "mesh.nInternalFaces() " << mesh.nInternalFaces() << endl;

label counter = 1;

forAll(patch, patchFaceI)
{
    // Info << "patchFaceI=" << patchFaceI << endl;

    point patchFaceCenterI = patchFaceCentres[patchFaceI];
    label startCellIndex = faceCells[patchFaceI];
    cell startCell = cells[startCellIndex];
    label startFaceLabel = patchStart + patchFaceI;
    label nextCellIndex = -1;

    // Info << "alpha1 " << startCellIndex << endl;
    // Info << "alpha1 " << alpha1[startCellIndex] << endl;

    // Monday: Start step by step on my own, perhaps start by one layer and
    // print out to get control, should also start easy and make my own code
    // based on Tims code, and not copying the whole thing
    for (label layer=0; layer < nLayers-1; layer++)
    {
	label opposingFace = startCell.opposingFaceLabel(startFaceLabel, faces);
	scalar distance =
                mag(faceCentres[opposingFace] - patchFaceCentres[patchFaceI]);


	if (owner[opposingFace] == startCellIndex)
        {
	    nextCellIndex = neighbour[opposingFace];
	}
	else
        {
	    nextCellIndex = owner[opposingFace];
	    // Info << "heihei" << endl;
	}
        
	    
	// // Info << "distance = " << opposingFace  << endl;
	// Info << "startCellIndex = " << owner[startFaceLabel]  << endl;
	// Info << "nextCellIndex = " << nextCellIndex  << endl;
	// Info << "alpha1 " << alpha1[nextCellIndex] << " layer " << layer << endl;

	if(alpha1[nextCellIndex] <= 0.9)
        {
	    // indexList[patchFaceI] = nextCellIndex;
	    indexList.append(nextCellIndex);
	    // Info << nextCellIndex << " Distance " << distance << " layer " << layer << endl; 
	    // Info << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	    // Info << indexList.size() << " " << patch.size()  << endl;
	    break;
	}


	// Info << "layer " << layer << endl;
	if(layer==nLayers-2)
	{
	    Info << "????????????????????" << " - find out what is wrong!" << endl;
	}
	
	startFaceLabel = opposingFace;
	startCellIndex = nextCellIndex;
	startCell = cells[startCellIndex];
	
    }
    // Info << "indexList " << indexList << endl;
    counter++;
    // if(counter>220)
    // 	break;

}

// // Info << indexList << endl;
// forAll(indexList, indx)
// {
//     // // Info << "indx " << indx << endl;
//     // Info << "alpha1[" << indexList[indx] << "] =" << alpha1[indexList[indx]] <<
//     // 	" deltafs: " << deltaFreeSurface[indexList[indx]] << " indexList[indx] "
//     // 	 << indexList[indx] << endl;
//     deltaFreeSurface[indexList[indx]] = 1;
// }

// Info << deltaFreeSurface.size() << endl;


forAll(deltaFreeSurface, i)
{
    bool flag = 0;
// Info << indexList << endl;
    forAll(indexList, indx)
    {
	if(i==indexList[indx])
	    flag = 1;
    }

    if(flag)
	deltaFreeSurface[i] = 1;
    else
	deltaFreeSurface[i] = 0;
}


// deltaFreeSurface = 
//     posPart(-0.25/aParabola*sqr(mag(gradAlpha1) - gradAlpha1Cr) + 1);



// // List of cells which are definitely close to the free surface
// DynamicList<label> interfaceCells(1000);

// // Set the function to 1 if above gradAlphaCr
// forAll(deltaFreeSurface, cellI)
// {
//     if (mag(gradAlpha1[cellI]) > gradAlpha1Cr.value())
//     {
//         deltaFreeSurface[cellI] = 1;
//     }
//     if (deltaFreeSurface[cellI] > 0.9)
//     {
//         interfaceCells.append(cellI);
//     }
// }

// interfaceCells.shrink();
// Info << interfaceCells.size() << nl;
// const volVectorField & C = mesh.C();

// forAll(deltaFreeSurface, i)
// {

//     if ((deltaFreeSurface[i] < SMALL) || (deltaFreeSurface[i] > 0.9))
//     {
//         continue;
//     }
//     else
//     {
//         scalar dist = 0;
//         bool flag = 1;

//         forAll (interfaceCells, cellI)
//         {

//             const vector vdist = C[interfaceCells[cellI]] - C[i];
//             dist = Foam::sqrt(vdist & vdist);

//             if (dist < interfaceThickness.value())
//             {
//                 flag = 0;
//                 break;
//             }
//         }
//         // No cell marked as interface found in 4*average(cellSize) proximity
//         // Set deltaFunction to 0 for the cell
//         if (flag)
//         {
//             deltaFreeSurface[i] = 0;
//         }
//     }
// }

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
