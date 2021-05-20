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
DynamicList<scalar> fractionDeltafs;

//indexList[0]=1;


label nLayers = 200;

label counter = 1;


forAll(patch, patchFaceI)
{
    // Info << "patchFaceI=" << faceNormals[patchFaceI] << endl;
    // Info << "g=" << g.value() << endl;

    // Only crawl in vertical direction
    label crawlingDir = -1;
    forAll(g.value(),i)
    {
	// Info << "i=" << g.value()[i] << endl;
	if(g.value()[i] < 0)
	    crawlingDir = i;
    }

    if(faceNormals[patchFaceI][crawlingDir] >= 0)
    {
	// Info << "Do only crawl in vertical direction " <<
	//     faceNormals[patchFaceI] << endl;
	continue;
    }

    point patchFaceCenterI = patchFaceCentres[patchFaceI];
    label startCellIndex = faceCells[patchFaceI];
    cell startCell = cells[startCellIndex];
    label startFaceLabel = patchStart + patchFaceI;
    label nextCellIndex = -1;
    label fsCellIndex = -1;
    
    // Info << "alpha1 " << startCellIndex << endl;
    // Info << "alpha1 " << alpha1[startCellIndex] << endl;

    bool flag = 0;
    

    scalar distanceFS = -1;
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
        
	    
	if(alpha1[nextCellIndex] <= 0.0001)
        {
	    // Here one should go down again
	    // indexList[patchFaceI] = nextCellIndex;
	    fsCellIndex = nextCellIndex;
	    distanceFS = distance;
	    // Info << nextCellIndex << " Distance " << distance << " layer " << layer << endl; 
	    // Info << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	    // Info << indexList.size() << " " << patch.size()  << endl;
	    flag = 1;
	    // break;
	}

	// Info << "layer " << layer << endl;
	if(layer==nLayers-2)
	{
	    Info << "????????????????????" << " - find out what is wrong!" << endl;
	}

	if(flag){
	    // indexList.append(fsCellIndex);
	    startFaceLabel = opposingFace;
	    //startCellIndex = nextCellIndex;
	    // startCell = cells[startCellIndex];
	    break;
	}
	
	startFaceLabel = opposingFace;
	startCellIndex = nextCellIndex;
	startCell = cells[startCellIndex];

    }

    // Starting a new crawling loop, starting from h100, and crawling down until
    // we reach h90 - so make that as the first next step
    // scalar distanceFS = -1;

    bool h90flag = 0;
    label startCellh100 = startFaceLabel; //startCellIndex;
    
	
    for (label layer=0; layer < nLayers-1; layer++)
    {
    	label opposingFace = startCell.opposingFaceLabel(startFaceLabel, faces);
    	scalar distance =
                mag(faceCentres[startCellh100] - faceCentres[opposingFace]);

    	if (owner[opposingFace] == startCellIndex)
        {
    	    nextCellIndex = neighbour[opposingFace];
    	}
    	else
        {
    	    nextCellIndex = owner[opposingFace];
    	    // Info << "heihei" << endl;
    	}

    	if(alpha1[nextCellIndex] >= 0.1) // find h90
        {
	    if(h90flag)
	    {
		// Info << "heiheiheihei" << endl;
		scalar dfs = -1/interfaceThickness.value()*distance + 1; // interfaceThickness
		// Info << "dfs: " << dfs  << " distance " << distance <<
		//     " 1/phiEnt.value()*distance " <<  distance/phiEnt.value() <<  endl;
		if(dfs > 1e-5){ // ): distance < phiEnt.value() 
		    fractionDeltafs.append(dfs);
		    indexList.append(nextCellIndex);
		}else
		    break;
		// Info << "distance from h100: " << distance << " " << faceCentres[startCellh100] <<  " " << faceCentres[opposingFace] << endl; 
	    }else
	    {
		// Info << "ID: " << nextCellIndex  << endl;
		// Info << "distance from h100: " << distance << " " << faceCentres[startCellh100] <<  " " << faceCentres[opposingFace] << endl; 
		//Info << "heiheiheihei" << endl;
		fsCellIndex = nextCellIndex;
		indexList.append(fsCellIndex);
		fractionDeltafs.append(1);
		h90flag = 1;
	    }
    	    //distanceFS = distance;
    	    // Info << nextCellIndex << " Distance " << distance << " layer " << layer << endl; 
    	    // Info << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    	    // Info << indexList.size() << " " << patch.size()  << endl;
    	    //break;
	    // if(distance > phiEnt.value()){ //0.01){ //alpha1[nextCellIndex] > 0.5)
	    // 	// Info << "heiheiheihei" << endl;
	    // 	break;
	    // }
    	}
	startFaceLabel = opposingFace;
	startCellIndex = nextCellIndex;
	startCell = cells[startCellIndex];

    }
    
    // Info << "indexList " << indexList << endl;
    counter++;
    // if(counter>10)
    // 	break;

}


forAll(deltaFreeSurface, i)
{
    bool flag = 0;
// Info << indexList << endl;
    scalar myIndx = -1;
    forAll(indexList, indx)
    {
	if(i==indexList[indx]){
	    flag = 1;
	    myIndx = indx;
	}
    }

    if(flag){
	deltaFreeSurface[i] = 1*fractionDeltafs[myIndx];
    }
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

forAll(dUdn,i)
{
    if(mesh.C()[i][0] < 0.7){
	// Info << "cellCentres[i]" << mesh.C()[i] << endl;
	deltaFreeSurface[i] = 0;
    }
}


// Compute the surface roughness
a = 0.5*Foam::pow(0.09, 0.75)
                   *Foam::pow(turbulence->k(), 1.5)
                   /turbulence->epsilon();

// Compute the source term for the alphag equation
alphagSource = a / phiEnt * dUdn * deltaFreeSurface;

scalar dbba = 0.05;

// Bound the source term in space using appropriate criteria
alphagSource *= pos(1.0 - alphag);
alphagSource *= pos(alpha1 - dbba);
alphagSource *= pos(UnHat);


// Info << "kCritical" << kCritical << endl;
// if (kCritical.value() > 0)
// {
//     alphagSource *= pos(turbulence->k() - kCritical);
// }
// else
// {
//     Pt = turbulence->k()*rho;
//     Pd = rho*mag(g & nHat)*a + sigma/a;  
//     alphagSource *= pos(Pt - Pd);
// }

Info<< "Air entrainment source: "
    << "  Min(alphagSource) = " << min(alphagSource).value()
    << "  Max(alphagSource) = " << max(alphagSource).value()
    << endl;


/*              ASSEMBLE AND SOLVE THE ALPHAG EQUATION                        */

scalarField & penaltySourceField = penaltySource.primitiveFieldRef();
penaltySourceField = -neg(alpha1 - dbba)/runTime.deltaTValue();


// *** Solve momentum equation for Ug
if (momentumType == "slipModel")
{
    #include "slipVelocity.H"
}
else if (momentumType == "off")
{
    Ug = U;
}

// *** Print the velocity values
Info << "Velocities: "
<< "Min(U) = " << min(mag(U)).value()
<< " Max(U) = " << max(mag(U)).value()
<< " Min(Ug) = " << min(mag(Ug)).value()
<< " Max(Ug) = " << max(mag(Ug)).value()
<< endl;

phiUg = linearInterpolate(Ug) & mesh.Sf();

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
