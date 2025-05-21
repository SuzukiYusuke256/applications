#include "fvCFD.H"

int main(int argc, char *argv[])
{
	// Read
	// omajinai
	Foam::argList args(argc,argv);
	if(not args.checkRootCase()) Foam::FatalError.exit();

	Foam::Time runTime(Foam::Time::controlDictName, args);

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

	// configuration
	IOdictionary mmdDict
    (
        IOobject
        (
            "myManualDecomposeDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

	// read the dictionary entries
    const vector baseRegionMin(mmdDict.get<vector>("baseRegionMin")); // Minimum corner of baseRegion
	const vector baseRegionMax(mmdDict.get<vector>("baseRegionMax"));  // Maximum corner of baseRegion
	const vector fineRegionMin(mmdDict.get<vector>("fineRegionMin"));  // Minimum corner of fineRegion
	const vector fineRegionMax(mmdDict.get<vector>("fineRegionMax"));  // Maximum corner of fineRegion

	const vector Nb(mmdDict.get<vector>("baseRegionDivision"));  // division number of baseRegion
	const vector Nf(mmdDict.get<vector>("fineRegionDivision"));  // division number of fineRegion

	// Write
	labelList procIds(mesh.cells().size());

	// Define parameters for baseRegion (outer) and fineRegion (inner)
	// const vector baseRegionMin(0.0, 0.0, 0.0);  // Minimum corner of baseRegion
	// const vector baseRegionMax(1.0, 1.0, 1.0);  // Maximum corner of baseRegion
	// const vector fineRegionMin(0.25, 0.25, 0.25);  // Minimum corner of fineRegion
	// const vector fineRegionMax(0.75, 0.75, 0.75);  // Maximum corner of fineRegion

	// Number of divisions for each region in each direction
	// const label Nbx = 2;  // Divisions of baseRegion in x direction
	// const label Nby = 2;  // Divisions of baseRegion in y direction
	// const label Nbz = 2;  // Divisions of baseRegion in z direction
	// const label Nfx = 2;  // Divisions of fineRegion in x direction
	// const label Nfy = 2;  // Divisions of fineRegion in y direction
	// const label Nfz = 2;  // Divisions of fineRegion in z direction

	// Process each cell
	forAll(mesh.cells(), cid) {
		const vector &C = mesh.C()[cid];
		
		// Debug output
		Info << "num : " << cid << " x : " << C[0] << " y : " << C[1] << " z : " << C[2] << endl;
		
		// Check if point is inside fineRegion
		bool insideFineRegion = 
			C[0] >= fineRegionMin[0] && C[0] <= fineRegionMax[0] &&
			C[1] >= fineRegionMin[1] && C[1] <= fineRegionMax[1] &&
			C[2] >= fineRegionMin[2] && C[2] <= fineRegionMax[2];
		
		if (insideFineRegion) {
			// Point is inside fineRegion
			// Directly calculate the region ID based on the cell position
			scalar normalizedX = (C[0] - fineRegionMin[0]) / (fineRegionMax[0] - fineRegionMin[0]);
			scalar normalizedY = (C[1] - fineRegionMin[1]) / (fineRegionMax[1] - fineRegionMin[1]);
			scalar normalizedZ = (C[2] - fineRegionMin[2]) / (fineRegionMax[2] - fineRegionMin[2]);
			
			// Calculate the indices directly
			label ix = static_cast<label>(normalizedX * Nf[0]);
			label iy = static_cast<label>(normalizedY * Nf[1]);
			label iz = static_cast<label>(normalizedZ * Nf[2]);
			
			// Ensure indices are within valid range (alternative to std::min)
			if (ix == Nf[0]) ix = Nf[0] - 1;
			if (iy == Nf[1]) iy = Nf[1] - 1;
			if (iz == Nf[2]) iz = Nf[2] - 1;
			
			// Compute linear index for the fineRegion
			procIds[cid] = ix + iy * Nf[0] + iz * Nf[0] * Nf[2];
		}
		else {
			// Point is inside baseRegion but outside fineRegion
			scalar normalizedX = (C[0] - baseRegionMin[0]) / (baseRegionMax[0] - baseRegionMin[0]);
			scalar normalizedY = (C[1] - baseRegionMin[1]) / (baseRegionMax[1] - baseRegionMin[1]);
			scalar normalizedZ = (C[2] - baseRegionMin[2]) / (baseRegionMax[2] - baseRegionMin[2]);
			
			// Calculate the indices directly
			label ix = static_cast<label>(normalizedX * Nb[0]);
			label iy = static_cast<label>(normalizedY * Nb[1]);
			label iz = static_cast<label>(normalizedZ * Nb[2]);
			
			// Ensure indices are within valid range
			if (ix == Nb[0]) ix = Nb[0] - 1;
			if (iy == Nb[1]) iy = Nb[1] - 1;
			if (iz == Nb[2]) iz = Nb[2] - 1;
			
			// Compute linear index for the baseRegion with offset
			procIds[cid] = ix + iy * Nb[0] + iz * Nb[0] * Nb[0] + Nf[0] * Nf[1] * Nf[2];
		}
	}

	labelIOList cellDecomposition
	(
		IOobject
		(
			"cellDecomposition",
			runTime.constant(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE,
			IOobject::NO_REGISTER
		),
		procIds
	);

	cellDecomposition.write();

	return 0;
}
