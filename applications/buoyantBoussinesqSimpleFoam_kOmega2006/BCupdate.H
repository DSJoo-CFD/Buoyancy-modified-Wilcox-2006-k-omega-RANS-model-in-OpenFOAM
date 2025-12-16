
// omega boundary condition

forAll(omega.boundaryField(), patchI){

    if (omega.boundaryField()[patchI].type() == "fixedValue")
    {
        // const fvPatch& patch = mesh.boundary()[patchI];
        // const labelUList& faceCells = patch.faceCells();

        forAll(omega.boundaryField()[patchI],faceI)
        {
            omega.boundaryFieldRef()[patchI][faceI] = 6.0*nu.value()/Cbeta0.value()/sqr(yWall[faceI]); // Wilcox (2006)
            // omega.boundaryFieldRef()[patchI][faceI] = 60.0/0.075*nu.value()/sqr(yWall[faceI]); // Mentor (1994)
        }

        omega.boundaryFieldRef()[patchI].updateCoeffs();
        // Info << "Patch " << mesh.boundary()[patchI].name() << " : omega BC updated " << endl;
    }
};

