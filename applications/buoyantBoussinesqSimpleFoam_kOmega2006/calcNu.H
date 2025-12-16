// For convenience, the reference length and temperature difference scales are set to one.
{

dimensionedScalar Ret_avg("Ret_avg",dimless,gSum(mesh.V()*nut.primitiveField())/nu.value()/gSum(mesh.V()));
dimensionedScalar Tavg("Tavg",dimTemperature,gSum(mesh.V()*T.primitiveField())/gSum(mesh.V()));
Nu = (T*U+thf-nu/Pr*gradT);
dimensionedVector Nuavg("Nuavg",dimless,gSum(mesh.V()*Nu.primitiveField())/gSum(mesh.V())/nu.value()*Pr.value());


const scalarField& Tf = T.primitiveField();
scalar Tmax = *std::max_element(Tf.begin(), Tf.end());
const scalarField& nutf = nut.primitiveField();
scalar nutmax = *std::max_element(nutf.begin(), nutf.end());

Info<< "\n" << endl;
Info<< "Ret_avg =" << Ret_avg.value()/Cbetas.value() << endl;
Info<< "Ret_max =" << nutmax/nu.value()/Cbetas.value() << endl;
Info<< "T_avg =" << Tavg.value() << endl;
Info<< "T_max =" << Tmax << endl;
Info<< "Nu_avg =" << Nuavg.value() << "\n" << endl;


// Calculate kinetic boundary-layer thickness and temperature at that location
scalar dnu = 1000.0;
scalar dnu2 = 0.0;
int dnu_indice=-1;
int dnu_indice2=-1;
forAll(T, i)
{
    // scalar dnu0 = nut.internalField()[i]/Prt.value()/(nu.value()/Pr.value());
    scalar dnu0 = nut.internalField()[i]/nu.value();

    if (1.0<dnu0 && T.primitiveField()[i]<Tavg.value()){
        if (dnu0 < dnu){
            dnu_indice=i;
            dnu=dnu0;
        }
    } else if (dnu0<1.0 && T.primitiveField()[i]<Tavg.value()){
        if (dnu2 < dnu0){
            dnu_indice2=i;
            dnu2=dnu0;
        }
    }
}
if (dnu_indice > -1 && dnu_indice2 > -1){
    scalar BLthickness=yWall[dnu_indice] + (scalar(1.0)-dnu)/(dnu2-dnu)*(yWall[dnu_indice2]-yWall[dnu_indice]);
    scalar T_BL=T.primitiveField()[dnu_indice] + (scalar(1.0)-dnu)/(dnu2-dnu)*(T.primitiveField()[dnu_indice2]-T.primitiveField()[dnu_indice]);
    Info<< "Kinetic BL thinkness =" << BLthickness << endl;
    Info<< "T at kinetic BL =" << T_BL  << endl;
}

// Calculate thermal boundary-layer thickness and temperature at that location
dnu = 1000.0;
dnu2 = 0.0;
dnu_indice=-1;
dnu_indice2=-1;
forAll(T, i)
{
    scalar dnu0 = nut.internalField()[i]/(nu.value()/Pr.value())/0.89;
    // scalar dnu0 = nut.internalField()[i]/nu.value();

    if (1.0<dnu0 && T.primitiveField()[i]<Tavg.value()){
        if (dnu0 < dnu){
            dnu_indice=i;
            dnu=dnu0;
        }
    } else if (dnu0<1.0 && T.primitiveField()[i]<Tavg.value()){
        if (dnu2 < dnu0){
            dnu_indice2=i;
            dnu2=dnu0;
        }
    }
}
if (dnu_indice > -1 && dnu_indice2 > -1){
    scalar BLthickness=yWall[dnu_indice] + (scalar(1.0)-dnu)/(dnu2-dnu)*(yWall[dnu_indice2]-yWall[dnu_indice]);
    scalar T_BL=T.primitiveField()[dnu_indice] + (scalar(1.0)-dnu)/(dnu2-dnu)*(T.primitiveField()[dnu_indice2]-T.primitiveField()[dnu_indice]);
    Info<< "Thermal BL thinkness =" << BLthickness << endl;
    Info<< "T at thermal BL =" << T_BL << "\n" << endl;
}


// Calculate Nu on the boundaries

forAll(T.boundaryField(), patchI){

    if (T.boundaryField()[patchI].type() == "fixedValue")
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const labelUList& faceCells = patch.faceCells();

        scalar Tsum1(0.0);
        scalar Fsum(0.0);

        forAll(T.boundaryField()[patchI],faceI)
        {
            Tsum1 = Tsum1 + (T[faceCells[faceI]]-T.boundaryFieldRef()[patchI][faceI])/yWall[faceI]*mesh.boundary()[patchI].magSf()[faceI];
            Fsum = Fsum + mesh.boundary()[patchI].magSf()[faceI];
        }
        Info << "Patch " << mesh.boundary()[patchI].name() << " : Nu = " << Tsum1/Fsum << endl; 
    }
};



forAll(U.boundaryField(), patchI){

    if ((U.boundaryField()[patchI].type() == "noSlip")||(U.boundaryField()[patchI].type() == "fixedValue"))
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        const labelUList& faceCells = patch.faceCells();

        vector Usum1(0.0, 0.0, 0.0);
        scalar Fsum(0.0);

        forAll(U.boundaryField()[patchI],faceI)
        {
            Usum1 = Usum1 + (U[faceCells[faceI]]-U.boundaryFieldRef()[patchI][faceI])/yWall[faceI]*mesh.boundary()[patchI].magSf()[faceI];
            Fsum = Fsum + mesh.boundary()[patchI].magSf()[faceI];
        }

        // scalar Re_tau(mag(Usum1)/Fsum);
        // scalar Re_tau(sqrt(mag(Usum1)/Fsum));
        Info << "Patch " << mesh.boundary()[patchI].name() << " : Re_tau = " << std::sqrt(mag(Usum1)/Fsum/nu.value()) << endl; 
    }
};





dimensionedScalar fG_B_avg1("fG_B_avg1",dimless,gSum(mesh.V()*G.primitiveField())/gSum(mesh.V()));
dimensionedScalar fG_B_avg2("fG_B_avg2",dimless,gSum(mesh.V()*Gb.primitiveField())/gSum(mesh.V()));
dimensionedScalar fG_B_avg4("fG_B_avg4",dimless,gSum(mesh.V()*Cbetas*omega.primitiveField()*k.primitiveField())/gSum(mesh.V()));

Info<< "k production shear : " << fG_B_avg1.value()/fG_B_avg4.value() << endl;
Info<< "k production buoyancy : " << fG_B_avg2.value()/fG_B_avg4.value() << endl;

volScalarField Usqr(U&U);
dimensionedScalar Usqr_avg(gSum(mesh.V()*Usqr.primitiveField())/gSum(mesh.V()));
dimensionedScalar Rey(sqrt(Usqr_avg));

dimensionedVector U_avg(gSum(mesh.V()*U.primitiveField())/gSum(mesh.V()));

Info<< "Re_avg =" << U_avg.value()/nu.value()  << endl;
Info<< "Re_max =" << Rey.value()/nu.value() << "\n" << endl;


}


