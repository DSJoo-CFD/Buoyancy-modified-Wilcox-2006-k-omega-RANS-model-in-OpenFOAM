
// Wall distance
wallDist y(mesh);
const volScalarField& yWall = y.y();
volScalarField WallDist_y
(
    IOobject
    (
        "WallDist_y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    yWall
);
