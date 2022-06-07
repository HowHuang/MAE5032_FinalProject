#ifndef _EXPLICIT_H_
#define _EXPLICIT_H_

#include <petsc.h>
#include <hdf5.h>

typedef struct InputPara
{
    PetscScalar dt;     PetscScalar dl;     PetscScalar rho;
    PetscScalar c;      PetscScalar k;      PetscScalar f;
}InputPara;

typedef struct Bound
{
    PetscScalar ub;     PetscScalar hb;     PetscScalar tb;
    PetscScalar ut;     PetscScalar ht;     PetscScalar tt;
    PetscScalar ul;     PetscScalar hl;     PetscScalar tl;
    PetscScalar ur;     PetscScalar hr;     PetscScalar tr;
}Bound; // ~ for a point

typedef struct IterMaterial
{
    PetscScalar W;      PetscScalar E;
    PetscScalar N;      PetscScalar S;
    PetscScalar P;      PetscScalar b;
}IterMaterial;

enum Location
{
    RightTop    = 1,    TopSide     = 2,
    LeftTop     = 3,    LeftSide    = 4,
    LeftBottom  = 5,    BottomSide  = 6,
    RightBottom = 7,    RightSide   = 8,
    Internal    = 9
};

void CalIterationMaterial(InputPara* IP, Bound* bound, IterMaterial* IM, enum Location loc);

#endif
