// 

#ifndef _PLICIT_H_
#define _PLICIT_H_

#include <petsc.h>
#include <hdf5.h>

// 用来存储全局物理参数的结构体。
typedef struct InputPara
{
    PetscScalar dt;     PetscScalar dl;     PetscScalar rho;
    PetscScalar c;      PetscScalar k;      PetscScalar f;
}InputPara;

// 用来存储某一边界点的边界条件，u代表问题，h代表热流，t代表边界条件类型
typedef struct Bound
{
    PetscScalar ub;     PetscScalar hb;     PetscScalar tb;
    PetscScalar ut;     PetscScalar ht;     PetscScalar tt;
    PetscScalar ul;     PetscScalar hl;     PetscScalar tl;
    PetscScalar ur;     PetscScalar hr;     PetscScalar tr;
}Bound; // ~ for a point

// 用来存储系数矩阵和右手向量b中的元素值
typedef struct IterMaterial
{
    PetscScalar W;      PetscScalar E;
    PetscScalar N;      PetscScalar S;
    PetscScalar P;      PetscScalar b;
}IterMaterial;

// 用来标记该点的位置
typedef enum Location
{
    RightTop    = 1,    TopSide     = 2,
    LeftTop     = 3,    LeftSide    = 4,
    LeftBottom  = 5,    BottomSide  = 6,
    RightBottom = 7,    RightSide   = 8,
    Internal    = 9
}Location;

// 计算稀疏矩阵和右手向量元素值的函数
void ExplicitIterationMaterial(InputPara* IP, Bound* bound, IterMaterial* IM, enum Location loc);
void ImplicitIterationMaterial(InputPara* IP, Bound* bound, IterMaterial* IM, enum Location loc);

#endif
