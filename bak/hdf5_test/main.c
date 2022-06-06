#include <petsc.h>
#include <hdf5.h>

int main()
{
    const int n = 10;

    PetscScalar     *g_b, *g_t, *g_l, *g_r;
    PetscScalar     *h_b, *h_t, *h_l, *h_r;

    PetscMalloc4(n, &g_b, n, &g_t, n, &g_l, n, &g_r);
    PetscMalloc4(n, &h_b, n, &h_t, n, &h_l, n, &h_r);

    hid_t   file_id, group_id, dataset_id;
    file_id=H5Fopen("g_fixed.hdf5",H5F_ACC_RDONLY,H5P_DEFAULT);
    group_id=H5Gopen(file_id,"/boundary",H5P_DEFAULT);
    dataset_id=H5Dopen(group_id,"g_top",H5P_DEFAULT);

    H5Dread(dataset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,g_t);
    int i = 0;
    for(i=0;i<n;++i)
    {
        printf("g_t[%d]=%g \n",i,g_t[i]);
    }
    return 0;
}