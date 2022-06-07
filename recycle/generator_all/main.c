static char help[] = "Generate the  input data. \n \
            -g : all of boundary are given g \n\
            -h : all of boundary are given h\n \
            -gl: specify the left face of g \n \
            -gr: specify the left face of g \n \
            -gt: specify the left face of g \n \
            -gb: specify the left face of g \n \
            -hl: specify the left face of h \n \
            -hr: specify the left face of h \n \
            -ht: specify the left face of h \n \
            -hb: specify the left face of h \n \
            *g and h cann't be given for the same point\n *\
            \
            -f : heat supply\n \
            -u0 : the initial value of u\n \
            -n : the scale of the problem\n \
            -dl : the space resolution of the problem, it should be 1/n\n\
            -dt : the time resolution of the problem\n \
            -k : the conductivity\n \
            -rho : the density\n \
            -c : heat capacity";

#include <petsc.h>

int main(int argc, char** argv)
{
    PetscInt        i, n, maxIts;
    PetscScalar     g, gl, gr, gt, gb;
    PetscScalar     h, hl, hr, ht, hb;
    PetscScalar     u0;
    PetscScalar     f, dl, dt, k, rho, c;
    PetscMPIInt     rank;
    MPI_Comm        comm;
    PetscErrorCode  ierr;
    PetscViewer     viewer;
    PetscChar       fname[PETSC_MAX_PATH_LEN]="default.hdf5";

    Vec             g_b,g_t,g_l,g_r;
    Vec             h_b,h_t,h_l,h_r;
    Vec             t_b,t_t,t_l,t_r;
    Vec             u_0;

    Vec             paras;    
    PetscInt        ip[7];
    PetscScalar     vp[7];

    g=gl=gr=gt=gb=0; //specify -1 for random input
    h=hl=hr=ht=hb=0;

    PetscInitialize(&argc, &argv, (char*) 0, help);  
    comm = PETSC_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);CHKERRMPI(ierr);

    PetscOptionsGetInt(NULL,NULL, "-n", &n, NULL);

    PetscOptionsGetScalar(NULL,NULL,"-u0",&u0,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-g", &g, NULL);
    PetscOptionsGetScalar(NULL,NULL,"-h", &h, NULL);

    PetscOptionsGetScalar(NULL,NULL,"-gl",&gl,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-gr",&gr,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-gt",&gt,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-gb",&gb,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-hl",&hl,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-hr",&hr,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-ht",&ht,NULL);
    PetscOptionsGetScalar(NULL,NULL,"-hb",&hb,NULL);

    PetscOptionsGetScalar(NULL,NULL,"-dt",  &dt,  NULL);
    PetscOptionsGetScalar(NULL,NULL,"-dl",  &dl,  NULL);
    PetscOptionsGetScalar(NULL,NULL,"-rho", &rho, NULL);
    PetscOptionsGetScalar(NULL,NULL,"-c",   &c,   NULL); 
    PetscOptionsGetScalar(NULL,NULL,"-k",   &k,   NULL);
    PetscOptionsGetScalar(NULL,NULL,"-f",   &f,   NULL);    

    PetscOptionsGetString(NULL,NULL,"-fname",fname,sizeof(fname),NULL);

    VecCreate(comm,&paras);
    VecSetSizes(paras,PETSC_DECIDE,7);
    VecSetFromOptions(paras);

    VecCreate(comm,&g_b);
    VecSetSizes(g_b,PETSC_DECIDE,n);
    VecSetFromOptions(g_b);
    VecDuplicate(g_b,&h_b);
    VecDuplicate(g_b,&t_b);

    VecDuplicate(g_b,&g_t);
    VecDuplicate(g_b,&h_t);
    VecDuplicate(g_b,&t_t);

    VecDuplicate(g_b,&g_l);
    VecDuplicate(g_b,&h_l);
    VecDuplicate(g_b,&t_l);

    VecDuplicate(g_b,&g_r);
    VecDuplicate(g_b,&h_r);
    VecDuplicate(g_b,&t_r);

    VecCreate(comm,&u_0);
    VecSetSizes(u_0,PETSC_DECIDE,n*n);
    VecSetFromOptions(u_0);

    if(dl*n!=1.0)
    {
        PetscPrintf(comm,"n and dl should be reciprocals of each other.\n");
        return -1;
    }
    
    if(rank==0)
    {   
        for(i=0;i<7;i++)
            ip[i]=i;
        vp[0]=dt;   vp[1]=dl;   vp[2]=rho;
        vp[3]=c;    vp[4]=k;    vp[5]=f;
        vp[6]=(PetscScalar)n;
        VecSetValues(paras,7,ip,vp,INSERT_VALUES);
        printf("Setting finished for parameters. dt:  %g\n",dt);
        printf("Setting finished for parameters. dl:  %g\n",dl);
        printf("Setting finished for parameters. rho: %g\n",rho);
        printf("Setting finished for parameters. c:   %g\n",c);
        printf("Setting finished for parameters. k:   %g\n",k);
        printf("Setting finished for parameters. f:   %g\n",f);
        printf("Setting finished for parameters. n:   %d\n",n);
    }
    if(u0<0)
    {
        PetscPrintf(comm,"Temperature(K) can not be negative.\n");
        return -1;
    }
    else
    {
        VecSet(u_0,u0);
        PetscPrintf(comm,"Setting finished for the initial u0: %g.\n",u0);
    }

    if(g==0&&h==0)
    {
        if((gl==0&&hl==0)||(gr==0&&hr==0)||(gb==0&&hb==0)||(gt==0&&ht==0))
        {
            PetscPrintf(comm,"You need to specify h or g for the boundary.\n");
            return -1;
        }
            
        else if ((gl&&hl)||(gr&&hr)||(gb&&hb)||(gt&&ht))
        {
             PetscPrintf(comm,"You can not specify g and h at same point.\n");
             return -1;
        }
        else
        {
            if(gl!=0)
            {
                VecSet(g_l,gl);
                VecSet(h_l,0.0);
                VecSet(t_l,1.0);
                PetscPrintf(comm,"Setting finished for the left face of g: %g.\n",gl);
            }
            else
            {
                VecSet(g_l,0.0);
                VecSet(h_l,hl);
                VecSet(t_l,2.0);
                PetscPrintf(comm,"Setting finished for the left face of h: %g.\n",hl);
            }

            if(gr!=0)
            {
                VecSet(g_r,gr);
                VecSet(h_r,0.0);
                VecSet(t_r,1.0);
                PetscPrintf(comm,"Setting finished for the right face of g: %g.\n",gr);
            }
            else
            {
                VecSet(g_r,0.0);
                VecSet(h_r,hr);
                VecSet(t_r,2.0);
                PetscPrintf(comm,"Setting finished for the right face of h: %g.\n",hr);
            }            

            if(gt!=0)
            {
                VecSet(g_t,gt);
                VecSet(h_t,0.0);
                VecSet(t_t,1.0);
                PetscPrintf(comm,"Setting finished for the top face of g: %g.\n",gt);
            }
            else
            {
                VecSet(g_t,0.0);
                VecSet(h_t,ht);
                VecSet(t_t,2.0);
                PetscPrintf(comm,"Setting finished for the top face of h: %g.\n",ht);
            }    

            if(gb!=0)
            {
                VecSet(g_b,gb);
                VecSet(h_b,0.0);
                VecSet(t_b,1.0);
                PetscPrintf(comm,"Setting finished for the bottom face of g: %g.\n",gb);
            }
            else
            {
                VecSet(g_b,0.0);
                VecSet(h_b,hb);
                VecSet(t_b,2.0);
                PetscPrintf(comm,"Setting finished for the bottom face of h: %g.\n",hb);
            }                    
        }
    }
    else if(g!=0&&h!=0)
    {
        PetscPrintf(comm,"You can not specify g and h at same point.\n");
        return -1;
    }
    else if(g!=0)
    {
        VecSet(g_b,g);  VecSet(h_b,0.0);    VecSet(t_b,1.0);
        VecSet(g_t,g);  VecSet(h_t,0.0);    VecSet(t_t,1.0);
        VecSet(g_l,g);  VecSet(h_l,0.0);    VecSet(t_l,1.0);
        VecSet(g_r,g);  VecSet(h_r,0.0);    VecSet(t_r,1.0);
        PetscPrintf(comm,"Setting finished for g of all the boundary: %g.\n",g);

    }
    else 
    {
        VecSet(g_b,0.0);  VecSet(h_b,h);    VecSet(t_b,2.0);
        VecSet(g_t,0.0);  VecSet(h_t,h);    VecSet(t_t,2.0);
        VecSet(g_l,0.0);  VecSet(h_l,h);    VecSet(t_l,2.0);
        VecSet(g_r,0.0);  VecSet(h_r,h);    VecSet(t_r,2.0);
        PetscPrintf(comm,"Setting finished for h of all the boundary: %g.\n",h);
    }

    PetscViewerHDF5Open(PETSC_COMM_WORLD,fname,FILE_MODE_WRITE,&viewer);
    PetscViewerHDF5PushGroup(viewer,"/boundary");

    PetscObjectSetName((PetscObject)g_b,"g_bottom");
    VecView(g_b,viewer);
    PetscObjectSetName((PetscObject)h_b,"h_bottom");
    VecView(h_b,viewer);
    PetscObjectSetName((PetscObject)t_b,"t_bottom");
    VecView(t_b,viewer);

    PetscObjectSetName((PetscObject)g_t,"g_top");
    VecView(g_t,viewer);   
    PetscObjectSetName((PetscObject)h_t,"h_top");
    VecView(h_t,viewer);  
    PetscObjectSetName((PetscObject)t_t,"t_top");
    VecView(t_t,viewer);

    PetscObjectSetName((PetscObject)g_l,"g_left");
    VecView(g_l,viewer);  
    PetscObjectSetName((PetscObject)h_l,"h_left");
    VecView(h_l,viewer);  
    PetscObjectSetName((PetscObject)t_l,"t_left");
    VecView(t_l,viewer);

    PetscObjectSetName((PetscObject)g_r,"g_right");
    VecView(g_r,viewer);  
    PetscObjectSetName((PetscObject)h_r,"h_right");
    VecView(h_r,viewer);  
    PetscObjectSetName((PetscObject)t_r,"t_right");
    VecView(t_r,viewer);
    PetscViewerHDF5PopGroup(viewer);

    PetscViewerHDF5PushGroup(viewer,"/Init");
    PetscObjectSetName((PetscObject)u_0,"u_0");
    VecView(u_0,viewer);    
    PetscViewerHDF5PopGroup(viewer);

    PetscViewerHDF5PushGroup(viewer,"/Parameters");
    PetscObjectSetName((PetscObject)paras,"paras");
    VecView(paras,viewer);    
    PetscViewerHDF5PopGroup(viewer);
    
    PetscPrintf(comm,"HDF5 file %s has been written.\n",fname);

    PetscFinalize();
    
    return 0;
}