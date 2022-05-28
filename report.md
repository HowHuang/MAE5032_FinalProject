Penghan CHEN (陈鹏瀚) Hao HUANG (黄灏)

#### Problem

<img src="/Users/haohuang/Documents/3 PhD study/2 Courses/HPC/Homework/MAE5032_FinalProject/fig1 2022-05-28 14.53.52.png" alt="截屏2022-05-28 14.53.52" style="zoom: 77%;" />

#### Solutions

We learned the solution methods from the videos on [bilibili](https://www.bilibili.com/video/BV13A411B7n9?spm_id_from=333.337.search-card.all.click) and the books from class MAE5032-HPC-by_Ju-LIU (https://ju-liu.github.io/)

![截屏2022-05-28 16.23.01](/Users/haohuang/Documents/3 PhD study/2 Courses/HPC/Homework/MAE5032_FinalProject/fig2 2022-05-28 16.23.01.png)

$u_{l,m}^t$ is denoted as temperature u in node l,m at time t.

**For internal node**, the $u_{l,m}^{t+1}$ should be depended on five data at time t (i.e., $u_{l,m}^{t}$, $u_{l+1,m}^{t}$, $u_{l-,m}^{t}$, $u_{l,m+1}^{t}$, $u_{l,m-1}^{t}$) spatially and temporally, whose relationship can be decribed by the following,
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t})-\kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}}-\dfrac{u_{l,m}^{t}-u_{l-1,m}^{t}}{\delta x_{l-1}})-\kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}}-\dfrac{u_{l,m}^{t}-u_{l,m-1}^{t}}{\delta x_{m-1}})
$$
where $\delta$ is actually the same as $\Delta$ for internal node.

**For boundary node**, as there is $\kappa\dfrac{\partial u}{\partial x}n_x+\kappa\dfrac{\partial u}{\partial y}n_y=h$, we considered this is a Neumann problem

