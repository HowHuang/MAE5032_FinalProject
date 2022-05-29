Penghan CHEN (陈鹏瀚) Hao HUANG (黄灏)

#### Problem

<img src="/Users/haohuang/Documents/3 PhD study/2 Courses/HPC/Homework/MAE5032_FinalProject/fig1 2022-05-28 14.53.52.png" alt="截屏2022-05-28 14.53.52" style="zoom: 77%;" />

#### Solutions

We learned the solution methods from the videos on [bilibili](https://www.bilibili.com/video/BV13A411B7n9?spm_id_from=333.337.search-card.all.click) and the books from class MAE5032-HPC by [Ju-LIU](https://ju-liu.github.io/)

<img src="/Users/haohuang/Documents/3 PhD study/2 Courses/HPC/Homework/MAE5032_FinalProject/fig2 2022-05-28 16.23.01.png" alt="截屏2022-05-28 16.23.01" style="zoom:80%;" />

$u_{l,m}^t$ is denoted as temperature u in node l,m at time t.

##### (1) For each node, the $u_{l,m}^{t+1}$ should be depended on five data at time t (i.e., $u_{l,m}^{t}$, $u_{l+1,m}^{t}$, $u_{l-,m}^{t}$, $u_{l,m+1}^{t}$, $u_{l,m-1}^{t}$) spatially and temporally.

**The basic equation** should be
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) = \kappa((\dfrac{\partial u}{\partial x})_{l+1,m}^{t\rightarrow t+1} - (\dfrac{\partial u}{\partial x})_{l-1,m}^{t\rightarrow t+1})\Delta y\Delta t + \kappa((\dfrac{\partial u}{\partial y})_{l,m+1}^{t\rightarrow t+1} - (\dfrac{\partial u}{\partial y})_{l,m-1}^{t\rightarrow t+1})\Delta x\Delta t
$$
As there is $\kappa\dfrac{\partial u}{\partial x}n_x+\kappa\dfrac{\partial u}{\partial y}n_y=h$, we considered this is a Neumann problem.
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) = (\kappa(\dfrac{\partial u}{\partial x})_{l+1,m}^{t\rightarrow t+1} + h_{l-1,m}^{t\rightarrow t+1})\Delta y\Delta t + (\kappa(\dfrac{\partial u}{\partial y})_{l,m+1}^{t\rightarrow t+1} + h_{l,m-1}^{t\rightarrow t+1})\Delta x\Delta t
$$
**In explicit format**, the relationship can be decribed by the following.
$$
\begin{aligned}
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) &= \kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}} +  h_{l-1,m}^{t})\Delta y\Delta t + \kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}} + h_{l,m-1}^{t})\Delta x\Delta t\\
\dfrac{\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t})}{\Delta t} &= \kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}} +  h_{l-1,m}^{t})\Delta y + \kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}} + h_{l,m-1}^{t})\Delta x\\
\end{aligned}
$$
$$
a_{l,m}^{t+1} u_{l,m}^{t+1} = a_{l+1,m}^{t}u_{l+1,m}^{t} + \kappa h_{l-1,m}^{t}\Delta y + a_{l,m+1}^{t}u_{l,m+1}^{t} + \kappa h_{l,m-1}^{t}\Delta x + a_{l,m}^{t}u_{l,m}^{t}
$$

where
$$
\begin{aligned}
a_{l+1,m}^{t} &= \dfrac{\kappa\Delta y}{\delta x_{l+1}}\\
a_{l,m+1}^{t} &= \dfrac{\kappa\Delta x}{\delta y_{m+1}}\\
a_{l,m}^{t} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t} - \dfrac{\kappa\Delta y}{\delta x_{l+1}} - \dfrac{\kappa\Delta x}{\delta y_{m+1}}\\
a_{l,m}^{t+1} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t} = a_{l+1,m}^{t} + a_{l,m+1}^{t}\\
\end{aligned}
$$
Note that $\delta$ is actually the same as $\Delta$ for internal nodes, and $\delta$ is the half of $\Delta$ for boundary nodes.

**In implicit format**, the relationship can be decribed by the following.
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) = \kappa(\dfrac{u_{l+1,m}^{t+1}-u_{l,m}^{t+1}}{\delta x_{l+1}} - \dfrac{u_{l,m}^{t+1}-u_{l-1,m}^{t+1}}{\delta x_{l-1}})\Delta y\Delta t + \kappa(\dfrac{u_{l,m+1}^{t+1} - u_{l,m}^{t+1}}{\delta y_{m+1}} - \dfrac{u_{l,m}^{t+1}-u_{l,m-1}^{t+1}}{\delta x_{m-1}})\Delta x\Delta t\\
\dfrac{\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t})}{\Delta t} = \kappa(\dfrac{u_{l+1,m}^{t+1}-u_{l,m}^{t+1}}{\delta x_{l+1}} +  h_{l-1,m}^{t+1})\Delta y + \kappa(\dfrac{u_{l,m+1}^{t+1}-u_{l,m}^{t+1}}{\delta y_{m+1}} + h_{l,m-1}^{t+1})\Delta x\\
$$

$$
a_{l,m}^{t+1} u_{l,m}^{t+1}=a_{l+1,m}^{t+1}u_{l+1,m}^{t+1} + \kappa h_{l-1,m}^{t+1}\Delta y + a_{l,m+1}^{t+1}u_{l,m+1}^{t+1} + \kappa h_{l,m-1}^{t+1}\Delta x + a_{l,m}^{t}u_{l,m}^{t}
$$

where
$$
\begin{aligned}
a_{l+1,m}^{t+1} &= \dfrac{\kappa\Delta y}{\delta x_{l+1}}\\
a_{l,m+1}^{t+1} &= \dfrac{\kappa\Delta x}{\delta y_{m+1}}\\
a_{l,m}^{t} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t}\\
a_{l,m}^{t+1} - a_{l+1,m}^{t+1} - a_{l,m+1}^{t+1} &= a_{l,m}^{t}\\
\end{aligned}
$$




##### (2) For boundary condition, as there is $\kappa\dfrac{\partial u}{\partial x}n_x+\kappa\dfrac{\partial u}{\partial y}n_y=h$, we considered this is a Neumann problem. 



















### Targets (Course Video 15)

#### 编写程序

- [ ] 程序可以重启（需要hfd5来帮我们）
- [ ] 防御性写法，能处理异常，有断言，有详细的注释
- [ ] 编译的时候不能有任何warning
- [ ] 有makfile或者cmake，可以在太乙上直接编译
- [ ] 有profiling，程序的哪里耗时最多
- [ ] 利用github进行版本控制和团队合作
- [ ] 用开源软件可视化，二维的话需要用VTK或者Paraview
- [ ] 技术性的报告，最好包括LaTex的使用

#### 测试

- [ ] 方法的稳定性
- [ ] 误差分析
- [ ] 并行效率（强可扩展性、弱可扩展性、Petsc不同求解器的影响）



































