Penghan CHEN (陈鹏翰) Hao HUANG (黄灏)

### Problem

<img src="./fig1 2022-05-28 14.53.52.png" style="zoom: 77%;" />

### Solutions

We learned the solution methods from the [videos](https://www.bilibili.com/video/BV13A411B7n9?spm_id_from=333.337.search-card.all.click) on bilibili and the books from class MAE5032-HPC by [Ju-LIU](https://ju-liu.github.io/)

<img src="./fig2 2022-05-28 16.23.01.png" style="zoom:80%;" />

$u_{l,m}^t$ is denoted as temperature u in node l,m at time t.

#### (1) For each internal node, the $u_{l,m}^{t+1}$ should be depended on five data at time t (i.e., $u_{l,m}^{t}$, $u_{l+1,m}^{t}$, $u_{l-,m}^{t}$, $u_{l,m+1}^{t}$, $u_{l,m-1}^{t}$) spatially and temporally.

**The basic equation** should be
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) = \kappa((\dfrac{\partial u}{\partial x})_{l+1,m}^{t\rightarrow t+1} - (\dfrac{\partial u}{\partial x})_{l-1,m}^{t\rightarrow t+1})\Delta y\Delta t + \kappa((\dfrac{\partial u}{\partial y})_{l,m+1}^{t\rightarrow t+1} - (\dfrac{\partial u}{\partial y})_{l,m-1}^{t\rightarrow t+1})\Delta x\Delta t
$$
**In explicit format**, the relationship can be decribed by the following.
$$
\dfrac{\rho c \Delta x \Delta y}{\Delta t}(u_{l,m}^{t+1}-u_{l,m}^{t}) &= \kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}} -  \dfrac{u_{l,m}^{t}-u_{l-1,m}^{t}}{\delta x_{l-1}})\Delta y + \kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}} - \dfrac{u_{l,m}^{t}-u_{l,m-1}^{t}}{\delta y_{m+1}})\Delta x\\
$$
$$
a_{l,m}^{t+1} u_{l,m}^{t+1} = a_{l+1,m}^{t}u_{l+1,m}^{t} + a_{l-1,m}^{t}u_{l-1,m}^{t} + a_{l,m+1}^{t}u_{l,m+1}^{t} + a_{l,m-1}^{t}u_{l,m-1}^{t} + a_{l,m}^{t}u_{l,m}^{t}
$$

where
$$
\begin{aligned}
a_{l+1,m}^{t} &= \dfrac{\kappa\Delta y}{\delta x_{l+1}} = \dfrac{\kappa\Delta y}{\Delta x} = a_{l-1,m}^{t} = \dfrac{\kappa\Delta y}{\delta x_{l-1}}\\
a_{l,m+1}^{t} &= \dfrac{\kappa\Delta x}{\delta y_{m+1}} = \dfrac{\kappa\Delta x}{\Delta y} = a_{l,m-1}^{t} = \dfrac{\kappa\Delta x}{\delta y_{m-1}}\\
a_{l,m}^{t} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t}  - \dfrac{\kappa\Delta y}{\delta x_{l+1}} - \dfrac{\kappa\Delta y}{\delta x_{l-1}} - \dfrac{\kappa\Delta x}{\delta y_{m+1}} - \dfrac{\kappa\Delta x}{\delta y_{m-1}}\\
&= \dfrac{\rho c \Delta x \Delta y}{\Delta t}  - 2\dfrac{\kappa\Delta y}{\Delta x} - 2\dfrac{\kappa\Delta x}{\Delta y}\\
a_{l,m}^{t+1} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t} = 2a_{l+1,m}^{t} + 2a_{l,m+1}^{t} + a_{l,m}^{t}\\
\end{aligned}
$$
Note that $\delta$ is actually the same as $\Delta$ for internal nodes.

**In implicit format**, the relationship can be decribed by the following.
$$
\dfrac{\rho c \Delta x \Delta y}{\Delta t}(u_{l,m}^{t+1}-u_{l,m}^{t}) &= \kappa(\dfrac{u_{l+1,m}^{t+1}-u_{l,m}^{t+1}}{\delta x_{l+1}} - \dfrac{u_{l,m}^{t+1}-u_{l-1,m}^{t+1}}{\delta x_{l-1}})\Delta y + \kappa(\dfrac{u_{l,m+1}^{t+1} - u_{l,m}^{t+1}}{\delta y_{m+1}} - \dfrac{u_{l,m}^{t+1}-u_{l,m-1}^{t+1}}{\delta x_{m-1}})\Delta x\\
$$

$$
a_{l,m}^{t+1} u_{l,m}^{t+1} - a_{l+1,m}^{t+1}u_{l+1,m}^{t+1} - a_{l-1,m}^{t+1}u_{l-1,m}^{t+1} - a_{l,m+1}^{t+1}u_{l,m+1}^{t+1} - a_{l,m-1}^{t+1}u_{l,m-1}^{t+1} = a_{l,m}^{t}u_{l,m}^{t}
$$

where
$$
\begin{aligned}
a_{l+1,m}^{t+1} &= \dfrac{\kappa\Delta y}{\delta x_{l+1}} = \dfrac{\kappa\Delta y}{\Delta x}\\
a_{l-1,m}^{t+1} &= \dfrac{\kappa\Delta y}{\delta x_{l-1}} = \dfrac{\kappa\Delta y}{\Delta x}\\
a_{l,m+1}^{t+1} &= \dfrac{\kappa\Delta x}{\delta y_{m+1}} = \dfrac{\kappa\Delta x}{\Delta y}\\
a_{l,m-1}^{t+1} &= \dfrac{\kappa\Delta x}{\delta y_{m-1}} = \dfrac{\kappa\Delta x}{\Delta y}\\
a_{l,m}^{t} &= \dfrac{\rho c \Delta x \Delta y}{\Delta t}\\
a_{l,m}^{t+1} &= a_{l+1,m}^{t+1} + a_{l-1,m}^{t+1} + a_{l,m+1}^{t+1} + a_{l,m-1}^{t+1} +a_{l,m}^{t}\\
\end{aligned}
$$

#### (2) Considering the boundary condition

When $l,m\in \Gamma$, data are given that $u=g$ or $\kappa\dfrac{\partial u}{\partial x}n_x+\kappa\dfrac{\partial u}{\partial y}n_y=h$.

Starting from only considering $h$ from the left and down sides, the basic equation can be decribed as the following.
$$
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) = (\kappa(\dfrac{\partial u}{\partial x})_{l+1,m}^{t\rightarrow t+1} + h_{l-1,m}^{t\rightarrow t+1})\Delta y\Delta t + (\kappa(\dfrac{\partial u}{\partial y})_{l,m+1}^{t\rightarrow t+1} + h_{l,m-1}^{t\rightarrow t+1})\Delta x\Delta t
$$
**In explicit format**, 
$$
\begin{aligned}
\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t}) &= \kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}} +  h_{l-1,m}^{t})\Delta y\Delta t + \kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}} + h_{l,m-1}^{t})\Delta x\Delta t\\
\dfrac{\rho c \Delta x \Delta y(u_{l,m}^{t+1}-u_{l,m}^{t})}{\Delta t} &= \kappa(\dfrac{u_{l+1,m}^{t}-u_{l,m}^{t}}{\delta x_{l+1}} +  h_{l-1,m}^{t})\Delta y + \kappa(\dfrac{u_{l,m+1}^{t}-u_{l,m}^{t}}{\delta y_{m+1}} + h_{l,m-1}^{t})\Delta x\\
a_{l,m}^{t+1} u_{l,m}^{t+1} &= a_{l+1,m}^{t}u_{l+1,m}^{t} + \kappa h_{l-1,m}^{t}\Delta y + a_{l,m+1}^{t}u_{l,m+1}^{t} + \kappa h_{l,m-1}^{t}\Delta x + a_{l,m}^{t}u_{l,m}^{t}\\
\end{aligned}
$$
We described this in general format.
$$
\begin{aligned}
a_{l,m}^{t+1} u_{l,m}^{t+1} &= a_{l+1,m}^{t}u_{l+1,m}^{t} + a_{l-1,m}^{t}u_{l-1,m}^{t} + a_{l,m+1}^{t}u_{l,m+1}^{t} + a_{l,m-1}^{t}u_{l,m-1}^{t} \\
&+ \kappa h_{l+1,m}^{t}\Delta y + \kappa h_{l-1,m}^{t}\Delta y + \kappa h_{l,m+1}^{t}\Delta x + \kappa h_{l,m-1}^{t}\Delta x + a_{l,m}^{t}u_{l,m}^{t}
\end{aligned}
$$
**In implicit format**, similarly we get that
$$
\begin{aligned}
a_{l,m}^{t+1}u_{l,m}^{t+1} &= a_{l+1,m}^{t+1}u_{l+1,m}^{t+1} + a_{l-1,m}^{t+1}u_{l-1,m}^{t+1} + a_{l,m+1}^{t+1}u_{l,m+1}^{t+1} + a_{l,m-1}^{t+1}u_{l,m-1}^{t+1} \\
&+ \kappa h_{l+1,m}^{t+1}\Delta y + \kappa h_{l-1,m}^{t+1}\Delta y + \kappa h_{l,m+1}^{t+1}\Delta x + \kappa h_{l,m-1}^{t+1}\Delta x + a_{l,m}^{t}u_{l,m}^{t}\\
\end{aligned}
$$

##### (a) If prescribed $u=g$ for the boundary nodes in $\Omega$

$h_{l+1,m}=h_{l-1,m}=h_{l,m+1}=h_{l,m-1}=0$;

when $l+\Delta x=1$, $u_{l+1,m}=g$;

when $l-\Delta x=0$, $u_{l-1,m}=g$;

when $m+\Delta y=0$, $u_{l,m+1}=g$;

when $m-\Delta y=0$, $u_{l,m-1}=g$.

##### (b) If prescribed $\kappa\dfrac{\partial u}{\partial x}n_x+\kappa\dfrac{\partial u}{\partial y}n_y=h$ for the boundary nodes in $\Omega$

when $l+1=1$ , $h_{l+1,m}=-h$ and $a_{l+1,m}=h_{l-1,m}=h_{l,m+1}=h_{l,m-1}=0$; 

when $l-1=0$ , $h_{l-1,m}=h$ and $a_{l-1,m}=h_{l+1,m}=h_{l,m+1}=h_{l,m-1}=0$; 

when $m+1=1$,  $h_{l,m+1}=-h$ and $a_{l,m+1}=h_{l+1,m}=h_{l-1,m}=h_{l,m-1}=0$; 

when $m-1=0$,  $h_{l,m-1}=h$ and $a_{l,m-1}=h_{l+1,m}=h_{l-1,m}=h_{l,m+1}=0$.



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



