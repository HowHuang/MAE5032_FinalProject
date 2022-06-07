# MAE5032_FinalProject

MAE5032课程，二维数值传热软件

在本地调试请使用该命令格式，在太乙上使用请使用专属脚本。


```bash
mpirun -np X ./program [FUNCTION] [PARAMETERS]
```

- `FUNCTION`为要使用的功能，目前包括
  - generator:  生成输入数据，包括边界条件，初始温度，其他常量等。
  - explicit:  显式迭代计算，可以指定最大迭代次数，可以重启
  - implicit:  隐式迭代计算，可以指定最大迭代次数，可以重启
- `PARAMETERS`功能要用到的输入参数，具体见下面部分的介绍。
## Generator

需要指定FUNCTIONS为generator

```bash
mpirun -np 4 ./main.out generator [PARAMETERS]
```

**PARAMETERS**

| PARAMETERS | 说明                                                     |
| ---------- | -------------------------------------------------------- |
| -fname     | 指定文件名                                               |
| -n         | 指定矩阵规模，为dl的倒数                                 |
| -g         | 指定所有边界给定g，使用给定的g的值                       |
| -h         | 指定所有边界给定h，使用给定的h的值h                      |
| -gl        | 指定左边界给定g，其他边同理(gr, gt, gb)                  |
| -hl        | 指定左边界给定h，其他边同理(hr, ht, hb)                  |
|            | 不同的边可以指定h或者g，例如 -gl 10 -hr 20 -ht 30 -gb 35 |
| -dt        | 指定时间步长(s)                                          |
| -dl        | 指定空间步长                                             |
| -rho       | 指定密度                                                 |
| -c         | 指定heat capacity                                        |
| -k         | 指定conductivity                                         |
| -f         | 指定heat supply                                          |
| -u0        | 指定初始平面温度，需要大于0                              |



## Explicit

需要指定FUNCTIONS为explicit

```
mpirun -np 4 ./main.out explicit [PARAMETERS]
```

**PARAMETERS**

| PARAMETERS | 说明                             |
| ---------- | -------------------------------- |
| -fname     | 指定存储信息的HDF5文件名         |
| -maxItsW   | 多次重启任务中，最大迭代次数     |
| -maxIts    | 该次任务中，累计最大迭代次数     |
| -restart   | 0代表非重启型任务，1代表重启任务 |



## Implicit

需要指定FUNCTIONS为implicit

```bash
mpirun -np 4 ./main.out imxplicit [PARAMETERS]
```

**PARAMETERS**

| PARAMETERS | 说明                             |
| ---------- | -------------------------------- |
| -fname     | 指定存储信息的HDF5文件名         |
| -maxItsW   | 多次重启任务中，最大迭代次数     |
| -maxIts    | 该次任务中，累计最大迭代次数     |
| -restart   | 0代表非重启型任务，1代表重启任务 |

