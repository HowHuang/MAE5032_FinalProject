# MAE5032_FinalProject
## Gnerator

需要指定FUNCTIONS为generator

```
mpirun -np 4 ./main.out generator [PARAMETERS]
```

**PARAMETERS**

| PARAMETERS | 说明                                                     |
| ---------- | -------------------------------------------------------- |
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
| -u0        | 指定初始平面温度                                         |
| -fname     | 指定文件名                                               |
| -n         | 指定矩阵规模，为dl的倒数                                 |



## Explicit

需要指定FUNCTIONS为explicit

```
mpirun -np 4 ./main.out explicit [PARAMETERS]
```

**PARAMETERS**

| PARAMETERS | 说明                             |
| ---------- | -------------------------------- |
| -maxItsW   | 多次重启任务中，最大迭代次数     |
| -maxIts    | 该次任务中，累计最大迭代次数     |
| -restart   | 0代表非重启型任务，1代表重启任务 |
| -fname     | 指定存储信息的HDF5文件名         |

