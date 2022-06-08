# HTS2D(MAE5032_FinalProject)
![4fe036dab33d34df8826e296532eed5](https://perhaps-1306015279.cos.ap-guangzhou.myqcloud.com/4fe036dab33d34df8826e296532eed5.jpg)
HTS2D使用git进行版本控制，并上传到github仓库，使用以下命令克隆至本地

```bash
git clone https://github.com/PerhapsChen/MAE5032_FinalProject
```

克隆完成后，进入build文件夹

```bash
cd build
```

加载相关的依赖后使用 `make` 进行构建得到 `HST2D`可执行文件。

```bash
module load purge
module load mpi/intel/2018.4

make HTS2D
```

在build文件夹下，`HST2D`的基本使用格式

```
mpirun -np X ./HST2D [FUNCTION] [PARAMETERS] [OTHERS]
```

其中 X 表示使用的核心数

- `FUNCTION`为要使用的功能，目前包括
  - generator:  生成输入数据，包括边界条件，初始温度，其他常量等。
  - explicit:  显式迭代计算，可以指定最大迭代次数，可以重启
  - implicit:  隐式迭代计算，可以指定最大迭代次数，可以重启
- `PARAMETERS`功能要用到的输入参数，具体见下面部分的介绍。
- `OTHERS`表示其他flag，比如 `-log_view`等

推荐使用太乙脚本运行

### 3.1 generator

generator可以生成计算所需要的HDF5文件，文件中记载了包括边界条件，各物理常量值，每次迭代得到的结果以及迭代次数等信息。

使用generator需要按照以下格式

```
mpirun -np X ./HTS2D generator [PARAMETERS] [OTHERS]
```

参数说明

| PARAMETERS | 说明                                                     |
| ---------- | -------------------------------------------------------- |
| -fname     | 指定文件名                                               |
| -n         | 指定矩阵规模，为空间步长dl的倒数                         |
| -g         | 指定所有边界给定g，使用给定的g的值，不能为负             |
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

例如，使用单核处理器，生成时间步长为1，空间步长为0.1，规模为10，上边界、左边界温度值为200，有边界、下边界热流密度为100，密度为300，head capacity 为1000，conductivity 为10 ，heat supply 为0，初始温度平面为50的文件名为default.hdf5的命令为

```bash
mpirun -np 1 ./HTS2D generator -fname ../data/default.hdf5 -n 10 -dl 0.1 -dt 1 -gl 200 -gt 200 -hr 100 -hb 100 -rho 300 -c 1000 -k 10 -f 0 -u0 50
```

![image-20220608183941646](https://perhaps-1306015279.cos.ap-guangzhou.myqcloud.com/image-20220608183941646.png)

生成的文件如上图所示，其中u_0中记录了每次时间步长迭代得到的温度值， Parameters记录了相关常量物理参数，boundary记录了边界条件，g表示温度，h表示热流密度，t指出该点的边界类型（1代表狄利克雷边界，2代表诺依曼边界条件）

### 3.2 Explicit

根据generator生成的数据，我们可以用来选择显式还是隐式计算每次时间步长迭代得到的温度值

使用以下命令选择使用显式计算

```
mpirun -np X ./HTS2D explicit [PARAMETERS] [OTHERS]
```

参数说明

| PARAMETERS | 说明                               |
| ---------- | ---------------------------------- |
| -fname     | 指定存储信息的HDF5文件路径及文件名 |
| -maxItsW   | 多次重启任务中，最大迭代次数       |
| -maxIts    | 该次任务中，累计最大迭代次数       |
| -restart   | 0代表非重启型任务，1代表重启任务   |

例如，我们希望使用4个处理器内核对3.1中生成的数据进行迭代计算，设置该任务总迭代次数最多为1000次，该次迭代次数为100次（不希望一次性计算完毕），首次运行时候，参考以下命令

```bash
mpirun -np 4 ./HTS2D explicit -fname ../data/default.hdf5 -maxItsW 1000 -maxIts 100 -restart 0
```

即可在该HDF5文件中记录每次迭代的结果及迭代次数，存放在u_t的group中

![image-20220608185234090](https://perhaps-1306015279.cos.ap-guangzhou.myqcloud.com/image-20220608185234090.png)

如果我们希望重启任务，从上次迭代结束的次数重新计算至500次，则可以考虑以下命令

```bash
mpirun -np 4 ./HTS2D explicit -fname ../data/default.hdf5 -maxItsW 1000 -maxIts 500 -restart 1
```

![image-20220608185454934](https://perhaps-1306015279.cos.ap-guangzhou.myqcloud.com/image-20220608185454934.png)

### 3.3 Implicit

隐式计算与显式计算相同，只需要将`[FUNCTION]`项修改为 `implicit`

```bash
# 生成数据
mpirun -np 1 ./HTS2D generator -fname ../data/implicit_test.hdf5 -n 10 -dl 0.1 -dt 1 -g 200 -rho 300 -c 1000 -k 10 -f 0 -u0 50

# 首次迭代100次
mpirun -np 4 ./HTS2D implicit -fname ../data/implicit_test.hdf5 -maxItsW 1000 -maxIts 100 -restart 0

# 重启迭代至500次
mpirun -np 4 ./HTS2D implicit -fname ../data/implicit_test.hdf5 -maxItsW 1000 -maxIts 500 -restart 1
```