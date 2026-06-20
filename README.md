# GPS 伪距单点定位程序

本项目包含 3 个 Python 文件，用于读取 RINEX 格式的 GPS 观测文件和导航文件，并进行伪距单点定位解算。

## 文件说明

### `spp.py`

主程序，完成完整的伪距单点定位流程，包括：

- 读取 RINEX 2.x O 文件中的观测数据
- 读取 RINEX 2.x N 文件中的广播星历
- 根据广播星历计算卫星位置
- 计算卫星钟差、相对论钟差、TGD 改正
- 进行地球自转改正
- 使用最小二乘迭代解算接收机坐标和接收机钟差
- 输出 ECEF 坐标和大地坐标

### `rinex_obs_reader_original.py`

观测文件读取示例程序，用于读取 RINEX O 文件。

主要功能：

- 读取观测文件头信息
- 解析测站名、概略坐标、观测类型等内容
- 按历元读取卫星观测值
- 支持将读取结果保存为文本或 CSV 文件

### `rinex_nav_reader_original.py`

导航文件读取示例程序，用于读取 RINEX N 文件。

主要功能：

- 读取导航文件头信息
- 解析电离层参数、UTC 参数、闰秒等内容
- 按卫星读取广播星历参数
- 支持将读取结果保存为文本或 CSV 文件

## 使用方法

运行主程序：

```powershell
python spp.py your_obs_file.12o your_nav_file.12n
```

指定某一个历元进行定位：

```powershell
python spp.py your_obs_file.12o your_nav_file.12n --epoch 0
```

从零坐标开始迭代：

```powershell
python spp.py your_obs_file.12o your_nav_file.12n --epoch 0 --initial zero
```

批量解算多个历元并保存结果：

```powershell
python spp.py your_obs_file.12o your_nav_file.12n --all --max-epochs 20 --csv result.csv
```

## 数据说明

仓库中不包含 RINEX 测试数据。使用时需要自行准备：

- RINEX O 文件，即观测文件
- RINEX N 文件，即导航电文文件

可以使用课程提供的数据，也可以使用公开 GNSS 数据进行测试。

## 定位原理简介

伪距单点定位利用多颗 GPS 卫星的位置和接收机到卫星的伪距观测值，解算接收机位置。

每个历元中，至少需要 4 颗可用 GPS 卫星，因为需要同时解算 4 个未知数：

- 接收机 X 坐标
- 接收机 Y 坐标
- 接收机 Z 坐标
- 接收机钟差

由于观测方程是非线性的，程序采用迭代最小二乘法求解。迭代收敛后，可以得到接收机的 ECEF 坐标，并进一步转换为纬度、经度和大地高。
