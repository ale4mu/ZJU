### 项目作业——MNIST&神经网络推理

#### 1. 项目简介
通过 Python 导出 CNN 的网络参数至二进制文件，接着通过 C 语言搭建可用于正向推理的神经网络，对自己手写数字进行识别，输出识别结果。

#### 2. 总仓库结构
总仓库结构如下：
```bash
.
├── README.md
├── .gitignore
├── ass01
│   └── ...
├── mnist
│   ├── README.md 		
│   ├── Makefile 		
│   ├── images 			
│   │   └── test_0.bmp
│   │   └── ...
│   ├── src 			
│   │   └── ...
│   ├── ExportPara 		
│   │   └── export_cnn.ipynb
│   │   └── parameters_cnn.bin
│   └── data
│       └── t10k-images-idx3-ubyte
│       └── t10k-labels-idx1-ubyte
└── ...
```
#### 3.使用方法
- 进入mnist文件夹目录下，打开终端编译：```make```
- 测试特定数字如数字1的手写图片：```./test images/test_1.bmp```
- 测试0-9手写数字：```make run```
- 输出数据集data/t10k-images-idx3-ubyte的识别准确率：```./test```
- 删除生成文件：```make clean```

#### 4.完成内容
- [√] 实现基本要求
- [√] 调整网络和参数，提升c程序的推理效率
- 针对第二项，在 python 构建的CNN模型中减少了卷积层通道数（1 -> 8 -> 16 ）和全连接层参数（784 -> 32 -> 10），同时对权重参数进行了剪枝，可在C代码中跳过零值计算，加速推理。


