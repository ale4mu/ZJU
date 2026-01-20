### OOP Project：OBJECT SERIALIZATION

#### 1.简介
借助tinyxml2和base64库，实现常见类型和用户自定义类型的序列化与反序列化。
- [√] Implement a module to support binary serialization/deserialization
- [√] Implement a wrapper module of tinyxml2 to support XML serialization
- [√] Every module has its namespace
- [√] Both modules should at least support the serialization of arithmetic types, C++ string type, and STL containers
- [√] Both modules should provide a convenient mechanism to support the serialization of user-defined types
- [√] During testing, you should cover all the required types and options
- [√] Use binary-to-text encoding/decoding (base64) to implement a binary mode of XML serialization
- [√] Support the serialization of smart pointers
  
#### 2.总仓库结构
总仓库结构如下：
```bash

├── OBJECT SERIALIZATION
│   ├── README.md 		
│   ├── Makefile 		
│   ├── data 			
│   │   └──...
│   ├── src 			
│   │   └── test.cpp
│   ├── include		
│   │   └── Serializer.hpp
│   │   └── tinyxml2.h
│   │   └── base64.h
│   │   └── tinyxml2.cpp
│   │   └── base64.cpp
│   │   └──...
└── ...

```

#### 3.使用方法
- 进入OBJECT SERIALIZATION文件夹目录下，打开终端编译：```make```
- 输出测试样例：```./test```
- 查看生成的二进制文件和XML文件： ``` ls data```
- 删除生成文件：```make clean```