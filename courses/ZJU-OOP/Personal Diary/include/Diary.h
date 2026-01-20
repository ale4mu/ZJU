#ifndef Diary_H
#define Diary_H

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <iostream>

// 每篇日记的结构
struct Diary
{
    std::string date;    // 日期
    std::string content; // 内容
};

// 对日记本进行管理
class DiaryBook
{

private:
    std::string name;                                // 日记本的名字
    std::unordered_map<std::string, Diary *> diarys; // 存储所有日记指针的哈希表，方便快速查找
    std::vector<std::string> all_dates;              // 存储所有的日期

public:
    DiaryBook(const std::string &name); // 构造函数

    ~DiaryBook();

    bool load(); // 加载日记本

    bool save(); // 保存日记本

    void add_diary(const std::string &date, const std::string &content);

    void date_sort();

    void list_diary(const std::string &start_date, const std::string &end_date); // 列出指定日期范围内的日记
    void list_diary();

    void get_diary(const std::string &date); // 获取指定日期的日记

    bool remove_diary(const std::string &date); // 删除指定日期的日记
};

bool is_valid_date(const std::string &date); // 判断日期是否合法

#endif
