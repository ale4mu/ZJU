#include "Diary.h"
#include <fstream>

DiaryBook::DiaryBook(const std::string &_name) : name(_name) {}

DiaryBook::~DiaryBook()
{
    for (const auto &pair : diarys)
    {
        delete pair.second;
    }
}

bool DiaryBook::load()
{
    std::ifstream file(name);
    if (!file.is_open())
    {
        return false;
    }
    std::string date, content, line;
    while (std::getline(file, date))
    {
        if (!is_valid_date(date)) // 处理异常日期
        {
            return false;
        }
        content.clear();
        while (std::getline(file, line) && !line.empty()) // 读取日记直到遇到空行，空行是分隔
        {
            if (!content.empty())
            {
                content += "\n";
            }
            content += line;
        }

        diarys[date] = new Diary{date, content};
        all_dates.push_back(date);
    }
    file.close();
    return true;
}

bool DiaryBook::save() // 保存日记
{
    std::ofstream file(name);
    if (!file.is_open())
    {
        return false;
    }
    for (const auto &pair : diarys) // pair.first是日期，pair.second是指向日记结构体的指针
    {
        file << pair.second->date << std::endl;
        file << pair.second->content << std::endl;
        file << std::endl; // 用空行分隔
    }
    file.close();
    return true;
}

void DiaryBook::add_diary(const std::string &date, const std::string &content)
{
    Diary *diary = new Diary{date, content};
    auto it = diarys.find(date);
    if (it != diarys.end())
    {
        delete it->second;  // 删除原来的日记
        it->second = diary; // 更新为新的日记
    }
    else
    {
        diarys[date] = diary;
        all_dates.push_back(date);
    }
}

void DiaryBook::date_sort() // 排序
{
    std::sort(all_dates.begin(), all_dates.end());
}

void DiaryBook::list_diary(const std::string &start_date, const std::string &end_date)
{
    date_sort();
    auto start_it = std::lower_bound(all_dates.begin(), all_dates.end(), start_date); // 大于等于起始日期
    auto end_it = std::upper_bound(all_dates.begin(), all_dates.end(), end_date);     // 小于等于结束日期
    for (auto it = start_it; it != end_it; ++it)
    {
        std::string current_date = *it;
        std::cout << current_date << std::endl;
        std::cout << diarys[current_date]->content << std::endl;
        std::cout << std::endl; // 空行分隔
    }
}

void DiaryBook::list_diary()
{
    date_sort();
    for (auto it = all_dates.begin(); it != all_dates.end(); ++it)
    {
        std::string current_date = *it;
        std::cout << current_date << std::endl;
        std::cout << diarys[current_date]->content << std::endl;
        std::cout << std::endl; // 空行分隔
    }
}

void DiaryBook::get_diary(const std::string &date) // 输出指定日期的日记
{
    std::cout << date << std::endl;
    std::cout << diarys[date]->content << std::endl;
}

bool DiaryBook::remove_diary(const std::string &date) // 删除指定日期的日记
{
    auto it = diarys.find(date);
    if (it != diarys.end())
    {
        diarys.erase(it);
        auto date_it = std::find(all_dates.begin(), all_dates.end(), date);
        all_dates.erase(date_it);
    }
    else
    {
        std::cout << "The date does not exist!" << std::endl;
        return false;
    }
    return true;
}

bool is_valid_date(const std::string &date) // 日期格式为 YYYY-MM-DD
{

    if (date.size() != 10 || date[4] != '-' || date[7] != '-')
    {
        return false;
    }
    for (size_t i = 0; i < date.size(); ++i)
    {
        if (i == 4 || i == 7)
            continue;
        if (!isdigit(date[i]))
        {
            return false;
        }
    }
    return true;
}