#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <stdexcept>
using namespace std;

int Zeller(int year, int month, int day) // 根据Zeller公式计算某年某月某天是星期几
{
    if (month == 1 || month == 2) // 1月和2月被认为是上一年的13月和14月
    {
        month += 12;
        --year;
    }
    return (year % 100 + (year % 100) / 4 + year / 400 - 2 * (year / 100) + 26 * (month + 1) / 10 + day - 1) % 7; // Zeller公式
}

vector<vector<string>> Month(int year, int month) // 将一个月的日历存储在一个二维数组中
{
    vector<vector<string>> singleMonth(7, vector<string>(21, " ")); // 将二维数组初始化
    const vector<string> monthName = {"January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
    vector<int> monthDay = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; // 每个月份的天数
    if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0)               // 闰年2月29天
        monthDay[1] = 29;
    int numSpaces = (21 - monthName[month - 1].size()) / 2; // 计算月份名字前的空格数
    for (char c : monthName[month - 1])
    {
        singleMonth[0][numSpaces++] = c; // 存储月份名字
    }
    int firstDay = Zeller(year, month, 1); // 计算这个月的第一天是星期几
    int rows = 1;
    int index = 3 * firstDay; // 第一天的位置
    int i = 1;
    while (i <= monthDay[month - 1])
    {
        if (i < 10) // 个位数
        {
            singleMonth[rows][++index] = to_string(i);
            index += 2;
        }
        else // 两位数
        {
            singleMonth[rows][index++] = to_string(i / 10);
            singleMonth[rows][index] = to_string(i % 10);
            index += 2;
        }
        if (index >= 21) // 换行
        {
            ++rows;
            index = 0;
        }
        ++i;
    }
    return singleMonth;
}

void PrintnMonth(int year, int beginMonth, int n) // 一排打印n个月
{
    vector<vector<vector<string>>> row;
    for (int i = 0; i < n; ++i)
    {
        if (beginMonth + i > 12) // 超出12个月跳出循环
            break;
        row.push_back(Month(year, beginMonth + i)); // 生成这一排月份的二维数组
    }
    for (int i = 0; i < 7; ++i) // 每个月共7行
    {
        for (auto month : row)
        {
            for (auto s : month[i])
            {
                cout << s;
            }
            cout << " ";
        }
        cout << endl;
    }
}

void PrintCalendar(int year, int n) // 打印12个月
{
    if (n <= 0 || n > 12)
    {
        throw invalid_argument("The number of months you want to print in one row should be in the range of 1 to 12.");
    }
    else if (year <= 1582 || year > 9999) // Zeller适用于1583年以后
    {
        throw invalid_argument("The year should be in the range of 1583 to 9999.");
    }
    for (int i = 1; i <= 12 / n + 1; ++i)
    {
        PrintnMonth(year, (i - 1) * n + 1, n);
    }
}

int main()
{
    int y;
    int n;
    cout << "Please input the year:";
    cin >> y;
    cout << "Please input the number of months you want to print in one row:";
    cin >> n;
    try
    {
        PrintCalendar(y, n);
    }
    catch (const invalid_argument &e)
    {
        cerr << "Error: " << e.what() << endl; // 抛出异常
    }
}