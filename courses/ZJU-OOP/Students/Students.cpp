#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <climits>
using namespace std;

struct Student
{
    int no;
    string name;
    float score1;
    float score2;
    float score3;
    float average;
};

class Students
{
private:
    vector<Student> students; // 学生信息：学号，姓名，成绩1，成绩2，成绩3，平均成绩
public:
    Students() {} // 构造函数

    void addStudent(int no, string name, float score1, float score2, float score3) // 插入学生信息
    {
        float average = (score1 + score2 + score3) / 3; // 3门课程的平均分
        students.push_back({no, name, score1, score2, score3, average});
    }

    void printStudents() // 打印学生信息
    {
        cout << "no\tname\tscore1\tscore2\tscore3\taverage" << endl; // 打印列索引
        for (const auto &student : students)                         // 遍历每一行，打印学生信息
        {
            cout << student.no << "\t";
            cout << student.name << "\t";
            cout << student.score1 << "\t";
            cout << student.score2 << "\t";
            cout << student.score3 << "\t";
            cout << student.average << endl;
        }
    }

    void outputcoursestats() // 输出每门课程的平均分、最高分、最低分
    {
        float sum1 = 0, sum2 = 0, sum3 = 0;
        float min1 = INT_MAX, min2 = INT_MAX, min3 = INT_MAX;
        float max1 = INT_MIN, max2 = INT_MIN, max3 = INT_MIN;

        for (const auto &student : students) // 遍历每个学生元组
        {
            sum1 += student.score1;
            sum2 += student.score2;
            sum3 += student.score3;

            min1 = min(min1, student.score1);
            min2 = min(min2, student.score2);
            min3 = min(min3, student.score3);

            max1 = max(max1, student.score1);
            max2 = max(max2, student.score2);
            max3 = max(max3, student.score3);
        }

        int studentnum = students.size();                                                                             // 学生数量
        cout << "\taverage\t" << sum1 / studentnum << "\t" << sum2 / studentnum << "\t" << sum3 / studentnum << endl; // 每门课程的平均分
        cout << "\tmin\t" << min1 << "\t" << min2 << "\t" << min3 << endl;                                            // 每门课程的最低分
        cout << "\tmax\t" << max1 << "\t" << max2 << "\t" << max3 << endl;                                            // 每门课程的最高分
    }
};

int main()
{
    cout << "请输入10个学生的学号、姓名、成绩1、成绩2、成绩3（用空格分隔）：" << endl;
    Students students1;
    int no;
    string name;
    float score1;
    float score2;
    float score3;
    for (int count = 1; count <= 10; ++count) // 输入10个学生元组
    {
        cin >> no >> name >> score1 >> score2 >> score3;
        students1.addStudent(no, name, score1, score2, score3); // 在表中添加该学生信息
    }
    // students1.addStudent(1, "K.Weng", 5, 5, 5);
    // students1.addStudent(2, "T.Dixon", 4, 3, 2);
    // students1.addStudent(3, "V.Chu", 3, 4, 5);
    // students1.addStudent(4, "L.Tson", 4, 3, 4);
    // students1.addStudent(5, "L.Lee", 3, 4, 3);
    // students1.addStudent(6, "I.Young", 4, 2, 5);
    // students1.addStudent(7, "K.Hiro", 4, 2, 1);
    // students1.addStudent(8, "G.Ping", 4, 4, 4);
    // students1.addStudent(9, "H.Gu", 2, 3, 4);
    // students1.addStudent(10, "J.Jon", 5, 4, 3);
    students1.printStudents();     // 打印学生信息
    students1.outputcoursestats(); // 输出每门课程的平均分、最高分、最低分
}
