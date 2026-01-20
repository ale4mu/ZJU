#include "zzl.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

const double PI = 3.14;
using namespace std;

struct student
{
    int id;
    string name;
    bool operator<(const student &s) // 重载小于
    {
        return id < s.id;
    }
};

ostream &operator<<(std::ostream &out, const student &s) // 重载输出
{
    return out << "(" << s.id << "," << s.name << ")";
}

class Shape
{
protected:
    double area, perimeter;

public:
    virtual ~Shape() {};
    virtual void calc_area() = 0;      // 计算面积
    virtual void calc_perimeter() = 0; // 计算周长
    virtual string name() const = 0;   // const 表示 name() 是一个 常量成员函数（const member function），它不会修改对象的任何成员变量。
    double get_area() const { return area; };
    double get_perimeter() const { return perimeter; };
    friend ostream &operator<<(ostream &, const Shape &); // 友元函数
};

ostream &operator<<(ostream &out, const Shape &s) // 重载输出运算符
{
    return out << "(" << s.name() << ":" << s.area << "," << s.perimeter << ")";
}

class Rectangle : public Shape
{
private:
    double w, h;

public:
    Rectangle(double w, double h) : w(w), h(h) {}
    void calc_area() override
    {
        area = w * h;
    }
    void calc_perimeter() override
    {
        perimeter = 2 * (w + h);
    }
    string name() const override
    {
        return "Rectangle";
    }
};

class Circle : public Shape
{
private:
    double r;

public:
    Circle(double r) : r(r) {}
    void calc_area() override
    {
        area = PI * r * r;
    }
    void calc_perimeter() override
    {
        perimeter = 2 * PI * r;
    }
    string name() const override
    {
        return "Circle";
    }
};

class Triangle : public Shape
{
private:
    double a, b, c;

public:
    Triangle(double a, double b, double c) : a(a), b(b), c(c) {}
    void calc_area() override
    {
        double p = (a + b + c) / 2;
        area = sqrt(p * (p - a) * (p - b) * (p - c));
    }
    void calc_perimeter() override
    {
        perimeter = a + b + c;
    }
    string name() const override
    {
        return "Triangle";
    }
};

bool less_shape_area(Shape *s1, Shape *s2) // 按面积升序
{
    return s1->get_area() < s2->get_area();
}

bool less_shape_perimeter(Shape *s1, Shape *s2) // 按周长升序
{
    return s1->get_perimeter() < s2->get_perimeter();
}

bool greater_shape_area(Shape *s1, Shape *s2) // 按面积降序
{
    return s1->get_area() > s2->get_area();
}

bool greater_shape_perimeter(Shape *s1, Shape *s2) // 按周长降序
{
    return s1->get_perimeter() > s2->get_perimeter();
}

void print_sorting_option(const std::string &option) // setw(20):设置输出字段宽度为20个字符 left设置输出内容左对齐
{
    std::cout << std::setw(20) << std::left << ("by " + option + ":");
}

void print_sorting_head(const std::string &elem_type) // 说明排序的类型
{
    std::cout << "\n===== sorting [" << elem_type << "] array =====\n";
}

int main()
{
    // 测试整数
    {
        print_sorting_head("int");
        vector<int> arr = {64, 25, 12, 23, 11};
        zzl::selection_sort(arr);
        zzl::print_array(arr);
    }
    // 测试double
    {
        print_sorting_head("double");
        vector<double> arr = {64.1, 25.2, 12.3, 23.4, 11.5};
        zzl::selection_sort(arr);
        zzl::print_array(arr);
    }
    // 测试字符串
    {
        print_sorting_head("string");
        vector<string> arr = {"hello", "world", "zju", "boys", "girls"};
        zzl::selection_sort(arr);
        zzl::print_array(arr);
    }
    // 测试结构体
    {
        print_sorting_head("Student");
        vector<student> arr = {{2011, "Newton"}, {2002, "Gauss"}, {2133, "Euler"}, {2076, "Riemann"}, {2054, "Mozi"}};
        print_sorting_option("id"); // 默认按id排序
        zzl::selection_sort(arr);
        zzl::print_array(arr);
        print_sorting_option("name"); // 按姓名排序
        zzl::selection_sort(arr, [](const student &a, const student &b)
                            { return a.name < b.name; });
        zzl::print_array(arr);
    }

    {
        print_sorting_head("Shape");
        vector<Shape *> arr = {new Rectangle(2, 3), new Rectangle(5, 5), new Circle(3), new Triangle(2, 5, 4)};
        for (Shape *s : arr)
        {
            s->calc_area();
            s->calc_perimeter();
        }
        print_sorting_option("area<");
        zzl::selection_sort(arr, less_shape_area);
        zzl::print_array(arr);
        print_sorting_option("perimeter<");
        zzl::selection_sort(arr, less_shape_perimeter);
        zzl::print_array(arr);
        print_sorting_option("area>");
        zzl::selection_sort(arr, greater_shape_area);
        zzl::print_array(arr);
        print_sorting_option("perimeter>");
        zzl::selection_sort(arr, greater_shape_perimeter);
        zzl::print_array(arr);
        // 防止内存泄漏
        for (Shape *s : arr)
            delete s;
    }
}
