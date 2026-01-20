#ifndef FUNCTION_HPP
#define FUNCTION_HPP

// Function.hpp
template <typename T> 
class Function
{
public:
    virtual T operator()(const T &x) const = 0; // 用来计算函数值
    virtual ~Function() = default;             
};

#endif // FUNCTION_HPP
