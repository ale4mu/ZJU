#pragma once
#include <iostream>
#include <vector>
#include <algorithm>

namespace zzl
{

    template <typename T, typename Compare = std::less<>>             // 可自定义排序，默认排序为从小到大
    void selection_sort(std::vector<T> &arr, Compare cmp = Compare()) // 排序
    {
        int n = arr.size();
        for (int i = 0; i < n - 1; ++i)
        {
            int min_index = i;
            for (int j = i + 1; j < n; ++j)
            {
                if (cmp(arr[j], arr[min_index]))
                {
                    min_index = j;
                }
            }
            if (min_index != i)
                std::swap(arr[i], arr[min_index]);
        }
    }

    template <typename T>
    void print_array(const std::vector<T> &arr) // 打印
    {
        for (auto member : arr)
        {
            std::cout << member << ' ';
        }
        std::cout << '\n';
    }

    template <typename T>
    void print_array(const std::vector<T *> &arr) // 打印
    {
        for (auto member : arr)
        {
            std::cout << *member << ' ';
        }
        std::cout << '\n';
    }
}