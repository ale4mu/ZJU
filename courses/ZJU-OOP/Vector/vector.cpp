#include <iostream>
#include <stdexcept>

template <class T>
class Vector
{
public:
    Vector();                   // creates an empty vector
    Vector(int size);           // creates a vector for holding 'size' elements
    Vector(const Vector &r);    // the copy ctor
    ~Vector();                  // destructs the vector
    T &operator[](int index);   // accesses the specified element without bounds checking
    T &at(int index);           // accesses the specified element, throws an exception of
                                // type 'std::out_of_range' when index <0 or >=m_nSize
    int size() const;           // return the size of the container
    void push_back(const T &x); // adds an element to the end
    void clear();               // clears the contents
    bool empty() const;         // checks whether the container is empty
private:
    void inflate();  // expand the storage of the container to a new capacity,
                     // e.g. 2*m_nCapacity
    T *m_pElements;  // pointer to the dynamically allocated storage
    int m_nSize;     // the number of elements in the container
    int m_nCapacity; // the total number of elements that can be held in the
                     // allocated storage
};

template <class T>
Vector<T>::Vector() // 默认构造函数
{
    m_nSize = 0;
    m_nCapacity = 0;
    m_pElements = nullptr;
}

template <class T>
Vector<T>::Vector(int size) // 构造函数
{
    m_nSize = 0;
    m_nCapacity = size;
    m_pElements = new T[m_nCapacity];
}

template <class T>
Vector<T>::Vector(const Vector &r) // 拷贝构造函数
{
    if (this == &r)
        return;
    m_nSize = r.m_nSize;
    m_nCapacity = r.m_nCapacity;
    m_pElements = new T[m_nCapacity];
    for (int i = 0; i < m_nSize; ++i)
    {
        m_pElements[i] = r.m_pElements[i];
    }
}

template <class T>
Vector<T>::~Vector() // 析构函数
{
    clear();
}

template <class T>
T &Vector<T>::operator[](int index) // 下标运算符，不考虑越界
{
    return m_pElements[index];
}

template <class T>
T &Vector<T>::at(int index) // 考虑越界
{
    if (index < 0 || index >= m_nSize)
    {
        throw std::out_of_range("index out of range");
    }
    return m_pElements[index];
}

template <class T>
int Vector<T>::size() const // 返回元素个数
{
    return m_nSize;
}

template <class T>
void Vector<T>::push_back(const T &x) // 在末尾添加元素
{
    if (m_nSize == m_nCapacity) // 扩容
    {
        inflate();
    }
    m_pElements[m_nSize] = x;
    ++m_nSize;
}

template <class T>
void Vector<T>::clear() // 清空
{
    delete[] m_pElements;
    m_nSize = 0;
    m_nCapacity = 0;
    m_pElements = nullptr;
}

template <class T>
bool Vector<T>::empty() const // 是否为空
{
    return m_nSize == 0;
}

template <class T>
void Vector<T>::inflate() // 扩容成原先容量的两倍
{
    int new_m_nCapacity = m_nCapacity;
    if (m_nCapacity == 0)
        new_m_nCapacity = 1;
    else
        new_m_nCapacity = m_nCapacity * 2;
    T *new_m_pElements = new T[new_m_nCapacity];
    for (int i = 0; i < m_nSize; ++i)
    {
        new_m_pElements[i] = m_pElements[i];
    }
    delete[] m_pElements;
    m_pElements = new_m_pElements;
    m_nCapacity = new_m_nCapacity;
}

int main()
{
    Vector<int> vec_int(5);     // 测试构造
    for (int i = 0; i < 5; ++i) // 测试push_back
    {
        vec_int.push_back(i);
    }
    std::cout << "size of vec_int: " << vec_int.size() << std::endl; // 测试size
    for (int i = 0; i < vec_int.size(); ++i)
    {
        std::cout << "vec_int[" << i << "] = " << vec_int[i] << std::endl; // 测试[]
    }
    std::cout << "test at() function: vec_int[2] = " << vec_int.at(2) << std::endl; // 测试at()
    vec_int.push_back(5);                                                           // 测试扩容
    std::cout << "size of vec_int after inflate: " << vec_int.size() << std::endl;
    std::cout << "vec_int[5] = " << vec_int.at(5) << std::endl;

    try
    {
        vec_int.at(10); // 测试越界
    }
    catch (const std::out_of_range &e)
    {
        std::cout << "test vec_int.at(6) which is out of range :" << e.what() << std::endl;
    }

    Vector<int> vec_copy(vec_int); //测试拷贝构造
    for(int i = 0; i < vec_copy.size(); ++i)
    {
        std::cout << "vec_copy[" << i << "] = " << vec_copy[i] << std::endl; 
    }
    vec_copy.push_back(100); //测试深拷贝
    std::cout << "size of vec_copy: " << vec_copy.size() << std::endl;
    std::cout << "size of vec_int after changing vec_copy: " << vec_int.size() << std::endl;

    vec_int.clear(); // 测试clear()
    std::cout << "size of vec_int after clear: " << vec_int.size() << std::endl;
    if (vec_int.empty()) // 测试empty()
        std::cout << "vec_int is empty" << std::endl;
}