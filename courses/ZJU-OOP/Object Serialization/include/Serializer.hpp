#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <list>
#include <utility>
#include <set>
#include <map>
#include <type_traits>
#include <tuple>
#include <memory>
#include <sstream>
#include "tinyxml2.h" //用于解析和生成xml
#include "base64.h"   //用于base64编码

// 宏定义，需要用户将其添加到自定义的结构体中，告诉程序哪些变量需要序列化
#define SERIALIZABLE(...)                      \
    static constexpr bool serializable = true; \
    static auto get_members()                  \
    {                                          \
        return std::make_tuple(__VA_ARGS__);   \
    }

// 接口类
template <typename T>
class Serializer
{
public:
    void serialize(const T &data, const std::string &filename) // 序列化
    {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs)
        {
            throw std::runtime_error("Cannot open the file for writing: " + filename);
        }
        serialize_help(data, ofs); // 调用辅助函数，用流的方法实现序列化
    }

    void deserialize(T &data, const std::string &filename) // 反序列化
    {
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs)
        {
            throw std::runtime_error("Cannot open the file for reading: " + filename);
        }
        deserialize_help(data, ifs); // 调用辅助函数，和序列化一一对应实现反序列化
    }

private:
    virtual void serialize_help(const T &data, std::ostream &os) = 0;
    virtual void deserialize_help(T &data, std::istream &is) = 0;
};

// 前向声明Binary中的函数模板
namespace Binary
{
    template <typename T>
    void serialize_stream(const T &data, std::ostream &os);
    template <typename T>
    void deserialize_stream(T &data, std::istream &is);

    template <typename T>
    void serialize_stream(const std::vector<T> &, std::ostream &);
    template <typename T>
    void deserialize_stream(std::vector<T> &, std::istream &);

    template <typename T>
    void serialize_stream(const std::list<T> &, std::ostream &);
    template <typename T>
    void deserialize_stream(std::list<T> &, std::istream &);

    template <typename T1, typename T2>
    void serialize_stream(const std::pair<T1, T2> &, std::ostream &);
    template <typename T1, typename T2>
    void deserialize_stream(std::pair<T1, T2> &, std::istream &);

    template <typename T>
    void serialize_stream(const std::set<T> &, std::ostream &);
    template <typename T>
    void deserialize_stream(std::set<T> &, std::istream &);

    template <typename K, typename V>
    void serialize_stream(const std::map<K, V> &, std::ostream &);
    template <typename K, typename V>
    void deserialize_stream(std::map<K, V> &, std::istream &);

    template <>
    void serialize_stream(const std::string &, std::ostream &);
    template <>
    void deserialize_stream(std::string &, std::istream &);

    template <typename T>
    void serialize_stream(const std::shared_ptr<T> &, std::ostream &);
    template <typename T>
    void deserialize_stream(std::shared_ptr<T> &, std::istream &);

    template <typename T>
    void serialize_stream(const std::unique_ptr<T> &, std::ostream &);
    template <typename T>
    void deserialize_stream(std::unique_ptr<T> &, std::istream &);

} // namespace Binary

namespace user_type
{

    // 检测用户自定义的类型是否定义了serializable 静态常量，由此推出是否支持序列化
    template <typename T, typename = void>
    struct has_serializable : std::false_type
    {
    };

    template <typename T>
    struct has_serializable<T, std::void_t<decltype(T::serializable)>>
        : std::true_type
    {
    };

    template <typename T>
    constexpr bool has_serializable_v = has_serializable<T>::value;

    // 序列化辅助函数
    template <typename T, typename Tuple, std::size_t... Is>
    void serialize_members(const T &data, std::ostream &os, Tuple &&members, std::index_sequence<Is...>)
    {
        using namespace Binary;
        // 使用折叠表达式，依次序列化每个成员
        (serialize_stream(data.*std::get<Is>(members), os), ...); // data.*std::get<Is>(members) 获取对象成员变量的值
        // 通过逗号运算符，将序列化操作绑定到每个成员变量
    }

    template <typename T, typename Tuple, std::size_t... Is>
    void deserialize_members(T &data, std::istream &is, Tuple &&members, std::index_sequence<Is...>)
    {
        using namespace Binary;
        (deserialize_stream(data.*std::get<Is>(members), is), ...); // 反序列化，与序列化同理
    }

    // 用户自定义类型的序列化实现
    template <typename T>
    std::enable_if_t<user_type::has_serializable_v<T>, void>
    serialize_user(const T &data, std::ostream &os)
    {
        auto members = T::get_members(); // 获取自定义类型的成员变量元组
        constexpr auto size = std::tuple_size_v<std::decay_t<decltype(members)>>;
        user_type::serialize_members(data, os, members, std::make_index_sequence<size>{});
    }

    // 用户自定义类型反序列化的实现
    template <typename T>
    std::enable_if_t<user_type::has_serializable_v<T>, void>
    deserialize_user(T &data, std::istream &is)
    {
        auto members = T::get_members();
        constexpr auto size = std::tuple_size_v<std::decay_t<decltype(members)>>;
        user_type::deserialize_members(data, is, members, std::make_index_sequence<size>{});
    }

    // 用户自定义类型的继承类
    template <typename T>
    class BinaryUserSerializer : public Serializer<T>
    {
    private:
        void serialize_help(const T &data, std::ostream &os) override
        {
            serialize_user(data, os);
        }
        void deserialize_help(T &data, std::istream &is) override
        {
            deserialize_user(data, is);
        }
    };
} // namespace user_type

namespace Binary
{
    using namespace user_type;

    // 常用类型的继承类
    template <typename T>
    class BinarySerializer : public Serializer<T>
    {
    private:
        void serialize_help(const T &data, std::ostream &os) override
        {
            serialize_stream(data, os);
        }
        void deserialize_help(T &data, std::istream &is) override
        {

            deserialize_stream(data, is);
        }
    };

    // std::string
    template <>
    void serialize_stream(const std::string &data, std::ostream &os)
    {
        uint64_t size = data.size();
        os.write(reinterpret_cast<const char *>(&size), sizeof(size));
        os.write(data.data(), size);
    }

    template <>
    void deserialize_stream(std::string &data, std::istream &is)
    {
        uint64_t size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size));
        data.resize(size);
        is.read(&data[0], size);
    }

    //容器序列化的实现都是通过递归调用基本类型来实现
    // std::vector
    template <typename T>
    void serialize_stream(const std::vector<T> &data, std::ostream &os)
    {
        uint64_t size = data.size();
        os.write(reinterpret_cast<const char *>(&size), sizeof(size));
        for (const auto &it : data)
        {
            serialize_stream(it, os);
        }
    }

    template <typename T>
    void deserialize_stream(std::vector<T> &data, std::istream &is)
    {
        uint64_t size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size)); // 从输入流中读取向量的大小
        data.resize(size);
        for (auto &it : data)
        {
            deserialize_stream(it, is); // 递归调用 deserialize_stream 函数读取单个元素
        }
    }

    // std::list
    template <typename T>
    void serialize_stream(const std::list<T> &data, std::ostream &os)
    {
        uint64_t size = data.size();
        os.write(reinterpret_cast<const char *>(&size), sizeof(size));
        for (const auto &it : data)
        {
            serialize_stream(it, os);
        }
    }

    template <typename T>
    void deserialize_stream(std::list<T> &data, std::istream &is)
    {
        uint64_t size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size));
        data.clear();
        for (uint64_t i = 0; i < size; ++i)
        {
            T item;
            deserialize_stream(item, is);
            data.push_back(item);
        }
    }

    // std::pair
    template <typename T1, typename T2>
    void serialize_stream(const std::pair<T1, T2> &data, std::ostream &os)
    {
        serialize_stream(data.first, os);
        serialize_stream(data.second, os);
    }

    template <typename T1, typename T2>
    void deserialize_stream(std::pair<T1, T2> &data, std::istream &is)
    {
        T1 first;
        T2 second;
        deserialize_stream(first, is);
        deserialize_stream(second, is);
        data = std::make_pair(first, second);
    }

    // std::set
    template <typename T>
    void serialize_stream(const std::set<T> &data, std::ostream &os)
    {
        uint64_t size = data.size();
        os.write(reinterpret_cast<const char *>(&size), sizeof(size));
        for (const auto &it : data)
        {
            serialize_stream(it, os);
        }
    }

    template <typename T>
    void deserialize_stream(std::set<T> &data, std::istream &is)
    {
        uint64_t size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size));
        data.clear();
        for (uint64_t i = 0; i < size; ++i)
        {
            T item;
            deserialize_stream(item, is);
            data.insert(item);
        }
    }

    // std::map
    template <typename K, typename V>
    void serialize_stream(const std::map<K, V> &data, std::ostream &os)
    {
        uint64_t size = data.size();
        os.write(reinterpret_cast<const char *>(&size), sizeof(size));
        for (const auto &pair : data)
        {
            serialize_stream(pair, os);
        }
    }

    template <typename K, typename V>
    void deserialize_stream(std::map<K, V> &data, std::istream &is)
    {
        uint64_t size;
        is.read(reinterpret_cast<char *>(&size), sizeof(size));
        data.clear();
        for (uint64_t i = 0; i < size; ++i)
        {
            std::pair<K, V> pair;
            deserialize_stream(pair, is);
            data.insert(pair);
        }
    }

    // std::shared_ptr
    template <typename T>
    void serialize_stream(const std::shared_ptr<T> &ptr, std::ostream &os)
    {
        bool is_null = (ptr == nullptr);
        serialize_stream(is_null, os); // 写入 null 标志
        if (!is_null)
        {
            serialize_stream(*ptr, os); // 序列化实际对象
        }
    }

    template <typename T>
    void deserialize_stream(std::shared_ptr<T> &ptr, std::istream &is)
    {
        bool is_null;
        deserialize_stream(is_null, is);
        if (!is_null)
        {
            ptr = std::make_shared<T>();
            deserialize_stream(*ptr, is); // 反序列化实际对象
        }
        else
        {
            ptr = nullptr;
        }
    }

    // std::unique_ptr
    template <typename T>
    void serialize_stream(const std::unique_ptr<T> &ptr, std::ostream &os)
    {
        bool is_null = (ptr == nullptr);
        serialize_stream(is_null, os); // 写入 null 标志
        if (!is_null)
        {
            serialize_stream(*ptr, os); // 序列化实际对象
        }
    }

    template <typename T>
    void deserialize_stream(std::unique_ptr<T> &ptr, std::istream &is)
    {
        bool is_null;
        deserialize_stream(is_null, is);
        if (!is_null)
        {
            ptr = std::make_unique<T>();
            deserialize_stream(*ptr, is); // 反序列化实际对象
        }
        else
        {
            ptr = nullptr;
        }
    }

    // 核心实现函数，判断是自定义类型还是常用类型，分别实现序列化
    template <typename T>
    void serialize_stream(const T &data, std::ostream &os)
    {
        if constexpr (has_serializable_v<T>)
        {
            user_type::serialize_user(data, os); // 用户自定义类型的序列化实现
        }
        else
        {
            os.write(reinterpret_cast<const char *>(&data), sizeof(T)); // 基本类型的序列化实现
        }
    }

    template <typename T>
    void deserialize_stream(T &data, std::istream &is)
    {
        if constexpr (has_serializable_v<T>)
        {
            user_type::deserialize_user(data, is); // 用户自定义类型的反序列化实现
        }
        else
        {
            is.read(reinterpret_cast<char *>(&data), sizeof(T)); // 基本类型的反序列化实现
        }
    }
} // namespace Binary

namespace XML
{
    // 模式：base64编码或者文本
    enum class Mode
    {
        text,
        binary
    };

    static Mode mode = Mode::text; // 默认为文本模式
    void set_mode(Mode m)
    {
        mode = m;
    }

    // 辅助函数，用于二进制序列化并编码为 Base64
    template <typename T>
    std::string binary_to_base64(const T &value)
    {
        std::ostringstream oss(std::ios::binary);
        Binary::serialize_stream(value, oss);
        std::string data = oss.str();
        return base64_encode(data); // 使用 base64 编码
    }

    // 辅助函数，用于从 Base64 解码并二进制反序列化
    template <typename T>
    void base64_to_binary(const std::string &data, T &value)
    {
        std::string new_data = base64_decode(data);
        std::istringstream iss(new_data, std::ios::binary);
        Binary::deserialize_stream(value, iss);
    }

    // 核心序列化函数
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const T &value, const char *name);

    // 核心反序列化函数
    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, T &value, const char *name);

    // 基本类型序列化
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const T &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name); // 创建一个新的 XML 节点，名称为name
        element->SetAttribute("val", value);                     //"val = value""
        parent->InsertEndChild(element);
    }

    // 基本类型反序列化
    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, T &value, const char *name)
    {
        const tinyxml2::XMLElement *element = nullptr;
        if (!name || name[0] == '\0') // 如果为deserialize_value(parent,value,"")则使用此分支
        {
            element = parent;
        }
        else
        {
            element = parent->FirstChildElement(name);
        }
        if (!element)
            return;
        element->QueryAttribute("val", &value); //<name,val=...>
    }

    // std::string
    template <>
    void serialize_value<std::string>(tinyxml2::XMLElement *parent, const std::string &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        element->SetAttribute("val", value.c_str());
        parent->InsertEndChild(element);
    }

    template <>
    void deserialize_value<std::string>(const tinyxml2::XMLElement *parent, std::string &value, const char *name)
    {
        const tinyxml2::XMLElement *element = nullptr;
        if (!name || name[0] == '\0') // 如果为deserialize_value(parent,value,"")则使用此分支
        {
            element = parent;
        }
        else
        {
            element = parent->FirstChildElement(name);
        }
        const char *str = element->Attribute("val");
        if (str)
        {
            value = str;
        }
    }

    //以下都是通过递归调用基本类型的序列化来实现
    // std::pair
    template <typename T1, typename T2>
    void serialize_value(tinyxml2::XMLElement *parent, const std::pair<T1, T2> &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        parent->InsertEndChild(element);
        serialize_value(element, value.first, "first");
        serialize_value(element, value.second, "second");
    }

    template <typename T1, typename T2>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::pair<T1, T2> &value, const char *name)
    {
        if (auto *element = parent->FirstChildElement(name))
        {
            deserialize_value(element, value.first, "first");
            deserialize_value(element, value.second, "second");
        }
    }

    // std::vector
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const std::vector<T> &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        element->SetAttribute("size", static_cast<unsigned int>(value.size())); // 记入向量大小
        for (auto &v : value)                                                   // 递归调用
        {
            serialize_value(element, v, "item");
        }
        parent->InsertEndChild(element);
    }

    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::vector<T> &value, const char *name)
    {
        value.clear();
        if (auto *element = parent->FirstChildElement(name))
        {
            unsigned int size = 0;
            if (element->QueryAttribute("size", &size) == tinyxml2::XML_SUCCESS)
            {
                value.reserve(size);
            }
            // 遍历每个item
            for (auto *child = element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                T it;
                deserialize_value(child, it, "");
                value.push_back(std::move(it));
            }
        }
    }

    // std::list，与vector同理
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const std::list<T> &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        for (auto &v : value) // 递归调用
        {
            serialize_value(element, v, "item");
        }
        parent->InsertEndChild(element);
    }

    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::list<T> &value, const char *name)
    {
        value.clear();
        if (auto *element = parent->FirstChildElement(name))
        {
            for (auto *child = element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                T it;
                deserialize_value(child, it, "");
                value.push_back(std::move(it));
            }
        }
    }

    // std::set，与vector同理
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const std::set<T> &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        for (auto &v : value)
        {
            serialize_value(element, v, "item");
        }
        parent->InsertEndChild(element);
    }

    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::set<T> &value, const char *name)
    {
        value.clear();
        if (auto *element = parent->FirstChildElement(name))
        {
            for (auto *child = element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                T it;
                deserialize_value(child, it, "");
                value.insert(std::move(it));
            }
        }
    }

    // std::map
    template <typename K, typename V>
    void serialize_value(tinyxml2::XMLElement *parent, const std::map<K, V> &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        for (const auto &pair : value)
        {
            auto *entry = element->GetDocument()->NewElement("entry");
            serialize_value(entry, pair.first, "key");    // 递归序列化键
            serialize_value(entry, pair.second, "value"); // 递归序列化值
            element->InsertEndChild(entry);
        }
        parent->InsertEndChild(element);
    }

    template <typename K, typename V>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::map<K, V> &value, const char *name)
    {
        value.clear();
        if (auto *element = parent->FirstChildElement(name))
        {
            for (auto *child = element->FirstChildElement(); child != nullptr; child = child->NextSiblingElement())
            {
                K key;
                V val;
                deserialize_value(child, key, "key");
                deserialize_value(child, val, "value");
                value[key] = val;
            }
        }
    }

    // std::shared_ptr
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const std::shared_ptr<T> &ptr, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        parent->InsertEndChild(element);

        if (ptr == nullptr)
        {
            element->SetAttribute("nullptr", "true"); // 指针为空
        }
        else
        {
            element->SetAttribute("nullptr", "false");
            serialize_value(element, *ptr, "object"); // 序列化指向的对象
        }
    }

    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::shared_ptr<T> &ptr, const char *name)
    {
        if (auto *element = parent->FirstChildElement(name))
        {
            const char *is_null = element->Attribute("nullptr");
            if (is_null && std::string(is_null) == "true") // 指针为空
            {
                ptr = nullptr;
            }
            else
            {
                ptr = std::make_shared<T>();
                deserialize_value(element, *ptr, "object");
            }
        }
        else
        {
            ptr = nullptr;
        }
    }

    // std::unique_ptr，与std::shared_ptr同理
    template <typename T>
    void serialize_value(tinyxml2::XMLElement *parent, const std::unique_ptr<T> &ptr, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        parent->InsertEndChild(element);

        if (ptr == nullptr)
        {
            element->SetAttribute("nullptr", "true");
        }
        else
        {
            element->SetAttribute("nullptr", "false");
            serialize_value(element, *ptr, "object");
        }
    }

    template <typename T>
    void deserialize_value(const tinyxml2::XMLElement *parent, std::unique_ptr<T> &ptr, const char *name)
    {
        if (auto *element = parent->FirstChildElement(name))
        {
            const char *is_null = element->Attribute("nullptr");
            if (is_null && std::string(is_null) == "true")
            {
                ptr = nullptr;
            }
            else
            {
                ptr = std::make_unique<T>();
                deserialize_value(element, *ptr, "object");
            }
        }
        else
        {
            ptr = nullptr;
        }
    }

    // 用户自定义类型序列化
    template <typename T>
    std::enable_if_t<user_type::has_serializable_v<T>, void>
    serialize_user(tinyxml2::XMLElement *parent, const T &value, const char *name)
    {
        auto *element = parent->GetDocument()->NewElement(name);
        parent->InsertEndChild(element);

        auto members = T::get_members();
        constexpr size_t member_count = std::tuple_size_v<decltype(members)>; // 需要序列化的变量数量

        // 使用索引展开处理所有成员
        [&]<size_t... Is>(std::index_sequence<Is...>)
        {
            (serialize_value(element, value.*std::get<Is>(members),
                             ("member" + std::to_string(Is)).c_str()),
             ...); // 对每个索引Is执行一次递归的serialize_value调用
        }(std::make_index_sequence<member_count>{});
    }

    // 用户自定义类型反序列化
    template <typename T>
    std::enable_if_t<user_type::has_serializable_v<T>, void>
    deserialize_user(const tinyxml2::XMLElement *parent, T &value, const char *name)
    {
        if (auto *element = parent->FirstChildElement(name))
        {
            auto members = T::get_members();
            constexpr size_t member_count = std::tuple_size_v<decltype(members)>;

            [&]<size_t... Is>(std::index_sequence<Is...>)
            {
                (deserialize_value(element, value.*std::get<Is>(members),
                                   ("member" + std::to_string(Is)).c_str()),
                 ...);
            }(std::make_index_sequence<member_count>{});
        }
    }

    // 顶层序列化接口
    template <typename T>
    void serialize_xml(const T &data, const std::string &type, const std::string &filename)
    {
        tinyxml2::XMLDocument doc;
        auto *root = doc.NewElement("serialization");
        doc.InsertFirstChild(root);

        if (mode == Mode::binary) // base64编码模式
        {
            std::string val = binary_to_base64(data); // 直接转为base64编码
            auto *element = root->GetDocument()->NewElement(type.c_str());
            element->SetAttribute("binary", "true");
            element->SetAttribute("val", val.c_str());
            root->InsertEndChild(element);
        }
        else
        {
            if constexpr (user_type::has_serializable_v<T>) // 自定义类型
                serialize_user(root, data, type.c_str());
            else
                serialize_value(root, data, type.c_str());
        }
        if (doc.SaveFile(filename.c_str()) != tinyxml2::XML_SUCCESS) // 保存xml文件
        {
            throw std::runtime_error("Can not save the XML file: " + filename);
        }
    }

    // 顶层反序列化接口
    template <typename T>
    void deserialize_xml(T &data, const std::string &type, const std::string &filename)
    {
        tinyxml2::XMLDocument doc;
        if (doc.LoadFile(filename.c_str()) != tinyxml2::XML_SUCCESS)
        {
            throw std::runtime_error("Can not load the XML file: " + filename);
        }

        if (auto *root = doc.FirstChildElement("serialization"))
        {
            auto *element = root->FirstChildElement(type.c_str());
            const char *binary = element->Attribute("binary");
            if (binary && std::string(binary) == "true") // base64编码模式
            {
                const char *val = element->Attribute("val");
                if (val)
                {
                    base64_to_binary(std::string(val), data);
                }
            }
            else
            {
                if constexpr (user_type::has_serializable_v<T>)
                    deserialize_user(root, data, type.c_str());
                else
                    deserialize_value(root, data, type.c_str());
            }
        }
    }
} // namespace XML