#include "Serializer.hpp"


struct UserDefinedType
{
    int idx;
    std::string name;
    std::vector<double> data;

    SERIALIZABLE(&UserDefinedType::idx,
                 &UserDefinedType::name,
                 &UserDefinedType::data);
};

struct Person
{
    SERIALIZABLE(&Person::name, &Person::age)
    std::string name;
    int age;
};
int main()
{
    //测试思路是先将数据序列化，再反序列化，如果能正确接受说明测试通过
    //int
    std::cout << "Tesing Binary: int..." << std::endl;
    int n0 = 33, n1;
    Binary::BinarySerializer<int> test1;
    test1.serialize(n0, "data/n.data");
    test1.deserialize(n1, "data/n.data");
    std::cout <<  n1 << std::endl;
    std::cout << "------------------------" << std::endl;

    //string
    std::cout << "Tesing Binary: string..." << std::endl;
    std::string s1 = "hello", s2;
    Binary::BinarySerializer<std::string> test2;
    test2.serialize(s1, "data/s.data");
    test2.deserialize(s2, "data/s.data");
    std::cout << s2 << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector<int>
    std::cout << "Tesing Binary: vector<int>..." << std::endl;
    std::vector<int> v1 = {1, 2, 3, 4, 5}, v2;
    Binary::BinarySerializer<std::vector<int>> test3;
    test3.serialize(v1, "data/v.data");
    test3.deserialize(v2, "data/v.data");
    for (const auto &it : v2)
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector<string>
    std::cout << "Tesing Binary:vector<string>..." << std::endl;
    std::vector<std::string> vv1 = {"I", "LOVE", "U"}, vv2;
    Binary::BinarySerializer<std::vector<std::string>> test33;
    test33.serialize(vv1, "data/vv.data");
    test33.deserialize(vv2, "data/vv.data");
    for (const auto &it : vv2)
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //list<int>
    std::cout << "Tesing Binary: list<int>..." << std::endl;
    std::list<int> l1 = {1, 2, 3, 4, 5}, l2;
    Binary::BinarySerializer<std::list<int>> test4;
    test4.serialize(l1, "data/l.data");
    test4.deserialize(l2, "data/l.data");
    for (const auto &it : l2)
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //std::pair<int,int>
    std::cout << "Tesing Binary: pair<int,int>..." << std::endl;
    std::pair<int, int> p1 = {1, 2}, p2;
    Binary::BinarySerializer<std::pair<int, int>> test5;
    test5.serialize(p1, "data/p.data");
    test5.deserialize(p2, "data/p.data");
    std::cout << p2.first << " " << p2.second << std::endl;
    std::cout << "------------------------" << std::endl;

    //set<int>
    std::cout << "Tesing Binary: set<int>..." << std::endl;
    std::set<int> set1 = {1, 2, 3, 4, 5}, set2;
    Binary::BinarySerializer<std::set<int>> test6;
    test6.serialize(set1, "data/set.data");
    test6.deserialize(set2, "data/set.data");
    for (const auto &it : set2)
    {
        std::cout << it << " ";
    }
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //map<int,int>
    std::cout << "Tesing Binary: map<int,int>..." << std::endl;
    std::map<int, int> m1, m2;
    m1[1] = 2;
    m1[2] = 4;
    m1[3] = 8;
    Binary::BinarySerializer<std::map<int, int>> test7;
    test7.serialize(m1, "data/m.data");
    test7.deserialize(m2, "data/m.data");
    for (const auto &it : m2)
    {
        std::cout << it.first << " " << it.second << std::endl;
    }
    std::cout << "------------------------" << std::endl;

    //用户自定义类型
    std::cout << "Tesing Binary: UserDefineType..." << std::endl;
    UserDefinedType data = {42, "zhang", {1.1, 2.2, 3.3}};
    Binary::BinarySerializer<UserDefinedType> test8;
    test8.serialize(data, "data/user_data.data");
    UserDefinedType loaded;
    test8.deserialize(loaded, "data/user_data.data");
    std::cout << "Loaded: " << loaded.idx << ", " << loaded.name << "," << loaded.data[0] << " " << loaded.data[1] << " " << loaded.data[2] << std::endl;
    std::cout << "------------------------" << std::endl;

    //shared_ptr<int>
    std::cout << "Tesing Binary: share_ptr<int>..." << std::endl;
    auto ptr1 = std::make_shared<int>(1);
    Binary::BinarySerializer<std::shared_ptr<int>> ptr;
    ptr.serialize(ptr1, "data/ptr.data");
    std::shared_ptr<int> ptr2;
    ptr.deserialize(ptr2, "data/ptr.data");
    std::cout << *ptr2 << std::endl;
    std::cout << "------------------------" << std::endl;

    //unique_ptr<int>
    std::cout << "Tesing Binary: unique_ptr<int>..." << std::endl;
    auto ptr3 = std::make_unique<int>(1);
    Binary::BinarySerializer<std::unique_ptr<int>> _ptr;
    _ptr.serialize(ptr3, "data/ptr_unique.data");
    std::unique_ptr<int> ptr4;
    _ptr.deserialize(ptr4, "data/ptr_unique.data");
    std::cout << *ptr4 << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector<person>
    std::cout << "Tesing Binary: vector<Person>..." << std::endl;
    std::vector<Person> people = {{"Alice", 30}, {"Bob", 25}};
    Binary::BinarySerializer<std::vector<Person>> test9;
    test9.serialize(people, "data/people.data");
    std::vector<Person> people2;
    test9.deserialize(people2, "data/people.data");
    std::cout << people[0].name << " " <<  people[0].age << std::endl;
    std::cout << people[1].name << " " <<  people[1].age << std::endl;
    std::cout << "------------------------" << std::endl;



    //xml(text)
    //int
    std::cout << "Tesing XML(text): int..." << std::endl;
    int n =1,m;
    XML::serialize_xml(n, "int", "data/int.xml");
    XML::deserialize_xml(m, "int", "data/int.xml");
    std::cout << m << std::endl;
    std::cout << "------------------------" << std::endl;

    //string
    std::cout << "Tesing XML(text): string..." << std::endl;
    std::string str1 ="hello world",str2;
    XML::serialize_xml(str1, "string", "data/string.xml");
    XML::deserialize_xml(str2, "string", "data/string.xml");
    std::cout << str2 << std::endl;
    std::cout << "------------------------" << std::endl;

    //pair
    std::cout << "Tesing XML(text): pair..." << std::endl;   
    std::pair<int, double> pair0 = {2, 3.1}, pair1;
    XML::serialize_xml(pair0, "std_pair", "data/pair.xml");
    XML::deserialize_xml(pair1, "std_pair", "data/pair.xml");
    std::cout << pair1.first << " " << pair1.second << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector
    std::cout << "Tesing XML(text): vector<int>..." << std::endl;
    std::vector<int> vec1={1,2,3},vec2;
    XML::serialize_xml(vec1, "std_vector", "data/vector.xml");
    XML::deserialize_xml(vec2, "std_vector", "data/vector.xml");
    for(auto it:vec2)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector
    std::cout << "Tesing XML(text): vector<string>..." << std::endl;
    std::vector<std::string> vec_s1={"zhang","chen"},vec_s2;
    XML::serialize_xml(vec_s1, "std_vector", "data/vector_string.xml");
    XML::deserialize_xml(vec_s2, "std_vector", "data/vector_string.xml");
    std::cout<< vec_s2[0] << " " << vec_s2[1] << std::endl;
    std::cout << "------------------------" << std::endl;

    //list
    std::cout << "Tesing XML(text): list<int>..." << std::endl;
    std::list<int> list1={1,2,3},list2;
    XML::serialize_xml(list1, "std_list", "data/list.xml");
    XML::deserialize_xml(list2, "std_list", "data/list.xml");
    for(auto it:list2)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //set
    std::cout << "Tesing XML(text): set<int>..." << std::endl;
    std::set<int> set3={1,2,3},set4;
    XML::serialize_xml(set3, "std_set", "data/set.xml");
    XML::deserialize_xml(set4, "std_set", "data/set.xml");
    for(auto it:set4)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;
    
    //map
    std::cout << "Tesing XML(text): map<int,int>..." << std::endl;
    std::map<int,int> map1={{1,2},{3,4}},map2;
    XML::serialize_xml(map1, "std_map", "data/map.xml");
    XML::deserialize_xml(map2, "std_map", "data/map.xml");
    for(auto it:map2)
        std::cout<< it.first << " "<< it.second << std::endl;
    std::cout << "------------------------" << std::endl;

    //Person
    std::cout << "Tesing XML(text): Person..." << std::endl;
    Person person1{"Zhang", 18},person2;
    XML::serialize_xml(person1,"Person","data/person.xml");
    XML::deserialize_xml(person2,"Person","data/person.xml");
    std::cout << person2.name << " " << person2.age << std::endl;
    std::cout << "------------------------" << std::endl;

    //shared_ptr
    std::cout << "Tesing XML(text): shared_ptr..." << std::endl;
    auto ptr111 = std::make_shared<int>(1);
    std::shared_ptr<int> ptr222;
    XML::serialize_xml(ptr111,"std_shared_ptr","data/share_ptr.xml");
    XML::deserialize_xml(ptr222,"std_shared_ptr","data/share_ptr.xml");
    std::cout << *ptr222 << std::endl;
    std::cout << "------------------------" << std::endl;

    //unique_ptr
    std::cout << "Tesing XML(text): unique_ptr..." << std::endl;
    auto ptr333 = std::make_unique<int>(1);
    std::unique_ptr<int> ptr444;
    XML::serialize_xml(ptr333,"std_unique_ptr","data/unique_ptr.xml");
    XML::deserialize_xml(ptr444,"std_unique_ptr","data/unique_ptr.xml");
    std::cout << *ptr444 << std::endl;
    std::cout << "------------------------" << std::endl;

///////---------------------

    //xml(base64)
    XML::set_mode(XML::Mode::binary);
    //int
    std::cout << "Tesing XML(text): int..." << std::endl;
    int a=3,b;
    XML::serialize_xml(a, "int", "data/int_base64.xml");
    XML::deserialize_xml(b, "int", "data/int_base64.xml");
    std::cout << b << std::endl;
    std::cout << "------------------------" << std::endl;

    //string
    std::cout << "Tesing XML(base64): string..." << std::endl;
    std::string str11 ="hello world !",str22;
    XML::serialize_xml(str11, "string", "data/string_base64.xml");
    XML::deserialize_xml(str22, "string", "data/string_base64.xml");
    std::cout << str22 << std::endl;
    std::cout << "------------------------" << std::endl;

    //pair
    std::cout << "Tesing XML(base64): pair..." << std::endl;   
    std::pair<int, double> pair00 = {2, 3.1}, pair11;
    XML::serialize_xml(pair00, "std_pair", "data/pair_base64.xml");
    XML::deserialize_xml(pair11, "std_pair", "data/pair_base64.xml");
    std::cout << pair11.first << " " << pair11.second << std::endl;
    std::cout << "------------------------" << std::endl;

    //vector
    std::cout << "Tesing XML(base64): vector<int>..." << std::endl;
    std::vector<int> vec11={1,2,3},vec22;
    XML::serialize_xml(vec11, "std_vector", "data/vector_base64.xml");
    XML::deserialize_xml(vec22, "std_vector", "data/vector_base64.xml");
    for(auto it:vec22)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //list
    std::cout << "Tesing XML(base64): list<int>..." << std::endl;
    std::list<int> list11={1,2,3},list22;
    XML::serialize_xml(list11, "std_list", "data/list_base64.xml");
    XML::deserialize_xml(list22, "std_list", "data/list_base64.xml");
    for(auto it:list22)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;

    //set
    std::cout << "Tesing XML(base 64): set<int>..." << std::endl;
    std::set<int> set33={1,2,3},set44;
    XML::serialize_xml(set33, "std_set", "data/set_base64.xml");
    XML::deserialize_xml(set44, "std_set", "data/set_base64.xml");
    for(auto it:set44)
        std::cout<< it << " ";
    std::cout << std::endl;
    std::cout << "------------------------" << std::endl;
    
    //map
    std::cout << "Tesing XML(base64): map<int,int>..." << std::endl;
    std::map<int,int> map11={{1,2},{3,4}},map22;
    XML::serialize_xml(map11, "std_map", "data/map_base64.xml");
    XML::deserialize_xml(map22, "std_map", "data/map_base64.xml");
    for(auto it:map22)
        std::cout<< it.first << " "<< it.second << std::endl;
    std::cout << "------------------------" << std::endl;

    //Person
    std::cout << "Tesing XML(base64): Person..." << std::endl;
    Person person11{"Zhang", 18},person22;
    XML::serialize_xml(person11,"Person","data/person_base64.xml");
    XML::deserialize_xml(person22,"Person","data/person_base64.xml");
    std::cout << person22.name << " " << person22.age << std::endl;
    std::cout << "------------------------" << std::endl;

    //shared_ptr
    std::cout << "Tesing XML(base64): shared_ptr..." << std::endl;
    auto ptr11 = std::make_shared<int>(1);
    std::shared_ptr<int> ptr22;
    XML::serialize_xml(ptr11,"std_shared_ptr","data/share_ptr_base64.xml");
    XML::deserialize_xml(ptr22,"std_shared_ptr","data/share_ptr_base64.xml");
    std::cout << *ptr22 << std::endl;
    std::cout << "------------------------" << std::endl;

    //unique_ptr
    std::cout << "Tesing XML(base64): unique_ptr..." << std::endl;
    auto ptr33 = std::make_unique<int>(1);
    std::unique_ptr<int> ptr44;
    XML::serialize_xml(ptr33,"std_unique_ptr","data/unique_ptr_base64.xml");
    XML::deserialize_xml(ptr44,"std_unique_ptr","data/unique_ptr_base64.xml");
    std::cout << *ptr44 << std::endl;
    std::cout << "------------------------" << std::endl;
}