#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <random>
#include <algorithm>
#include <memory>
#include <ctime>

using namespace std;

class Room
{
private:
    string name;                         // 房间的名字
    bool monster = false;                // 默认无monster
    bool princess = false;               // 默认无princess
    unordered_map<string, string> exits; // 哈希表存房间的出口

public:
    Room() {}

    Room(string _name) : name(_name) {} // 构造函数

    ~Room() {}

    void description() // 描述房间信息
    {
        cout << "Welcome to the " << name << " . There are " << exits.size() << " exits: ";
        for (auto it = exits.begin(); it != exits.end(); ++it)
        {
            cout << it->first;
            if (next(it) != exits.end())
            {
                cout << ", ";
            }
        }
        cout << "." << endl;
    }

    void addExit(const string &direction, const string &room_name) // 为房间增加出口
    {
        exits[direction] = room_name;
    }

    bool has_monster() const
    {
        return monster;
    }

    bool has_princess() const
    {
        return princess;
    }

    void set_monster() // 设置monster
    {
        monster = true;
    }

    void set_princess() // 设置princess
    {
        princess = true;
    }

    const string &get_name() const
    {
        return name;
    }

    string get_exit(const string &direction) // 返回某个出口对应的房间
    {
        if (exits.find(direction) == exits.end())
        {
            return "";
        }
        return exits[direction];
    }
};

class Game
{
private:
    vector<unique_ptr<Room>> rooms;        // 房间的信息
    unordered_map<string, int> room_index; // 房间名称到索引的映射
    bool fixed_map = false;                // 是否固定地图
public:
    void set_fixed_map() // 固定地图
    {
        fixed_map = true;
    }

    bool is_fixed() const
    {
        return fixed_map;
    }

    void set_game() // 准备游戏，连接各个房间，随机生成monster和princess
    {
        vector<string> room_names = {"Lobby", "Kitchen", "Bedroom", "Living Room", "Bathroom", "Hallway", "Basement", "Toilet", "Study", "Attic", "Balcony", "Laundry"};
        for (const auto &name : room_names)
        {
            rooms.emplace_back(make_unique<Room>(name));
        }
        for (int i = 0; i < rooms.size(); ++i)
        {
            room_index[rooms[i]->get_name()] = i;
        }

        rooms[0]->addExit("north", "Kitchen");
        rooms[0]->addExit("east", "Bedroom");
        rooms[0]->addExit("south", "Living Room");
        rooms[0]->addExit("west", "Bathroom");
        rooms[0]->addExit("up", "Hallway");
        rooms[0]->addExit("down", "Basement");
        rooms[1]->addExit("south", "Lobby");
        rooms[2]->addExit("west", "Lobby");
        rooms[3]->addExit("north", "Lobby");
        rooms[4]->addExit("east", "Lobby");
        rooms[5]->addExit("down", "Lobby");
        rooms[5]->addExit("east", "Toilet");
        rooms[6]->addExit("up", "Lobby");
        rooms[7]->addExit("west", "Hallway");
        rooms[7]->addExit("east", "Laundry");
        rooms[11]->addExit("west", "Toilet");
        rooms[5]->addExit("north", "Study");
        rooms[8]->addExit("south", "Hallway");
        rooms[8]->addExit("up", "Attic");
        rooms[9]->addExit("down", "Study");
        rooms[5]->addExit("south", "Balcony");
        rooms[10]->addExit("north", "Hallway");
        rooms[7]->addExit("down", "Bedroom");
        rooms[2]->addExit("up", "Toilet");

        if (is_fixed())//固定地图下,Monster在Toilet,Princess在Attic
        {
            rooms[7]->set_monster();
            rooms[9]->set_princess();
        }
        else
        {
            mt19937 gen(static_cast<unsigned>(time(nullptr))); // 用时间生成随机数
            uniform_int_distribution<> dist(1, rooms.size() - 1);

            // 随机生成monster和princess
            int PrincessIndex, MonsterIndex;
            do
            {
                PrincessIndex = dist(gen);
                MonsterIndex = dist(gen);
            } while (PrincessIndex == MonsterIndex);
            rooms[MonsterIndex]->set_monster();
            rooms[PrincessIndex]->set_princess();
        }
    }

    bool check_game(int current_room, int flag) const // 判断游戏结束的的条件
    {
        if (current_room == 0 && flag == 1)
        {
            cout << "Win!!" << endl;
            return true;
        }
        if (rooms[current_room]->has_monster())
        {
            cout << "You have been eaten by the monster! Game over!" << endl;
            return true;
        }
        return false;
    }

    void play() // 开始游戏
    {
        int current_room = 0; // 开始位于lobby
        string input;
        int flag = 0; // 是否救出公主
        while (true)
        {
            rooms[current_room]->description(); // 描述房间
            if (check_game(current_room, flag)) // 判断游戏是否结束
                break;
            if (rooms[current_room]->has_princess()) // 遇见公主
            {
                cout << "You have saved the princess!" << endl;
                cout << "Princess: Thank u for saving me! Please take me to the lobby!" << endl;
                flag = 1;
            }
            cout << "Enter your command(e.g. 'go north'):";
            string go;
            string direction;
            cin >> go >> direction;

            if (go != "go") // 输入异常检测
            {
                cout << "Invalid command. Use 'go <direction>'." << endl;
                cin.clear();                                         // 清除错误状态
                cin.ignore(numeric_limits<streamsize>::max(), '\n'); // 忽略剩余输入
                continue;
            }
            string next_room = rooms[current_room]->get_exit(direction);

            if (next_room.empty()) // 不存在的房间
            {
                cout << "There is no exit in that direction!" << endl;
                continue;
            }

            if (room_index.find(next_room) != room_index.end()) // 找下一个房间的索引
            {
                current_room = room_index[next_room];
            }
        }
    }
};