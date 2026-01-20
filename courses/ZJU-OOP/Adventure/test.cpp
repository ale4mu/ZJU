#include "adventure.h"
#include <sstream>
#include <fstream>

int main()
{
    Game game;
    game.set_fixed_map(); // 固定地图，Monster在Toilet,Princess在Attic
    game.set_game();

    // 用isstringstream的方法模拟用户输入
    string input =
        "go up\n"    // Lobby -> Hallway
        "go north\n" // Hallway -> Study
        "go up\n"    // Study -> Attic (成功找到公主)
        "go down\n"  // Attic -> Study
        "go south\n" // Study -> Hallway
        "go down\n"; // Hallway -> Lobby (win)

    istringstream inputstream(input);
    streambuf *origCin = cin.rdbuf(inputstream.rdbuf()); // 将cin的缓冲区替换为inputsream的缓冲区

    game.play();
    cin.rdbuf(origCin); // 恢复原来的cin缓冲区
}