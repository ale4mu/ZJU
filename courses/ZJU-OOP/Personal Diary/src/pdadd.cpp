#include "Diary.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: pdadd <date>" << std::endl;
        exit(-1);
    }

    std::string date = argv[1];
    DiaryBook diarybook("diary.txt");
    diarybook.load(); // 加载日记本

    std::string content, line;
    while (true)
    {
        std::getline(std::cin, line);
        if (line == "." || std::cin.eof()) // 输入为.或EOF结束
            break;
        if (!content.empty())
            content += "\n";
        content += line;
    }

    diarybook.add_diary(date, content);
    diarybook.save();
}