#include "Diary.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: pdremove <date>" << std::endl;
        return -1;
    }
    std::string date = argv[1];
    DiaryBook diarybook("diary.txt");
    diarybook.load();

    if (!is_valid_date(date))
    {
        std::cout << "Date format: YYYY-MM-DD." << std::endl;
        return -1;
    }

    if (!diarybook.remove_diary(date))
        return -1;
    if (!diarybook.save()) // 保存更改
        return -1;

    return 0;
}