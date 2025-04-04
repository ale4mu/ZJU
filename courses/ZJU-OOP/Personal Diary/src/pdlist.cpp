#include "Diary.h"

int main(int argc, char *argv[])
{
    DiaryBook diarybook("diary.txt");
    diarybook.load();

    if (argc == 3)
    {
        std::string start = argv[1];
        std::string end = argv[2];
        if (!is_valid_date(start) || !is_valid_date(end)) // 处理异常日期
        {
            std::cout << "Invalid date format. Please use YYYY-MM-DD." << std::endl;
            exit(-1);
        }
        if (start > end)
        {
            std ::cout << "Start date should be less than end date!" << std::endl;
            exit(-1);
        }
        diarybook.list_diary(start, end);
    }
    else if (argc != 1)
    {
        std::cout << "Usage: pdlist <start_date> <end_date>" << std::endl;
        exit(-1);
    }
    else
    {
        diarybook.list_diary();
    }
}