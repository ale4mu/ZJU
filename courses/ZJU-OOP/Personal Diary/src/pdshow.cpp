#include "Diary.h"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage: pdshow <date>" << std::endl;
        exit(-1);
    }
    std::string date = argv[1];
    DiaryBook diarybook("diary.txt");
    diarybook.load();
    if (!is_valid_date(date))
    {
        std::cout << "Date format: YYYY-MM-DD." << std::endl;
        exit(-1);
    }
    diarybook.get_diary(date);
}