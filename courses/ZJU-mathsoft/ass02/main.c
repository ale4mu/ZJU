#include <ncurses.h>         ///< 用于创建文本用户界面
#include <unistd.h>          ///< 用于 usleep 函数
#include <stdlib.h>          ///< 用于 atoi 和 malloc 等函数
#include <string.h>          ///< 用于字符串操作
#include "Conway.h"          ///< 包含游戏逻辑和数据结构的定义
#include <time.h>            ///< 用于获取当前时间
#include "stb_image_write.h" ///< 用于保存图像

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        printf("Usage: %s width height rle_file[|rand]\n", argv[0]);
        return 1;
    }

    boundary = 1; ///< 设置边界模式为循环边界。

    initscr();             ///< 初始化 curses 库。
    noecho();              ///< 关闭输入回显。
    curs_set(FALSE);       ///< 隐藏光标。
    cbreak();              ///< 禁止行缓冲，使得 getch 立即返回。
    nodelay(stdscr, TRUE); ///< 设置非阻塞模式。

    init_cells(atoi(argv[1]), atoi(argv[2])); ///< 初始化细胞状态。

    if (argc < 4)
        set_rand_cells(); ///< 生成随机模式。
    else
        load_pattern(10, 10, argv[3]); ///< 加载指定的模式文件。

    int width = world.width;   ///< 游戏区域的宽度。
    int height = world.height; ///< 游戏区域的高度。
    int **cell = world.cell;   ///< 游戏区域的细胞状态。

    while (1)
    {
        clear();        ///< 清除屏幕。
        draw_cells();   ///< 绘制细胞。
        refresh();      ///< 刷新屏幕，将绘制的细胞显示出来。
        usleep(10000);  ///< 睡眠一段时间，控制游戏速度。
        update_cells(); ///< 更新细胞状态。
        usleep(10000);  ///< 睡眠一段时间，控制游戏速度。

        int ch = getch(); ///< 检测按键输入。
        if (ch == 'q')    ///< 按下 'q' 退出程序。
        {
            break;
        }
        else if (ch == 's') ///< 按下 's' 保存当前游戏状态为图像。
        {
            save_image(width, height, cell); ///< 保存图像。
        }
        update_cells(); ///< 更新细胞状态。
        usleep(10000);  ///< 睡眠一段时间，控制游戏速度。
    }

    endwin();     ///< 结束 curses 模式。
    free_cells(); ///< 释放内存。

    return 0; ///< 程序正常退出。
}