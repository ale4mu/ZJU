#include "Conway.h"
#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct World world;

int boundary = 1;

void init_cells(int _w, int _h)
{
    world.cell = (int **)malloc(_w * sizeof(int *)); // 为行分配内存
    for (int i = 0; i < _w; i++)
        world.cell[i] = (int *)malloc(_h * sizeof(int)); // 为每一行分配列
    world.width = _w;
    world.height = _h;
    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < _h; j++)
        {
            world.cell[i][j] = 0; // 初始的死世界
        }
    }
}

void set_rand_cells()
{
    for (int i = 0; i < world.width; i++)
    {
        for (int j = 0; j < world.height; j++)
        {
            world.cell[i][j] = rand() % 2;
        }
    }
}

void load_pattern(int _x, int _y, char *rle_file)
{
    FILE *file = fopen(rle_file, "r");
    if (file == NULL)
    {
        exit(1);
    }
    int w, h;
    int k = 0;
    char pattern[1000];
    int width = world.width;
    int height = world.height;
    int **cell = world.cell;
    if (file != NULL)
    {
        char line[256];
        while (fgets(line, sizeof(line), file))
        {
            if (line[0] == '#')
            {
                continue;
            }
            else if (line[0] == 'x')
            {
                sscanf(line, "x = %d, y = %d", &w, &h);
            }
            else if (line[0] == 'b' ||
                     line[0] == 'o' ||
                     line[0] == '!' ||
                     line[0] == '$' ||
                     isdigit(line[0]))
            {
                for (int i = 0; i < strlen(line); i++)
                {
                    pattern[k] = line[i];
                    k++;
                }
            }
        }
    }
    else
    {
        printf("File not found\n");
    }
    pattern[k] = '\0';
    int px = _x;
    int py = _y;

    for (int i = 0; i < strlen(pattern); i++)
    {
        int count = 0;
        if (isdigit(pattern[i]))
        {
            while (isdigit(pattern[i]))
            {
                count = count * 10 + (pattern[i] - '0');
                i++;
            }
        }
        if (pattern[i] == 'b')
        {
            if (count == 0)
                count = 1;
            for (int j = 0; j < count; j++)
            {
                cell[py][px] = 0; // cell[x][y]
                py++;
            }
        }
        else if (pattern[i] == 'o')
        {
            if (count == 0)
                count = 1;
            for (int j = 0; j < count; j++)
            {
                cell[py][px] = 1; // cell[x][y]
                py++;
            }
        }
        else if (pattern[i] == '$')
        {
            px++;
            py = _y;
        }
        else if (pattern[i] == '!')
        {
            break;
        }
    }

    fclose(file);
}

void draw_cells()
{
    for (int i = 0; i < world.width; i++)
    {
        for (int j = 0; j < world.height; j++)
        {
            mvaddch(j, i, world.cell[i][j] ? '*' : ' ');
        }
    }
}

int count_neighbors(int x, int y)
{
    int count = 0;
    int width = world.width;
    int height = world.height;
    int **cell = world.cell;

    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            if (i == 0 && j == 0)
                continue; // 直接跳过自身

            if (boundary == 0) // 固定边界
            {
                int nx = x + i;
                int ny = y + j;

                if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                {
                    continue;
                }

                count += cell[nx][ny];
            }
            else if (boundary == 1) // 循环边界
            {
                int nx = (x + i + width) % width;
                int ny = (y + j + height) % height;

                count += cell[nx][ny];
            }
        }
    }

    return count;
}

void update_cells()
{
    int **new_cell;
    int **cell = world.cell;
    int width = world.width;
    int height = world.height;
    new_cell = (int **)malloc(width * sizeof(int *));
    for (int i = 0; i < width; i++)
        new_cell[i] = (int *)malloc(height * sizeof(int));

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            int count = count_neighbors(i, j);
            if (count < 2 || count > 3)
            {
                new_cell[i][j] = 0;
            }
            else if (count == 3)
            {
                new_cell[i][j] = 1;
            }
            else
            {
                new_cell[i][j] = cell[i][j];
            }
        }
    }

    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            cell[i][j] = new_cell[i][j];
        }
    }

    for (int i = 0; i < width; i++)
        free(new_cell[i]);
    free(new_cell);
}

void free_cells()
{
    for (int i = 0; i < world.width; i++)
    {
        free(world.cell[i]);
    }
    free(world.cell);
}

void save_image(int width, int height, int **cell)
{
    // 创建一个图像数组，每个像素用3个字节表示（RGB）
    unsigned char *image = (unsigned char *)malloc(width * height * 3);

    // 填充图像数据
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            int index = (y * width + x) * 3;
            if (cell[y][x] == 1) // 活细胞为白色
            {
                image[index] = 255;     // R
                image[index + 1] = 255; // G
                image[index + 2] = 255; // B
            }
            else // 死细胞为黑色
            {
                image[index] = 0;     // R
                image[index + 1] = 0; // G
                image[index + 2] = 0; // B
            }
        }
    }

    // 获取当前时间作为文件名
    time_t now = time(NULL);
    char filename[100];
    strftime(filename, sizeof(filename), "image_%Y%m%d_%H%M%S.bmp", localtime(&now));

    // 保存图像
    stbi_write_bmp(filename, width, height, 3, image);

    // 释放内存
    free(image);
}