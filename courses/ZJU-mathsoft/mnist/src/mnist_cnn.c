#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
// 参数维度常量
#define w_conv1_1 8
#define w_conv1_2 1
#define w_conv1_3 3
#define w_conv1_4 3
#define b_conv1_len 8
#define w_conv2_1 16
#define w_conv2_2 8
#define w_conv2_3 3
#define w_conv2_4 3
#define b_conv2_len 16
#define w_fc1_1 32
#define w_fc1_2 784
#define b_fc1_len 32
#define w_fc2_1 10
#define w_fc2_2 32
#define b_fc2_len 10
#define IMG_W 28

uint32_t read_big_endian(FILE *fp);

uint8_t **read_mnist_images(const char *filename, uint32_t *num_images, uint32_t *rows, uint32_t *cols);

uint8_t *read_mnist_labels(const char *filename, uint32_t *num_labels);

void mnist_forward(float *x, float *W1, float *b1, float *W2, float *b2, float *W3, float *b3, float *W4, float *b4, float *y2);

void readbmp(const char *filename, float *img);

void conv2d(float *input, float *W, float *b, int in_channels, int out_channels, int img_size, float *output);

void relu(float *x, int size);

void linear(float *x, float *W, float *b, int in_dim, int out_dim, float *y);

void max_pool2d(float *input, int channels, int size, float *output);

int main(int argc, char *argv[])
{
    assert(sizeof(float) == 4); // 通不过则报错

    FILE *fp = fopen("ExportPara/parameters_cnn.bin", "rb");
    if (!fp)
    {
        printf("无法打开参数文件！\n");
        return 1;
    }

    float *b_conv1 = (float *)malloc(sizeof(float) * b_conv1_len);
    float *w_conv1 = (float *)malloc(sizeof(float) * w_conv1_1 * w_conv1_2 * w_conv1_3 * w_conv1_4);
    float *b_conv2 = (float *)malloc(sizeof(float) * b_conv2_len);
    float *w_conv2 = (float *)malloc(sizeof(float) * w_conv2_1 * w_conv2_2 * w_conv2_3 * w_conv2_4);
    float *w_fc1 = (float *)malloc(sizeof(float) * w_fc1_1 * w_fc1_2);
    float *b_fc1 = (float *)malloc(sizeof(float) * b_fc1_len);
    float *w_fc2 = (float *)malloc(sizeof(float) * w_fc2_1 * w_fc2_2);
    float *b_fc2 = (float *)malloc(sizeof(float) * b_fc2_len);

    if (!w_conv1 || !b_conv1 || !w_conv2 || !b_conv2 || !w_fc1 || !b_fc1 || !w_fc2 || !b_fc2)
    {
        printf("内存分配失败！\n");
        fclose(fp);
        return 2;
    }

    //依次读取参数
    size_t nread = fread(w_conv1, sizeof(float), w_conv1_1 * w_conv1_2 * w_conv1_3 * w_conv1_4, fp);
    if (nread != w_conv1_1 * w_conv1_2 * w_conv1_3 * w_conv1_4)
    {
        printf("w_convv1读取失败\n");
        return 3;
    }
    nread = fread(b_conv1, sizeof(float), b_conv1_len, fp);
    if (nread != b_conv1_len)
    {
        printf("b_con1读取失败\n");
        return 4;
    }
    nread = fread(w_conv2, sizeof(float), w_conv2_1 * w_conv2_2 * w_conv2_3 * w_conv2_4, fp);
    if (nread != w_conv2_1 * w_conv2_2 * w_conv2_3 * w_conv2_4)
    {
        printf("w_conv2读取失败\n");
        return 5;
    }
    nread = fread(b_conv2, sizeof(float), b_conv2_len, fp);
    if (nread != b_conv2_len)
    {
        printf("b_conv2读取失败\n");
        return 6;
    }
    nread = fread(w_fc1, sizeof(float), w_fc1_1 * w_fc1_2, fp);
    if (nread != w_fc1_1 * w_fc1_2)
    {
        printf("w_fc1读取失败\n");
        return 7;
    }
    nread = fread(b_fc1, sizeof(float), b_fc1_len, fp);
    if (nread != b_fc1_len)
    {
        printf("b_fc1读取失败\n");
        return 8;
    }
    nread = fread(w_fc2, sizeof(float), w_fc2_1 * w_fc2_2, fp);
    if (nread != w_fc2_1 * w_fc2_2)
    {
        printf("w_fc2读取失败\n");
        return 9;
    }
    nread = fread(b_fc2, sizeof(float), b_fc2_len, fp);
    if (nread != b_fc2_len)
    {
        printf("b_fc2读取失败\n");
        return 10;
    }
    fclose(fp);

    if (argc == 2)
    {
        // 将图像数据转换为浮点数
        float img[IMG_W * IMG_W];
        float y2[10]; // 输出数据

        readbmp(argv[1], img);

        mnist_forward(img, w_conv1, b_conv1, w_conv2, b_conv2, w_fc1, b_fc1, w_fc2, b_fc2, y2);
        int number = 0;
        for (int i = 0; i < 10; ++i)
        {
            if (y2[number] < y2[i])
                number = i;
        }

        printf("The number is %d\n", number);
    }
    else
    {
        const char *image_file = "data/t10k-images-idx3-ubyte";
        const char *label_file = "data/t10k-labels-idx1-ubyte";

        uint32_t num_images, rows, cols, num_labels;
        uint8_t **images = read_mnist_images(image_file, &num_images, &rows, &cols);
        uint8_t *labels = read_mnist_labels(label_file, &num_labels);

        int correct = 0;
        for (int i = 0; i < num_images; ++i)
        {
            float y2[10]; // 输出数据
            // 将图像数据转换为浮点数
            float img[IMG_W * IMG_W];
            for (int j = 0; j < IMG_W * IMG_W; ++j)
            {
                img[j] = (images[i][j] / 255.0f - 0.1307f) / 0.3081f; // 标准化，与Python中一致
            }

            mnist_forward(img, w_conv1, b_conv1, w_conv2, b_conv2, w_fc1, b_fc1, w_fc2, b_fc2, y2);
            int number = 0;
            for (int i = 0; i < 10; ++i)
            {
                if (y2[number] < y2[i])
                    number = i;
            }
            if (number == labels[i])
                correct++;
        }

        printf("The correction is: %f\n", correct * 1.0 / num_images);
    }
    // 释放内存
    free(w_conv1);
    free(b_conv1);
    free(w_conv2);
    free(b_conv2);
    free(w_fc1);
    free(b_fc1);
    free(w_fc2);
    free(b_fc2);

    return 0;
}

void readbmp(const char *filename, float *img)
{
    FILE *f = fopen(filename, "rb");
    if (!f)
    {
        printf("Can not open file %s\n", filename);
        return;
    }

    // 读取文件头（14字节）
    uint8_t file_header[14];
    fread(file_header, 1, 14, f);

    // 读取信息头（40字节）
    uint8_t info_header[40];
    fread(info_header, 1, 40, f);

    // 对齐
    int width = info_header[4] | (info_header[5] << 8) | (info_header[6] << 16) | (info_header[7] << 24);
    int height = info_header[8] | (info_header[9] << 8) | (info_header[10] << 16) | (info_header[11] << 24);
    uint16_t bit_count = *(uint16_t *)&info_header[14];
    if (bit_count != 24)
    {
        printf("Only 24-bit BMP is supported.\n");
        fclose(f);
    }

    // 获取像素数据偏移
    uint32_t pixel_offset = *(uint32_t *)&file_header[10];
    fseek(f, pixel_offset, SEEK_SET);

    // 每行字节数4字节对齐
    int row_bytes = ((width * 3 + 3) / 4) * 4;
    uint8_t *row_buf = malloc(row_bytes);

    // BMP像素自下而上存储
    for (int y = height - 1; y >= 0; --y)
    {
        fread(row_buf, 1, row_bytes, f);
        for (int x = 0; x < width; ++x)
        {
            uint8_t B = row_buf[x * 3 + 0];
            uint8_t G = row_buf[x * 3 + 1];
            uint8_t R = row_buf[x * 3 + 2];
            float gray = (0.299f * R + 0.587f * G + 0.114f * B) / 255.0f;
            gray = (gray - 0.1307f) / 0.3081f;
            img[y * width + x] = gray;
        }
    }

    free(row_buf);
    fclose(f);
}

// 读取大端整数（4字节）
uint32_t read_big_endian(FILE *fp)
{
    uint8_t bytes[4];
    fread(bytes, sizeof(uint8_t), 4, fp);
    return (bytes[0] << 24) | (bytes[1] << 16) | (bytes[2] << 8) | bytes[3];
}

// 读取图像数据
uint8_t **read_mnist_images(const char *filename, uint32_t *num_images, uint32_t *rows, uint32_t *cols)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
    {
        perror("Failed to open image file");
        exit(1);
    }

    uint32_t magic = read_big_endian(fp);
    if (magic != 0x00000803)
    {
        fprintf(stderr, "Invalid image file magic number: %08X\n", magic);
        exit(1);
    }

    *num_images = read_big_endian(fp);
    *rows = read_big_endian(fp);
    *cols = read_big_endian(fp);

    uint8_t **images = malloc(*num_images * sizeof(uint8_t *));
    for (uint32_t i = 0; i < *num_images; ++i)
    {
        images[i] = malloc((*rows) * (*cols) * sizeof(uint8_t));
        fread(images[i], sizeof(uint8_t), (*rows) * (*cols), fp);
    }

    fclose(fp);
    return images;
}

// 读取标签数据
uint8_t *read_mnist_labels(const char *filename, uint32_t *num_labels)
{
    FILE *fp = fopen(filename, "rb");
    if (!fp)
    {
        perror("Failed to open label file");
        exit(1);
    }

    uint32_t magic = read_big_endian(fp);
    if (magic != 0x00000801)
    {
        fprintf(stderr, "Invalid label file magic number: %08X\n", magic);
        exit(1);
    }

    *num_labels = read_big_endian(fp);
    uint8_t *labels = malloc(*num_labels * sizeof(uint8_t));
    fread(labels, sizeof(uint8_t), *num_labels, fp);

    fclose(fp);
    return labels;
}

// relu，把小于0的置为0，非线性
void relu(float *x, int size)
{
    for (int i = 0; i < size; ++i)
        if (x[i] < 0)
            x[i] = 0;
}

// fc1和fc2
void linear(float *x, float *W, float *b, int in_dim, int out_dim, float *y)
{
    for (int i = 0; i < out_dim; ++i)
    {
        float acc = 0.0f;
        for (int j = 0; j < in_dim; ++j)
        {
            float val = W[j*out_dim + i];
            if(val != 0.0f)
                acc += val * x[j];
        }
        y[i] = acc + b[i];
    }
}

// 最大池化函数，每次池化尺寸减半
void max_pool2d(float *input, int channels, int size, float *output)
{
    for (int c = 0; c < channels; ++c)
    {
        for (int i = 0; i < size / 2; ++i)
        {
            for (int j = 0; j < size / 2; ++j)
            {
                float max_val = -INFINITY;
                for (int di = 0; di < 2; ++di)
                {
                    for (int dj = 0; dj < 2; ++dj)
                    {
                        float val = input[c * size * size + (2 * i + di) * size + (2 * j + dj)];
                        if (val > max_val)
                            max_val = val;
                    }
                }
                output[c * (size / 2) * (size / 2) + i * (size / 2) + j] = max_val;
            }
        }
    }
}

// 卷积函数，padding=1
void conv2d(float *input, float *W, float *b, int in_channels, int out_channels, int img_size, float *output)
{
    int pad = 1;
    int padded_size = img_size + 2 * pad;

    // 动态分配三维数组：in_channels x padded_size x padded_size
    float ***padded = (float ***)malloc(in_channels * sizeof(float **));
    for (int c = 0; c < in_channels; ++c)
    {
        padded[c] = (float **)malloc(padded_size * sizeof(float *));
        for (int i = 0; i < padded_size; ++i)
        {
            padded[c][i] = (float *)calloc(padded_size, sizeof(float)); // 初始化为0
        }
    }

    // 填充每个通道
    for (int c = 0; c < in_channels; ++c)
    {
        for (int i = 0; i < img_size; ++i)
        {
            for (int j = 0; j < img_size; ++j)
            {
                padded[c][i + pad][j + pad] = input[c * img_size * img_size + i * img_size + j];
            }
        }
    }

    // 执行卷积
    for (int oc = 0; oc < out_channels; ++oc)
    {
        for (int i = 0; i < img_size; ++i)
        {
            for (int j = 0; j < img_size; ++j)
            {
                float acc = 0.0f;
                for (int ic = 0; ic < in_channels; ++ic)
                {
                    for (int ki = 0; ki < 3; ++ki)
                    {
                        for (int kj = 0; kj < 3; ++kj)
                        {
                            int weight_idx = oc * (in_channels * 3 * 3) + ic * (3 * 3) + ki * 3 + kj;
                            float weight_val = W[weight_idx];
                            if(weight_val != 0.0f)
                            {
                                float val = padded[ic][i + ki][j + kj];
                                acc += W[weight_idx] * val;
                            }
                        }
                    }
                }
                output[oc * img_size * img_size + i * img_size + j] = acc + b[oc];
            }
        }
    }

    // 释放内存
    for (int c = 0; c < in_channels; ++c)
    {
        for (int i = 0; i < padded_size; ++i)
        {
            free(padded[c][i]);
        }
        free(padded[c]);
    }
    free(padded);
}
void mnist_forward(float *x, float *W1, float *b1, float *W2, float *b2, float *W3, float *b3, float *W4, float *b4, float *y2)
{
    // 中间变量
    float conv1_out[w_conv1_1][IMG_W][IMG_W];
    float conv2_out[w_conv2_1][IMG_W][IMG_W];
    float pool1_out[w_conv2_1][14][14];
    float pool2_out[w_conv2_1][7][7];
    float fc1_in[w_fc1_2];
    float fc1_out[w_fc1_1];

    // Step 1: Conv1 + ReLU
    conv2d(x, W1, b1, 1, w_conv1_1, IMG_W, (float *)conv1_out);
    relu((float *)conv1_out, w_conv1_1 * IMG_W * IMG_W);

    // Step 2: Conv2 + ReLU
    conv2d((float *)conv1_out, W2, b2, w_conv1_1, w_conv2_1, IMG_W, (float *)conv2_out);
    relu((float *)conv2_out, w_conv2_1 * IMG_W * IMG_W);

    // Step 3: MaxPool -> (B, 16, 14, 14)
    max_pool2d((float *)conv2_out, w_conv2_1, IMG_W, (float *)pool1_out);

    // Step 4: MaxPool -> (B, 16, 7, 7)
    max_pool2d((float *)pool1_out, w_conv2_1, 14, (float *)pool2_out);

    // Step 5: 展平
    int idx = 0;
    for (int c = 0; c < w_conv2_1; ++c)
    {
        for (int i = 0; i < 7; ++i)
        {
            for (int j = 0; j < 7; ++j)
            {
                fc1_in[idx++] = pool2_out[c][i][j];
            }
        }
    }

    // Step 6: FC1 + ReLU
    linear(fc1_in, W3, b3, w_fc1_2, w_fc1_1, fc1_out);
    relu(fc1_out, w_fc1_1);

    // Step 7: FC2
    linear(fc1_out, W4, b4, w_fc1_1, w_fc2_1, y2);
}