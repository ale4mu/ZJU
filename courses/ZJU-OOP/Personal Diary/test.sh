#!/bin/bash

# 初始化
clean() {
    echo "清理旧文件"
    rm -f diary.txt
    make clean > /dev/null
}

# 编译程序
build() {
    echo "编译程序"
    make || { echo "编译失败！"; exit 1; }
}

# 测试1：添加日记
test_pdadd() {
    echo -e "\n测试1：添加日记(pdadd)"
    echo "美好的一天" | ./pdadd 2022-06-02
    echo "不错的一天!" | ./pdadd 2023-03-12
    echo "青春结束了!" | ./pdadd 2023-06-08
    echo "心碎了!" | ./pdadd 2023-08-22
    echo "完成了OOP作业!" | ./pdadd 2025-04-04
    ./pdlist
    echo "内容已添加"
}

# 测试2：添加多行日记
test_pdadd_more() {
    echo -e "\n测试2：添加多行日记"
    echo -e "哈哈哈\n嘻嘻嘻\n啦啦啦" | ./pdadd 2023-06-27
    echo "多行内容已添加"
}

# 测试3：覆盖相同日期的日记
test_overwrite() {
    echo -e "\n测试3：覆盖日记(./pdadd 2023-03-12)"
    echo "值得回忆。" | ./pdadd 2023-03-12
    ./pdlist
    echo "内容已覆盖"
}

# 测试4：列出所有日记
test_list_all() {
    echo -e "\n测试4：按时间顺序列出所有日记(./pdlist)"
    ./pdlist
}

# 测试5：按日期范围列出日记
test_list_range() {
    echo -e "\n测试5：列出指定范围日记 (2023-03-01 到 2023-09-01)"
    ./pdlist 2023-03-01 2023-09-01
}

# 测试6：查看特定日记
test_show() {
    echo -e "\n测试6：查看日记 (2025-04-04)"
    ./pdshow 2025-04-04
}

# 测试7：删除日记
test_remove() {
    echo -e "\n测试7：删除日记 (2023-08-22)"
    ./pdremove 2023-08-22
    echo "删除后："
    ./pdlist
}

# 主流程
main() {
    clean
    build
    
    test_pdadd
    test_pdadd_more
    test_overwrite
    
    test_list_all
    test_list_range
    test_show
    test_remove

    echo -e "\n所有测试完成！最终日记文件内容："
    cat diary.txt
}

main