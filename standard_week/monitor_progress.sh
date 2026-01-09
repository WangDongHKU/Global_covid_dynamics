#!/bin/bash

# 进度监控脚本
# 使用方法: ./monitor_progress.sh [job_id]

PROGRESS_FILE="/home/dongw21/Global_COVID/standard_week/progress_log.txt"
LOG_FILE="/home/dongw21/Global_COVID/standard_week/progress_log.txt"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# 检查文件是否存在
if [ ! -f "$PROGRESS_FILE" ]; then
    echo -e "${RED}错误: 进度文件不存在: $PROGRESS_FILE${NC}"
    echo "请确保Stan模型正在运行..."
    exit 1
fi

echo -e "${CYAN}=== Stan模型进度监控器 ===${NC}"
echo -e "${YELLOW}进度文件: $PROGRESS_FILE${NC}"
echo -e "${YELLOW}按 Ctrl+C 退出监控${NC}"
echo ""

# 显示最后几行并持续监控
tail -f "$PROGRESS_FILE" | while read line; do
    # 根据内容类型添加颜色
    if [[ $line == *"错误"* ]] || [[ $line == *"Error"* ]]; then
        echo -e "${RED}$line${NC}"
    elif [[ $line == *"完成"* ]] || [[ $line == *"完成"* ]]; then
        echo -e "${GREEN}$line${NC}"
    elif [[ $line == *"开始"* ]] || [[ $line == *"开始"* ]]; then
        echo -e "${BLUE}$line${NC}"
    elif [[ $line == *"采样"* ]] || [[ $line == *"sampling"* ]]; then
        echo -e "${PURPLE}$line${NC}"
    elif [[ $line == *"时间"* ]] || [[ $line == *"time"* ]]; then
        echo -e "${YELLOW}$line${NC}"
    else
        echo "$line"
    fi
done
