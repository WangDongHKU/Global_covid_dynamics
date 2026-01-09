#!/bin/bash

# 简单进度查看脚本
PROGRESS_FILE="/home/dongw21/Global_COVID/standard_week/progress_log.txt"

echo "=== Stan模型进度查看 ==="
echo "进度文件: $PROGRESS_FILE"
echo ""

if [ ! -f "$PROGRESS_FILE" ]; then
    echo "❌ 进度文件不存在，模型可能还未开始运行"
    exit 1
fi

echo "📊 最新进度:"
echo "----------------------------------------"
tail -20 "$PROGRESS_FILE"
echo "----------------------------------------"

# 检查是否完成
if grep -q "=== Stan模型运行完成 ===" "$PROGRESS_FILE"; then
    echo ""
    echo "✅ 模型运行已完成！"
    
    # 显示总运行时间
    if grep -q "总运行时间:" "$PROGRESS_FILE"; then
        echo "⏱️  运行时间统计:"
        grep "总运行时间:" "$PROGRESS_FILE" | tail -1
    fi
else
    echo ""
    echo "🔄 模型正在运行中..."
    
    # 显示当前阶段
    echo "📋 当前阶段:"
    tail -5 "$PROGRESS_FILE" | grep -E "(开始|完成|采样|处理)" | tail -1
fi
