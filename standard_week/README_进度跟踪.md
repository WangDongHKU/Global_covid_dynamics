# Stan模型进度跟踪系统

## 概述
这个进度跟踪系统允许你在SLURM环境中监控Stan模型的运行进度，即使无法直接看到控制台输出。

## 文件说明

### 主要文件
- `standard_model_10.R` - 已添加进度跟踪功能的Stan模型脚本
- `progress_log.txt` - 进度日志文件（运行时会自动创建）
- `monitor_progress.sh` - 实时进度监控脚本
- `check_progress.sh` - 简单进度查看脚本
- `submit_stan_job.sh` - SLURM提交脚本示例

## 使用方法

### 1. 直接运行R脚本
```bash
cd /home/dongw21/Global_COVID/standard_week
Rscript standard_model_10.R
```

### 2. 使用SLURM提交任务
```bash
# 修改submit_stan_job.sh中的邮箱地址
sbatch submit_stan_job.sh
```

### 3. 监控进度

#### 实时监控（推荐）
```bash
./monitor_progress.sh
```
这会显示彩色的实时进度更新，包括：
- 🔵 蓝色：开始阶段
- 🟢 绿色：完成阶段  
- 🟣 紫色：采样阶段
- 🟡 黄色：时间信息
- 🔴 红色：错误信息

#### 简单查看
```bash
./check_progress.sh
```
这会显示最新的20行进度信息和当前状态。

#### 手动查看
```bash
# 查看最新进度
tail -f /home/dongw21/Global_COVID/standard_week/progress_log.txt

# 查看全部进度
cat /home/dongw21/Global_COVID/standard_week/progress_log.txt

# 搜索特定信息
grep "采样" /home/dongw21/Global_COVID/standard_week/progress_log.txt
```

## 进度信息说明

### 主要阶段
1. **初始化** - 加载库和设置环境
2. **数据加载** - 加载预训练模型、病例数据、移动性数据等
3. **模型编译** - 编译Stan模型
4. **数据准备** - 准备模型输入数据
5. **MCMC初始化** - 初始化MCMC链
6. **Stan采样** - 执行MCMC采样（最耗时的阶段）
7. **结果处理** - 提取和保存结果
8. **可视化** - 生成图表

### 采样进度
在采样阶段，系统会记录：
- 每个链的进度
- 当前阶段（warmup/sampling）
- 迭代次数和百分比
- 预计剩余时间

## 在SLURM环境中的使用

### 1. 提交任务
```bash
sbatch submit_stan_job.sh
```

### 2. 查看任务状态
```bash
squeue -u $USER
```

### 3. 监控进度
```bash
# 在另一个终端中运行
./monitor_progress.sh
```

### 4. 查看输出
```bash
# 查看SLURM输出
cat stan_job_<job_id>.out

# 查看错误日志
cat stan_job_<job_id>.err
```

## 自定义配置

### 修改进度文件路径
在`standard_model_10.R`中修改：
```r
progress_file <- "/your/custom/path/progress_log.txt"
```

### 修改进度记录频率
在采样回调函数中修改：
```r
if (iter %% 50 == 0) {  # 每50次迭代记录一次
```

### 添加更多进度点
在代码中添加：
```r
log_progress("你的自定义消息")
```

## 故障排除

### 进度文件不存在
- 确保模型正在运行
- 检查文件路径权限
- 确认工作目录正确

### 进度更新缓慢
- 检查磁盘空间
- 确认文件系统性能
- 考虑减少进度记录频率

### SLURM任务失败
- 检查资源请求是否合理
- 查看错误日志文件
- 确认依赖库已安装

## 性能优化建议

1. **减少进度记录频率** - 对于长时间运行的任务，可以减少记录频率
2. **使用SSD存储** - 将进度文件放在SSD上以提高写入速度
3. **定期清理** - 定期清理旧的进度文件
4. **监控磁盘空间** - 确保有足够的磁盘空间

## 示例输出

```
[2024-01-15 10:30:00] === Stan模型运行开始 ===
[2024-01-15 10:30:01] 初始化进度跟踪系统
[2024-01-15 10:30:02] 加载预训练模型数据...
[2024-01-15 10:30:05] 预训练模型数据加载完成
[2024-01-15 10:30:06] 编译Stan模型...
[2024-01-15 10:35:20] Stan模型编译完成
[2024-01-15 10:35:21] 开始Stan采样...
[2024-01-15 10:35:22] 链 1: warmup阶段 - 迭代 50 / 600 (8.3%)
[2024-01-15 10:35:23] 链 2: warmup阶段 - 迭代 50 / 600 (8.3%)
...
[2024-01-15 15:20:30] Stan采样完成！
[2024-01-15 15:20:31] 总采样时间: 4.75 小时
[2024-01-15 15:20:32] === Stan模型运行完成 ===
```
