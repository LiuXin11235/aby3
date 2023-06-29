
# ./Eval/basic_task_evaluate.sh;
# wait;

# target_time="22"

# while true; do
#   current_time=$(date -u +%H)
#   if [[ $current_time == $target_time ]]; then
#       # 当前时间达到目标时间，启动你的程序
#       ./Eval/hd_task_evaluate.sh
#       wait;

#       ./Eval/multi_steps_evaluate.sh
#       wait;
#       break  # 退出循环
#   else
#     echo "check time: "$current_time
#   fi
#   # 休眠一段时间后再次检查
#   sleep 1800;  # 每分钟检查一次
# done;

# ./Eval/basic_task_evaluate.sh
# wait;

./Eval/multi_steps_evaluate.sh
wait;

./Eval/hd_task_evaluate.sh
wait;