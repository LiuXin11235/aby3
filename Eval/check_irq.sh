#!/bin/bash

irq=$1

if [ -z "$irq" ]; then
    echo "请输入一个IRQ号码"
    exit 1
fi

if [ ! -e "/proc/irq/$irq/smp_affinity_list" ]; then
    echo "$irq DO NOT EXIST!"
    exit 1
fi

echo "(HARD) IRQ $irq is processing on CPU:"

declare -A cpu_interrupts1
declare -A cpu_interrupts2

awk -v irq=$irq '$1 ~ irq+":" { for(i=2; i<=NF; i++) print "CPU" i-2 ": " $i }' /proc/interrupts > /tmp/interrupts1

sleep 1  

awk -v irq=$irq '$1 ~ irq+":" { for(i=2; i<=NF; i++) print "CPU" i-2 ": " $i }' /proc/interrupts > /tmp/interrupts2

diff /tmp/interrupts1 /tmp/interrupts2
rm /tmp/interrupts1 /tmp/interrupts2
# for cpu in "${!cpu_interrupts1[@]}"; do
#     if [ "${cpu_interrupts1[$cpu]}" != "${cpu_interrupts2[$cpu]}" ]; then
#         echo "CPU $cpu: ${cpu_interrupts1[$cpu]} -> ${cpu_interrupts2[$cpu]}"
#     fi
# done