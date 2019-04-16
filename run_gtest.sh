#/bin/bash

sort () {
    for ((i=0; i <= $((${#arr_time[@]} - 2)); ++i))
    do
        for ((j=((i + 1)); j <= ((${#arr_time[@]} - 1)); ++j))
        do
            if [[ ${arr_time[i]} -gt ${arr_time[j]} ]]
            then
                # echo $i $j ${arr[i]} ${arr[j]}
                tmp=${arr_time[i]}
                arr_time[i]=${arr_time[j]}
                arr_time[j]=$tmp

                tmp=${arr_name[i]}
                arr_name[i]=${arr_name[j]}
                arr_name[j]=$tmp
            fi
        done
    done
}

./test/gtest/gtest 2>&1 | tee gtest_file_output

targets="$(grep -ni "ms)" gtest_file_output | awk '{print $5}' | tr --delete "(" | tr '\n' ' ')"
name_targets="$(grep -ni "ms)" gtest_file_output | awk '{print $4}' | tr --delete "(" | tr '\n' ' ')"

cntr=0
IFS=$' '
for j in $targets
do
    arr_time[$cntr]=$j
    ((cntr++))
done

cntr=0
IFS=$' '
for j in $name_targets
do
    arr_name[$cntr]=$j
    ((cntr++))
done

sort ${arr_time[@]} ${arr_name[@]}

for ((i=0; i < $((${#arr_time[@]})); ++i))
do
    echo "${arr_time[i]} - ${arr_name[i]}"
done

rm -f gtest_file_output
