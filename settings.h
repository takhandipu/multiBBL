#ifndef SETTINGS_H_
#define SETTINGS_H_
#include <bits/stdc++.h>

#include "config.h"
#include "convert.h"

using namespace std;

class Settings
{
public:
    uint64_t min_distance_cycle_count;
    uint64_t max_distance_cycle_count;
    uint64_t total_dyn_ins_count;
    uint64_t total_cycles;
    double dyn_ins_count_inc_ratio;
    double fan_out_ratio;
    uint64_t max_step_ahead;
    uint64_t max_prefetch_length;
    double min_ratio_percentage;
    uint64_t multiline_mode;
    double ipc;
    uint64_t insert_as_many_as_possible;
    uint64_t max_bbl_count;
    Settings(MyConfigWrapper &conf)
    {
        min_distance_cycle_count = string_to_u64(conf.lookup("min_distance_cycle_count","40"));
        max_distance_cycle_count = string_to_u64(conf.lookup("max_distance_cycle_count","200"));
        total_dyn_ins_count = string_to_u64(conf.lookup("total_dyn_ins_count","163607974"));
        total_cycles = string_to_u64(conf.lookup("total_cycles","225464000"));
        dyn_ins_count_inc_ratio = string_to_double(conf.lookup("dyn_ins_count_inc_ratio","0.01"));
        fan_out_ratio = string_to_double(conf.lookup("fan_out_ratio","0.5"));
        max_step_ahead = string_to_u64(conf.lookup("max_step_ahead", "10"));
        max_prefetch_length = string_to_u64(conf.lookup("max_prefetch_length", "8"));
        min_ratio_percentage = string_to_double(conf.lookup("min_ratio_percentage","0.99"));
        multiline_mode = string_to_u64(conf.lookup("multiline_mode", "1"));
        insert_as_many_as_possible = 1;//string_to_u64(conf.lookup("insert_as_many_as_possible", "0"));
        max_bbl_count = string_to_u64(conf.lookup("max_bbl_count", "1"));

        ipc = ((1.0*total_dyn_ins_count)/total_cycles);
    }
    uint64_t get_min_distance()
    {
        return 1;
        /*if(multiline_mode==1)//ASMDB static implementation
        {
            return 51;
        }
        else //if(multiline_mode==2) //OUR
        {
            return min_distance_cycle_count;
        }*/
    }
    uint64_t get_max_distance()
    {
        return 200;
        /*if(multiline_mode==1)//ASMDB static implementation
        {
            return 200;
        }
        else //if(multiline_mode==2) //OUR
        {
            return max_distance_cycle_count;
        }*/
    }
};
#endif //SETTINGS_H_
