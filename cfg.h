#ifndef CFG_H_
#define CFG_H_
#include <bits/stdc++.h>

using namespace std;

#include <boost/algorithm/string.hpp>

#include "config.h"
#include "log.h"
#include "gz_reader.h"
#include "convert.h"
#include "settings.h"

class BBL
{
public:
    uint64_t total_count;
    unordered_map<uint64_t,uint64_t> neighbor_count;
    uint32_t instrs;
    uint32_t bytes;
    uint64_t average_cycles;
    BBL(uint32_t instr_count, uint32_t byte_size)
    {
        total_count=0;
        instrs=instr_count;
        bytes=byte_size;
        average_cycles=instr_count;
    }
    void increment()
    {
        total_count+=1;
    }
    void add_neighbor(uint64_t neighbor_bbl)
    {
        if(neighbor_count.find(neighbor_bbl)==neighbor_count.end())
        {
            neighbor_count[neighbor_bbl]=1;
        }
        else
        {
            neighbor_count[neighbor_bbl]+=1;
        }
    }
    double get_neighbor_ratio(uint64_t neighbor_bbl)
    {
        if(total_count<1)panic("Neighbor ratio is called on a BBL with 0 dynamic count");
        double denominator = total_count;
        double numerator = 0;
        if(neighbor_count.find(neighbor_bbl)!=neighbor_count.end())numerator+=neighbor_count[neighbor_bbl];
        return numerator/denominator;
    }
    uint64_t get_count()
    {
        return total_count;
    }
};

class CFG
{
private:
    unordered_map<uint64_t,BBL *> bbl_infos;
    uint64_t static_size;
    uint64_t static_count;
    Settings *sim_settings;
    set<uint64_t> self_modifying_bbls;
public:
    CFG(MyConfigWrapper &conf, Settings &settings)
    {
        sim_settings = &settings;

        string self_modifying_bbl_info_path = conf.lookup("bbl_info_self_modifying", "/tmp/");
        if(self_modifying_bbl_info_path == "/tmp/")
        {
            //do nothing
            //since self modifying bbls are optional only for jitted codes
        }
        else
        {
            vector<string> raw_self_modifying_bbl_info_data;
            read_full_file(self_modifying_bbl_info_path, raw_self_modifying_bbl_info_data);
            for(uint64_t i = 0; i<raw_self_modifying_bbl_info_data.size(); i++)
            {
                string line = raw_self_modifying_bbl_info_data[i];
                boost::trim_if(line,boost::is_any_of(","));
                vector<string> parsed;
                boost::split(parsed,line,boost::is_any_of(",\n"),boost::token_compress_on);
                if(parsed.size()!=3)
                {
                    panic("A line on the self-modifying BBL info file does not have exactly three tuples");
                }
                uint64_t bbl_address = string_to_u64(parsed[0]);
                self_modifying_bbls.insert(bbl_address);
            }
            raw_self_modifying_bbl_info_data.clear();
        }
        
        string bbl_info_path = conf.lookup("bbl_info", "/tmp/");
        if(bbl_info_path == "/tmp/")
        {
            panic("BBL Info file is not present in the config");
        }
        vector<string> raw_bbl_info_data;
        read_full_file(bbl_info_path, raw_bbl_info_data);
        static_count = 0;
        static_size = 0;
        for(uint64_t i = 0; i<raw_bbl_info_data.size(); i++)
        {
            string line = raw_bbl_info_data[i];
            boost::trim_if(line,boost::is_any_of(","));
            vector<string> parsed;
            boost::split(parsed,line,boost::is_any_of(",\n"),boost::token_compress_on);
            if(parsed.size()!=3)
            {
                panic("A line on the BBL info file does not have exactly three tuples");
            }
            uint64_t bbl_address = string_to_u64(parsed[0]);
            uint32_t instrs = string_to_u32(parsed[1]);
            static_count+=instrs;
            uint32_t bytes = string_to_u32(parsed[2]);
            static_size+=bytes;
            if(bbl_infos.find(bbl_address)==bbl_infos.end())
            {
                bbl_infos[bbl_address]=new BBL(instrs, bytes);
            }
            else
            {
                panic("BBL info file includes multiple info entries for the same BBL start address");
            }
            parsed.clear();
        }
        raw_bbl_info_data.clear();
        
        vector<pair<uint64_t,uint64_t>> ordered_bbl_trace;
        string bbl_log_path = conf.lookup("bbl_log", "/tmp/");
        if(bbl_log_path == "/tmp/")
        {
            panic("BBL Log file is not present in the config");
        }
        vector<string> raw_bbl_log_data;
        read_full_file(bbl_log_path, raw_bbl_log_data);
        uint64_t prev[] = {0, 0};
        unordered_map<uint64_t,pair<uint64_t,uint64_t>> sum_counts;
        for(uint64_t i = 0; i< raw_bbl_log_data.size();i++)
        {
            string line = raw_bbl_log_data[i];
            boost::trim_if(line,boost::is_any_of(","));
            vector<string> parsed;
            boost::split(parsed,line,boost::is_any_of(",\n"),boost::token_compress_on);
            if(parsed.size()!=2)
            {
                panic("A line on the BBL log file does not have exactly two tuples");
            }
            uint64_t bbl_address = string_to_u64(parsed[0]);
            uint64_t cycle = string_to_u64(parsed[1]);
            if(prev[0]!=0)
            {
                if(sum_counts.find(prev[0])!=sum_counts.end())
                {
                    sum_counts[prev[0]].first+=cycle;
                    sum_counts[prev[0]].second++;
                }
                else
                {
                    sum_counts[prev[0]]=make_pair(cycle,1);
                }
            }
            prev[0]=prev[1];
            prev[1]=bbl_address;
            ordered_bbl_trace.push_back(make_pair(bbl_address,cycle));
            parsed.clear();
        }
        raw_bbl_log_data.clear();
        
        for(auto it: sum_counts)
        {
            uint64_t bbl_address = it.first;
            uint64_t sum = it.second.first;
            uint64_t count = it.second.second;
            if(bbl_infos.find(bbl_address)==bbl_infos.end())
            {
                panic("BBL log file includes a BBL that is not present in BBL info file");
            }
            else
            {
                bbl_infos[bbl_address]->average_cycles = 1 + ((sum - 1) / count);
            }
        }
        sum_counts.clear();
        
        uint64_t min_distance = sim_settings->get_min_distance();
        uint64_t max_distance = sim_settings->get_max_distance();
        
        for(uint64_t i = 0; i < ordered_bbl_trace.size(); i++)
        {
            uint64_t bbl_address = ordered_bbl_trace[i].first;
            if(bbl_infos.find(bbl_address)==bbl_infos.end())
            {
                panic("BBL log file includes a BBL that is not present in BBL info file");
            }
            else
            {
                bbl_infos[bbl_address]->increment();
            }
            uint64_t current_distance = 0;
            set<uint64_t> already_added_neighbors;
            for(uint64_t j = i+1; j < ordered_bbl_trace.size(); j++)
            {
                if(sim_settings->multiline_mode == 1)//ASMDB
                {
                    current_distance+=bbl_infos[ordered_bbl_trace[j-1].first]->instrs;
                }
                else // if(settings.multiline_mode == 2) OUR
                {
                    current_distance+=bbl_infos[ordered_bbl_trace[j-1].first]->average_cycles;//ordered_bbl_trace[j].second;
                }
                if(current_distance<min_distance)continue;
                if(current_distance>max_distance)break;
                if(already_added_neighbors.find(ordered_bbl_trace[j].first)!=already_added_neighbors.end())continue;
                bbl_infos[bbl_address]->add_neighbor(ordered_bbl_trace[j].first);
                already_added_neighbors.insert(ordered_bbl_trace[j].first);
            }
            already_added_neighbors.clear();
        }
        ordered_bbl_trace.clear();
    }
    bool is_self_modifying(uint64_t bbl_address)
    {
        return self_modifying_bbls.find(bbl_address) != self_modifying_bbls.end();
    }
    double get_fan_out(uint64_t predecessor_bbl_address, uint64_t successor_bbl_address)
    {
        if(bbl_infos.find(predecessor_bbl_address)==bbl_infos.end())
        {
            panic("Trying to calculate fan-out of a predecessor basic block not present in the trace");
        }
        return bbl_infos[predecessor_bbl_address]->get_neighbor_ratio(successor_bbl_address);
    }
    bool is_valid_candidate(uint64_t predecessor_bbl_address, uint64_t successor_bbl_address)
    {
        if(is_self_modifying(predecessor_bbl_address))return false;
        return get_fan_out(predecessor_bbl_address, successor_bbl_address) >= sim_settings->fan_out_ratio;
    }
    uint64_t get_bbl_execution_count(uint64_t bbl_address)
    {
        if(bbl_infos.find(bbl_address)==bbl_infos.end())
        {
            panic("Trying to get execution count for a basic block not present in the trace");
        }
        return bbl_infos[bbl_address]->get_count();
    }
    uint64_t get_static_count()
    {
        return static_count;
    }
    uint64_t get_static_size()
    {
        return static_size;
    }
    uint32_t get_bbl_instr_count(uint64_t bbl_address)
    {
        if(bbl_infos.find(bbl_address)==bbl_infos.end())
        {
            panic("Trying to get instruction count for a basic block not present in the bbl-info or trace");
        }
        return bbl_infos[bbl_address]->instrs;
    }
    uint64_t get_bbl_avg_cycles(uint64_t bbl_address)
    {
        return bbl_infos[bbl_address]->average_cycles;
    }
    uint64_t get_bbl_size(uint64_t bbl_address)
    {
        return bbl_infos[bbl_address]->bytes;
    }
    void print_modified_bbl_mappings(unordered_map<uint64_t,vector<uint64_t>> &bbl_to_prefetch_targets, uint8_t size_of_prefetch_inst, ofstream *out, unordered_map<uint64_t,uint64_t> &prev_to_new)
    {
        prev_to_new.clear();
        vector<uint64_t> sorted_bbl_addrs;
        for(auto it: bbl_infos)
        {
            sorted_bbl_addrs.push_back(it.first);
        }
        sort(sorted_bbl_addrs.begin(), sorted_bbl_addrs.end());

        uint64_t start_address = 0;
        for(uint32_t i=0; i<sorted_bbl_addrs.size(); i++)
        {
            uint64_t prev_bbl_address = sorted_bbl_addrs[i];
            uint64_t new_bbl_address = prev_bbl_address;
            if(start_address<=new_bbl_address)
            {
                // no offset due to upper bbls
                start_address = new_bbl_address;
            }
            else
            {
                // offset due to upper bbls
                new_bbl_address = start_address;
            }

            prev_to_new[prev_bbl_address] = new_bbl_address;
            if(out!=nullptr)*out<<prev_bbl_address<<" "<<new_bbl_address<<endl;

            start_address += get_bbl_size(prev_bbl_address);
            if(bbl_to_prefetch_targets.find(prev_bbl_address)!=bbl_to_prefetch_targets.end())
            {
                if(bbl_to_prefetch_targets[prev_bbl_address].size()>0)
                {
                    start_address+=size_of_prefetch_inst;
                }
            }
        }
    }
};
#endif //CFG_H_
