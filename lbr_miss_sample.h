#ifndef LBR_MISS_SAMPLE_H_
#define LBR_MISS_SAMPLE_H_

#include <bits/stdc++.h>

using namespace std;

#include <boost/algorithm/string.hpp>

#include "gz_reader.h"
#include "log.h"
#include "convert.h"
#include "cfg.h"

int measure_prefetch_length(uint64_t current, uint64_t next)
{
    return (next>>6) - (current>>6);
}

unsigned int my_abs(int a)
{
    if(a<0)return -a;
    return a;
}

string vec_to_string(vector<uint64_t> &v)
{
    stringstream ss;
    for(size_t i = 0; i < v.size(); ++i)
    {
        if(i != 0)ss << ",";
        ss << v[i];
    }
    string s = ss.str();
    return s;
}
void string_to_vec(string &s, vector<uint64_t> &v)
{
    v.clear();
    vector<string> current_record;
    boost::split(current_record,s,boost::is_any_of(","),boost::token_compress_on);
    for(size_t i=0;i<current_record.size();i++)v.push_back(string_to_u64(current_record[i]));
}

class Candidate
{
public:
    vector<uint64_t> predecessor_bbl_address;
    double covered_miss_ratio_to_predecessor_count;
    Candidate(vector<uint64_t> &addr, double c_ratio)
    {
        for(auto it: addr)predecessor_bbl_address.push_back(it);
        covered_miss_ratio_to_predecessor_count = c_ratio;
    }
    Candidate(string &addr, double c_ratio)
    {
        string_to_vec(addr,predecessor_bbl_address);
        covered_miss_ratio_to_predecessor_count = c_ratio;
    }
    bool operator < (const Candidate& right) const
    {
        return covered_miss_ratio_to_predecessor_count < right.covered_miss_ratio_to_predecessor_count;
    }
};

struct prefetch_benefit
{
    uint64_t addr;
    uint64_t prefetch_count;
    uint64_t covered_miss_count;
    prefetch_benefit(uint64_t a, uint64_t p_count, uint64_t cm_count) : prefetch_count(p_count), covered_miss_count(cm_count)
    {
        addr=a;
    }
    bool operator < (const prefetch_benefit& right) const
    {
        if(prefetch_count < 1)return false;
        if(right.prefetch_count < 1) return true;
        double left_ratio = (1.0 * covered_miss_count) / prefetch_count;
        double right_ratio = (1.0 * right.covered_miss_count) / right.prefetch_count;
        if (left_ratio == right_ratio)
        {
            if (covered_miss_count == right.covered_miss_count)
            {
                return prefetch_count < right.prefetch_count;
            }
            else
            {
                return covered_miss_count < right.covered_miss_count;
            }
        }
        else
        {
            return left_ratio < right_ratio;
        }
    }
};

class LBR_Sample
{
public:
    vector<pair<uint64_t,uint64_t>> current_lbr;
    LBR_Sample()
    {
        current_lbr.clear();
    }
    void push_pair(uint64_t prev_bbl, uint64_t cycle)
    {
        current_lbr.push_back(make_pair(prev_bbl, cycle));
    }
    void push(string lbr_pair)
    {
        vector<string> current_record;
        boost::split(current_record,lbr_pair,boost::is_any_of(";"),boost::token_compress_on);
        if(current_record.size() != 2)panic("LB-record does not contain exactly 2 entries\n");
        push_pair(string_to_u64(current_record[0]),string_to_u64(current_record[1]));
    }
    uint64_t size()
    {
        return current_lbr.size();
    }
    void get_candidates(vector<uint64_t> &candidates, Settings *settings, CFG *dynamic_cfg)
    {
        candidates.clear();
        uint64_t current_distance = 0;
        
        
        uint64_t min_distance = settings->get_min_distance();
        uint64_t max_distance = settings->get_max_distance();
        
        for(uint64_t i = 0; i<current_lbr.size(); i++)
        {
            if(settings->multiline_mode == 1)//ASMDB
            {
                current_distance+=dynamic_cfg->get_bbl_instr_count(current_lbr[i].first);
            }
            else // if(settings.multiline_mode == 2) OUR
            {
                current_distance+=dynamic_cfg->get_bbl_avg_cycles(current_lbr[i].first);//current_lbr[i].second;
            }
            if(current_distance>=min_distance)
            {
                if(current_distance>max_distance)return;
                candidates.push_back(current_lbr[i].first);
            }
        }
    }
    bool is_covered(vector<uint64_t> &predecessor_bbl_address, Settings *settings, CFG *dynamic_cfg)
    {
        uint64_t current_distance = 0;
        
        uint64_t min_distance = settings->get_min_distance();
        uint64_t max_distance = settings->get_max_distance();

        uint64_t index = 0;
        if(predecessor_bbl_address.size()==0)return false;
        
        for(uint64_t i = 0; i<current_lbr.size(); i++)
        {
            if(settings->multiline_mode == 1)//ASMDB
            {
                current_distance+=dynamic_cfg->get_bbl_instr_count(current_lbr[i].first);
            }
            else // if(settings.multiline_mode == 2) OUR
            {
                current_distance+=dynamic_cfg->get_bbl_avg_cycles(current_lbr[i].first);//current_lbr[i].second;
            }
            if(current_distance>=min_distance)
            {
                if(current_distance>max_distance)break;
                if(current_lbr[i].first==predecessor_bbl_address[index])index++;//return true;
                if(index>=predecessor_bbl_address.size())return true;
            }
        }
        return false;
    }
    bool is_covered(string &predecessor_bbl_address, Settings *settings, CFG *dynamic_cfg)
    {
        vector<uint64_t> tmp;
        string_to_vec(predecessor_bbl_address,tmp);
        return is_covered(tmp, settings, dynamic_cfg);
    }
    uint64_t count_accessed_cache_lines(Settings *settings, CFG *dynamic_cfg)
    {
        
        uint64_t min_distance = settings->get_min_distance();
        uint64_t max_distance = settings->get_max_distance();
        
        set<uint64_t> unique_cache_lines;
        uint64_t current_distance = 0;
        
        for(uint64_t i = 0; i<current_lbr.size(); i++)
        {
            if(settings->multiline_mode == 1)//ASMDB
            {
                current_distance+=dynamic_cfg->get_bbl_instr_count(current_lbr[i].first);
            }
            else // if(settings.multiline_mode == 2) OUR
            {
                current_distance+=dynamic_cfg->get_bbl_avg_cycles(current_lbr[i].first);//current_lbr[i].second;
            }
            if(current_distance>=min_distance)
            {
                if(current_distance>max_distance)break;
                uint64_t start = current_lbr[i].first;
                uint64_t end = start + dynamic_cfg->get_bbl_size(start);
                start>>=6;
                end>>=6;
                for(uint64_t j=start; j<=end; j++)unique_cache_lines.insert(j);
            }
        }
        
        return unique_cache_lines.size();
    }
    void clear()
    {
        current_lbr.clear();
    }
};

class LBR_List
{
public:
    vector<LBR_Sample *> all_misses;
    LBR_List()
    {
        all_misses.clear();
    }
    void push_miss_sample(LBR_Sample *current_sample)
    {
        all_misses.push_back(current_sample);
    }
    void push(vector<string> lbr_pairs)
    {
        LBR_Sample *current_sample = new LBR_Sample();
        for(uint64_t i=2/*omit first two items, cache line id and missed bbl address*/;i<lbr_pairs.size();i++)
        {
            current_sample->push(lbr_pairs[i]);
        }
        push_miss_sample(current_sample);
    }
    uint64_t size()
    {
        return all_misses.size();
    }
    void measure_prefetch_window_cdf(unordered_map<uint64_t,uint64_t> &table, Settings *settings, CFG *dynamic_cfg)
    {
        for(uint64_t i=0; i<all_misses.size(); i++)
        {
            uint64_t accessed_cache_line_count = all_misses[i]->count_accessed_cache_lines(settings, dynamic_cfg);
            if(table.find(accessed_cache_line_count)==table.end())table[accessed_cache_line_count]=1;
            else table[accessed_cache_line_count]+=1;
        }
    }
    bool get_next_candidate(string *result, Settings *settings, set<string> &omit_list, CFG *dynamic_cfg)
    {
        unordered_map<string,uint64_t> candidate_counts;
        set<string> candidates;
        for(uint64_t i=0;i<all_misses.size();i++)
        {
            vector<uint64_t> ordered_candidate;
            all_misses[i]->get_candidates(ordered_candidate, settings, dynamic_cfg);
            if(ordered_candidate.size()<1)continue;
            candidates.clear();
            for(uint64_t k=0;k<settings->max_bbl_count;k++)
            {
                for(uint64_t j=0;j<ordered_candidate.size()-k;j++)
                {
                    vector<uint64_t> tmp;
                    for(uint64_t l=j;l<j+k+1;l++)
                    {
                        tmp.push_back(ordered_candidate[l]);
                    }
                    candidates.insert(vec_to_string(tmp));
                    tmp.clear();
                }
            }
            for(auto it: candidates)
            {
                if(omit_list.find(it)!=omit_list.end())continue;
                if(candidate_counts.find(it)==candidate_counts.end())candidate_counts[it]=1;
                else candidate_counts[it]+=1;
            }
        }
        vector<Candidate> sorted_candidates;
        for(auto it: candidate_counts)
        {
            string tmp = it.first;
            vector<uint64_t> all;
            string_to_vec(tmp, all);
            if(all.size()<1)continue;
            reverse(all.begin(),all.end());
            uint64_t converged_bbl_count = dynamic_cfg->get_bbl_execution_count(all[0]);
            for(uint64_t i=1;i<all.size();i++)converged_bbl_count*=dynamic_cfg->get_next_out(all[i-1],all[i]);
            sorted_candidates.push_back(Candidate(all, ((1.0*it.second)/(converged_bbl_count))));
        }
        sort(sorted_candidates.begin(),sorted_candidates.end());
        reverse(sorted_candidates.begin(),sorted_candidates.end());
        if(sorted_candidates.size()<1)return false;
        *result = vec_to_string(sorted_candidates[0].predecessor_bbl_address);
        return true;
    }
    uint64_t remove_covered_misses(string next, Settings *settings, CFG *dynamic_cfg)
    {
        vector<LBR_Sample *> tmp;
        uint64_t total=0;
        for(uint64_t i =0; i<all_misses.size();i++)
        {
            if(all_misses[i]->is_covered(next, settings, dynamic_cfg))
            {
                total+=1;
                all_misses[i]->clear();
            }
            else
            {
                tmp.push_back(all_misses[i]);
            }
        }
        all_misses.clear();
        for(uint64_t i = 0; i< tmp.size(); i++)
        {
            all_misses.push_back(tmp[i]);
        }
        return total;
    }
};

class Miss_Map
{
private:
    ofstream *bbl_mapping_out = nullptr;
public:
    unordered_map<uint64_t,LBR_List *> missed_pc_to_LBR_sample_list;
    vector<uint64_t> miss_profile;
    Miss_Map()
    {
        missed_pc_to_LBR_sample_list.clear();
    }
    ~Miss_Map()
    {
        if(bbl_mapping_out!=nullptr)bbl_mapping_out->close();
        missed_pc_to_LBR_sample_list.clear();
        miss_profile.clear();
    }
    void set_bbl_mapping_out(string file_name)
    {
        bbl_mapping_out = new ofstream(file_name.c_str());
    }
    uint64_t size()
    {
        return missed_pc_to_LBR_sample_list.size();
    }
    uint64_t total_misses()
    {
        uint64_t total = 0;
        for(auto it: missed_pc_to_LBR_sample_list)
        {
            total+=it.second->size();
        }
        return total;
    }
    void push(string line)
    {
        boost::trim_if(line,boost::is_any_of(","));
        vector<string> parsed;
        boost::split(parsed,line,boost::is_any_of(",\n"),boost::token_compress_on);
        if(parsed.size()<3)panic("ICache Miss LBR Sample contains less than 3 branch records");
        uint64_t missed_cache_line_id = string_to_u64(parsed[0]);
        uint64_t missed_bbl_address = string_to_u64(parsed[1]);
        miss_profile.push_back(missed_bbl_address);
        if((missed_bbl_address>>6)!=missed_cache_line_id)panic("In access based next line prefetching system a icache miss cannot happen other than the start of a new BBL");
        if(missed_pc_to_LBR_sample_list.find(missed_bbl_address)==missed_pc_to_LBR_sample_list.end())
        {
            missed_pc_to_LBR_sample_list[missed_bbl_address]=new LBR_List();
        }
        missed_pc_to_LBR_sample_list[missed_bbl_address]->push(parsed);
    }
    void measure_prefetch_window_cdf(Settings *settings, CFG *dynamic_cfg, ofstream &out)
    {
        unordered_map<uint64_t,uint64_t> table;
        for(auto it: missed_pc_to_LBR_sample_list)
        {
            it.second->measure_prefetch_window_cdf(table, settings, dynamic_cfg);
        }
        uint64_t total = 0;
        vector<uint64_t> sorted_counts;
        for(auto it: table)
        {
            total+=it.second;
            sorted_counts.push_back(it.first);
        }
        sort(sorted_counts.begin(),sorted_counts.end());
        double value = 0.0;
        for(uint64_t i=0;i<sorted_counts.size(); i++)
        {
            value+=table[sorted_counts[i]];
            out<<sorted_counts[i]<<" "<<((100.0*value)/total)<<endl;
        }
    }
    void generate_prefetch_locations(Settings *settings, CFG *dynamic_cfg, ofstream &out, bool multiline_enabled)
    {
        unordered_map<string,vector<uint64_t>> bbl_to_prefetch_targets;
        
        vector<pair<uint64_t,uint64_t>> sorted_miss_pcs;
        for(auto it: missed_pc_to_LBR_sample_list)
        {
            sorted_miss_pcs.push_back(make_pair(it.second->size(), it.first));
        }
        
        sort(sorted_miss_pcs.begin(), sorted_miss_pcs.end());
        reverse(sorted_miss_pcs.begin(), sorted_miss_pcs.end());


        uint64_t total_miss_sample_count = total_misses();
        
        if(multiline_enabled)
        {
            unordered_map<uint64_t,unordered_map<uint64_t,uint64_t>> neighbors;
            for(uint64_t i = 0; i < miss_profile.size(); i++)
            {
                uint64_t missed_pc = miss_profile[i];
                if(neighbors.find(missed_pc)==neighbors.end())
                {
                    neighbors[missed_pc]=unordered_map<uint64_t,uint64_t>();
                }
                for(uint64_t j = 1; j<= settings->max_step_ahead; j++)
                {
                    if(i+j>=miss_profile.size())break;
                    uint64_t neighbor_missed_pc = miss_profile[i+j];
                    int prefetch_length = measure_prefetch_length(missed_pc, neighbor_missed_pc);
                    if(my_abs(prefetch_length)<=settings->max_prefetch_length)
                    {
                        if(prefetch_length==0)continue;
                        if(neighbors[missed_pc].find(prefetch_length)==neighbors[missed_pc].end())
                        {
                            neighbors[missed_pc][neighbor_missed_pc]=1;
                        }
                        else
                        {
                            neighbors[missed_pc][neighbor_missed_pc]+=1;
                        }
                    }
                }
            }
            vector<prefetch_benefit> multiline_sorted_missed_pcs;
            map<uint64_t,vector<uint64_t>> prefetch_lists;
            for(int i=0;i<sorted_miss_pcs.size();i++)
            {
                uint64_t current_pc = sorted_miss_pcs[i].second;
                uint64_t prefetch_count = sorted_miss_pcs[i].first;
                uint64_t covered_miss_count = prefetch_count;
                prefetch_lists[current_pc]=vector<uint64_t>();
                prefetch_lists[current_pc].push_back(current_pc);
                for(auto it: neighbors[current_pc])
                {
                    double ratio_percentage = (1.0 * it.second) / prefetch_count;
                    if(ratio_percentage>=settings->min_ratio_percentage)
                    {
                        covered_miss_count+=it.second;
                        prefetch_lists[current_pc].push_back(it.first);
                    }
                }
                multiline_sorted_missed_pcs.push_back(prefetch_benefit(current_pc, prefetch_count, covered_miss_count));
            }
            sort(multiline_sorted_missed_pcs.begin(), multiline_sorted_missed_pcs.end());
            reverse(multiline_sorted_missed_pcs.begin(), multiline_sorted_missed_pcs.end());
            
            uint64_t permissible_prefetch_count = settings->dyn_ins_count_inc_ratio * settings->total_dyn_ins_count;
            uint64_t current_running_count = 0;
            bool prefetch_count_finished = false;
            bbl_to_prefetch_targets.clear();
            uint64_t total_covered_miss_count=0;
            for(int i=0;i<multiline_sorted_missed_pcs.size();i++)
            {
                if(prefetch_count_finished)break;
                uint64_t missed_pc = multiline_sorted_missed_pcs[i].addr;
                uint64_t missed_count_for_this_pc = multiline_sorted_missed_pcs[i].prefetch_count;
                string prefetch_candidate;
                set<string> omit_list;
                while(missed_count_for_this_pc>0)
                {
                    bool candidate_found = missed_pc_to_LBR_sample_list[missed_pc]->get_next_candidate(&prefetch_candidate, settings, omit_list, dynamic_cfg);
                    if(!candidate_found)break;
                    
                    omit_list.insert(prefetch_candidate);
                    
                    //
                    vector<uint64_t> all;
                    string_to_vec(prefetch_candidate, all);
                    if(all.size()<1)continue;
                    reverse(all.begin(),all.end());
                    uint64_t converged_bbl_count = dynamic_cfg->get_bbl_execution_count(all[0]);
                    for(uint64_t k=1;k<all.size();k++)converged_bbl_count*=dynamic_cfg->get_next_out(all[k-1],all[k]);
                    //
                    
                    uint64_t inserted_prefetch_count = converged_bbl_count;//dynamic_cfg->get_bbl_execution_count(prefetch_candidate);
                    
                    if( bbl_to_prefetch_targets.find(prefetch_candidate)==bbl_to_prefetch_targets.end() )bbl_to_prefetch_targets[prefetch_candidate]=vector<uint64_t>();
                    else if(bbl_to_prefetch_targets[prefetch_candidate].size()>0)
                    {
                        //Check if all the previous prefetch targets and this target are nearby
                        //if yes, okay
                        //else, not okay
                        uint64_t first_inserted_target = bbl_to_prefetch_targets[prefetch_candidate][0];
                        int prefetch_length = measure_prefetch_length(missed_pc, first_inserted_target);
                        if(my_abs(prefetch_length)<=settings->max_prefetch_length)
                        {
                            //this is good to insert
                            //and there is no extra static or dynamic instruction insert
                            //benfit
                            inserted_prefetch_count = 0;
                        }
                        else
                        {
                            if(settings->insert_as_many_as_possible == 0)continue;
                        }
                    }
                    
                    if(settings->insert_as_many_as_possible == 0)
                    {
                        if(current_running_count+inserted_prefetch_count >= permissible_prefetch_count+100)
                        {
                            if(bbl_to_prefetch_targets[prefetch_candidate].size()==0)bbl_to_prefetch_targets.erase(prefetch_candidate);
                            continue;
                        }
                    }
                    
                    uint64_t covered_miss_count = 0;/*missed_pc_to_LBR_sample_list[missed_pc]->remove_covered_misses(prefetch_candidate, settings, dynamic_cfg);*/
                    /*for(uint64_t j=0;j<prefetch_lists[missed_pc].size();j++)
                    {
                        bbl_to_prefetch_targets[prefetch_candidate].push_back(prefetch_lists[missed_pc][j]);
                        covered_miss_count+=missed_pc_to_LBR_sample_list[prefetch_lists[missed_pc][j]]->remove_covered_misses(prefetch_candidate, settings, dynamic_cfg);
                        if(j==0)missed_count_for_this_pc-=covered_miss_count;
                    }*/
                    bbl_to_prefetch_targets[prefetch_candidate].push_back(missed_pc);
                    covered_miss_count+=missed_pc_to_LBR_sample_list[missed_pc]->remove_covered_misses(prefetch_candidate, settings, dynamic_cfg);
                    missed_count_for_this_pc-=covered_miss_count;
                    
                    current_running_count += inserted_prefetch_count;
                    total_covered_miss_count+=covered_miss_count;
                    
                    if(settings->insert_as_many_as_possible == 0)
                    {
                        if(current_running_count>=permissible_prefetch_count)
                        {
                            prefetch_count_finished=true;
                            break;
                        }
                    }
                }
            }
            uint64_t static_prefetch_count = 0;

            uint8_t size_of_prefetch_inst = (7+((settings->max_prefetch_length*1)/8));
            unordered_map<uint64_t,uint64_t> prev_to_new;
            //dynamic_cfg->print_modified_bbl_mappings(bbl_to_prefetch_targets, size_of_prefetch_inst, bbl_mapping_out, prev_to_new);

            vector<string> tmp;
            for(auto it:bbl_to_prefetch_targets)
            {
                tmp.push_back(it.first);
            }
            sort(tmp.begin(), tmp.end());
            for(string bbl_addr: tmp)
            {
                auto it = bbl_to_prefetch_targets[bbl_addr];
                if(it.size()<1)continue;
                //out<<prev_to_new[bbl_addr]<<" "<<it.size();
                /*for(auto set_it: it)
                {
                    out<<" "<<prev_to_new[set_it];
                }*/
                //out<<endl;
                static_prefetch_count++;
            }
            tmp.clear();
            prev_to_new.clear();
            cerr<<total_covered_miss_count<<","<<(100.0*total_covered_miss_count)/total_miss_sample_count<<","<<current_running_count<<","<<(100.0*current_running_count)/settings->total_dyn_ins_count<<endl;
            uint64_t total_count = dynamic_cfg->get_static_count();
            cerr<<"Static counts: "<<static_prefetch_count<<","<<total_count<<","<<(100.0*static_prefetch_count)/total_count<<endl;
            uint64_t total_size = dynamic_cfg->get_static_size();
            double static_prefetch_size = static_prefetch_count*size_of_prefetch_inst;
            cerr<<"Static size: "<<static_prefetch_size<<","<<total_size<<","<<(100.0*static_prefetch_size)/total_size<<endl;
        }
        else
        {
            uint64_t permissible_prefetch_count = settings->dyn_ins_count_inc_ratio * settings->total_dyn_ins_count;
            uint64_t current_running_count = 0;
            bool prefetch_count_finished = false;
            bbl_to_prefetch_targets.clear();
            uint64_t total_covered_miss_count=0;
            for(int i=0;i<sorted_miss_pcs.size();i++)
            {
                if(prefetch_count_finished)break;
                uint64_t missed_pc = sorted_miss_pcs[i].second;
                uint64_t missed_count_for_this_pc = sorted_miss_pcs[i].first;
                string prefetch_candidate;
                set<string> omit_list;
                while(missed_count_for_this_pc>0)
                {
                    bool candidate_found = missed_pc_to_LBR_sample_list[missed_pc]->get_next_candidate(&prefetch_candidate, settings, omit_list, dynamic_cfg);
                    if(!candidate_found)break;
                    
                    omit_list.insert(prefetch_candidate);
                    
                    //bool candidate_valid = dynamic_cfg->is_valid_candidate(prefetch_candidate, missed_pc);
                    //if(!candidate_valid)continue;

                    //
                    vector<uint64_t> all;
                    string_to_vec(prefetch_candidate, all);
                    if(all.size()<1)continue;
                    reverse(all.begin(),all.end());
                    uint64_t converged_bbl_count = dynamic_cfg->get_bbl_execution_count(all[0]);
                    for(uint64_t k=1;k<all.size();k++)converged_bbl_count*=dynamic_cfg->get_next_out(all[k-1],all[k]);
                    //
                    
                    uint64_t inserted_prefetch_count = converged_bbl_count;//dynamic_cfg->get_bbl_execution_count(prefetch_candidate);
                    if(current_running_count+inserted_prefetch_count >= permissible_prefetch_count+100)continue;
                    
                    if( bbl_to_prefetch_targets.find(prefetch_candidate)==bbl_to_prefetch_targets.end() )bbl_to_prefetch_targets[prefetch_candidate]=vector<uint64_t>();
                    else continue;
                    
                    bbl_to_prefetch_targets[prefetch_candidate].push_back(missed_pc);
                    current_running_count += inserted_prefetch_count;
                    uint64_t covered_miss_count = missed_pc_to_LBR_sample_list[missed_pc]->remove_covered_misses(prefetch_candidate, settings, dynamic_cfg);
                    missed_count_for_this_pc-=covered_miss_count;
                    total_covered_miss_count+=covered_miss_count;
                    
                    if(current_running_count>=permissible_prefetch_count)
                    {
                        prefetch_count_finished=true;
                        break;
                    }
                }
            }
            uint64_t static_prefetch_count = 0;

            unordered_map<uint64_t,uint64_t> prev_to_new;
            //dynamic_cfg->print_modified_bbl_mappings(bbl_to_prefetch_targets, 7, bbl_mapping_out, prev_to_new);
            for(auto it:bbl_to_prefetch_targets)
            {
                if(it.second.size()<1)continue;
                //out<<prev_to_new[it.first]<<" "<<it.second.size();
                /*for(auto set_it: it.second)
                {
                    out<<" "<<prev_to_new[set_it];
                }*/
                //out<<endl;
                static_prefetch_count++;
            }
            cerr<<total_covered_miss_count<<","<<(100.0*total_covered_miss_count)/total_miss_sample_count<<","<<current_running_count<<","<<(100.0*current_running_count)/settings->total_dyn_ins_count<<endl;
            uint64_t total_count = dynamic_cfg->get_static_count();
            cerr<<"Static counts: "<<static_prefetch_count<<","<<total_count<<","<<(100.0*static_prefetch_count)/total_count<<endl;
            uint64_t total_size = dynamic_cfg->get_static_size();
            cerr<<"Static size: "<<static_prefetch_count*7<<","<<total_size<<","<<(700.0*static_prefetch_count)/total_size<<endl;
        }
        
    }
};

class Miss_Profile
{
public:
    Miss_Map profile;
    Miss_Profile(MyConfigWrapper &conf, Settings &settings, CFG &dynamic_cfg, bool measure_prefetch_window_cdf=false)
    {
        string miss_log_path = conf.lookup("miss_log", "/tmp/");
        if(miss_log_path == "/tmp/")
        {
            panic("ICache Miss Log file is not present in the config");
        }
        vector<string> raw_icache_miss_lbr_sample_data;
        read_full_file(miss_log_path, raw_icache_miss_lbr_sample_data);
        for(uint64_t i = 0; i< raw_icache_miss_lbr_sample_data.size();i++)
        {
            string line = raw_icache_miss_lbr_sample_data[i];
            profile.push(line);
        }
        string output_name = "";
        output_name+=to_string(settings.multiline_mode)+".";
        output_name+=to_string(settings.min_distance_cycle_count)+".";
        output_name+=to_string(settings.max_distance_cycle_count)+".";
        output_name+=to_string(int(settings.dyn_ins_count_inc_ratio*100))+".";
        output_name+=to_string(int(settings.fan_out_ratio*100))+".";
        output_name+="txt";
        ofstream prefetch_out(output_name.c_str());
        profile.set_bbl_mapping_out(("bbl."+output_name));
        if(measure_prefetch_window_cdf)
        {
            ofstream prefetch_window_accessed_cache_line_count_cdf("prefetch_window_accessed_cache_line_count_cdf.txt");
            profile.measure_prefetch_window_cdf(&settings,&dynamic_cfg, prefetch_window_accessed_cache_line_count_cdf);
        }
        if(settings.multiline_mode == 1)//ASMDB
        {
            profile.generate_prefetch_locations(&settings,&dynamic_cfg, prefetch_out, false);
        }
        else if(settings.multiline_mode == 2)//OUR
        {
            profile.generate_prefetch_locations(&settings,&dynamic_cfg, prefetch_out, true);
        }
    }
};

#endif //LBR_MISS_SAMPLE_H_
