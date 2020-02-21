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

class Candidate
{
public:
    uint64_t predecessor_bbl_address;
    double covered_miss_ratio_to_predecessor_count;
    Candidate(uint64_t addr, double c_ratio)
    {
        predecessor_bbl_address = addr;
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
    prefetch_benefit(uint64_t a, uint64_t p_count, uint64_t cm_count) : addr(a), prefetch_count(p_count), covered_miss_count(cm_count) {}
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

struct benefit
{
    uint64_t missed_bbl_address;
    uint64_t predecessor_bbl_address;
    uint64_t covered_miss_counts;
    double fan_in;
    benefit (uint64_t m, uint64_t p, uint64_t c, double f) : missed_bbl_address(m), predecessor_bbl_address(p), covered_miss_counts(c), fan_in(f) {}
    bool operator < (const benefit& right) const
    {
        if(right.fan_in>fan_in)return true;
        else if(fan_in>right.fan_in)return false;
        else
        {
            if(right.covered_miss_counts>covered_miss_counts)return true;
            return false;
        }
    }
};

struct multi_bbl_benefit
{
    uint64_t missed_bbl_address;
    uint64_t predecessor_bbl_address;
    uint64_t predicate_bbl_address;
    uint64_t covered_miss_counts;
    double fan_in;
    uint64_t dynamic_prefetch_counts;
    multi_bbl_benefit (uint64_t m, uint64_t predecessor, uint64_t predicate, uint64_t c, double f, uint64_t d_p_counts) : missed_bbl_address(m), predecessor_bbl_address(predecessor), predicate_bbl_address(predicate), covered_miss_counts(c), fan_in(f), dynamic_prefetch_counts(d_p_counts) {}
    bool operator < (const multi_bbl_benefit& right) const
    {
        if(right.fan_in>fan_in)return true;
        else if(fan_in>right.fan_in)return false;
        else
        {
            if(right.covered_miss_counts>covered_miss_counts)return true;
            return false;
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
    void get_candidates(set<uint64_t> &candidates, Settings *settings, CFG *dynamic_cfg)
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
                candidates.insert(current_lbr[i].first);
            }
        }
    }
    void get_correlated_candidates(uint64_t top_candidate, set<uint64_t> &candidates, Settings *settings, CFG *dynamic_cfg)
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
                //candidates.insert(current_lbr[i].first);
                if(current_lbr[i].first == top_candidate)
                {
                    //insert previous 8 BBLs on the candidate list
                    for(uint64_t j=i+1; (j<current_lbr.size() && j<i+9); j++)
                    {
                        if(current_lbr[j].first == top_candidate)continue;
                        candidates.insert(current_lbr[j].first);
                    }
                    return;
                }
            }
        }
    }
    bool is_covered(uint64_t predecessor_bbl_address, Settings *settings, CFG *dynamic_cfg)
    {
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
                if(current_distance>max_distance)break;
                if(current_lbr[i].first==predecessor_bbl_address)return true;
            }
        }
        return false;
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
    bool get_next_candidate(uint64_t *result, Settings *settings, set<uint64_t> &omit_list, CFG *dynamic_cfg)
    {
        unordered_map<uint64_t,uint64_t> candidate_counts;
        set<uint64_t> candidates;
        for(uint64_t i=0;i<all_misses.size();i++)
        {
            all_misses[i]->get_candidates(candidates, settings, dynamic_cfg);
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
            sorted_candidates.push_back(Candidate(it.first, ((1.0*it.second)/(dynamic_cfg->get_bbl_execution_count(it.first)))));
        }
        sort(sorted_candidates.begin(),sorted_candidates.end());
        reverse(sorted_candidates.begin(),sorted_candidates.end());
        if(sorted_candidates.size()<1)return false;
        *result = sorted_candidates[0].predecessor_bbl_address;
        return true;
    }
    bool get_top_candidate(vector<uint64_t> &result, Settings *settings, CFG *dynamic_cfg)
    {
        result.clear();
        unordered_map<uint64_t,uint64_t> candidate_counts;
        set<uint64_t> candidates;
        for(uint64_t i=0;i<all_misses.size();i++)
        {
            all_misses[i]->get_candidates(candidates, settings, dynamic_cfg);
            for(auto it: candidates)
            {
                if(candidate_counts.find(it)==candidate_counts.end())candidate_counts[it]=1;
                else candidate_counts[it]+=1;
            }
        }
        vector<Candidate> sorted_candidates;
        for(auto it: candidate_counts)
        {
            sorted_candidates.push_back(Candidate(it.first, it.second));
        }
        sort(sorted_candidates.begin(),sorted_candidates.end());
        reverse(sorted_candidates.begin(),sorted_candidates.end());
        if(sorted_candidates.size()<1)return false;
        for(int i=0;i<sorted_candidates.size();i++)result.push_back(sorted_candidates[i].predecessor_bbl_address);
        sorted_candidates.clear();
        candidate_counts.clear();
        return true;
    }
    void get_candidate_counts(unordered_map<uint64_t,uint64_t> &candidate_counts, Settings *settings, set<uint64_t> &omit_list, CFG *dynamic_cfg)
    {
        set<uint64_t> candidates;
        for(uint64_t i=0;i<all_misses.size();i++)
        {
            all_misses[i]->get_candidates(candidates, settings, dynamic_cfg);
            for(auto it: candidates)
            {
                if(omit_list.find(it)!=omit_list.end())continue;
                if(candidate_counts.find(it)==candidate_counts.end())candidate_counts[it]=1;
                else candidate_counts[it]+=1;
            }
        }
    }
    uint64_t remove_covered_misses(uint64_t next, Settings *settings, CFG *dynamic_cfg)
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
    bool get_top_correlated_candidate(uint64_t *result, uint64_t *covered_miss_count, uint64_t next, Settings *settings, CFG *dynamic_cfg)
    {
        unordered_map<uint64_t,uint64_t> candidate_counts;
        set<uint64_t> candidates;
        for(uint64_t i=0;i<all_misses.size();i++)
        {
            all_misses[i]->get_correlated_candidates(next, candidates, settings, dynamic_cfg);
            for(auto it: candidates)
            {
                if(it==next)continue;
                if(candidate_counts.find(it)==candidate_counts.end())candidate_counts[it]=1;
                else candidate_counts[it]+=1;
            }
        }
        vector<Candidate> sorted_candidates;
        for(auto it: candidate_counts)
        {
            sorted_candidates.push_back(Candidate(it.first, it.second));
        }
        sort(sorted_candidates.begin(),sorted_candidates.end());
        reverse(sorted_candidates.begin(),sorted_candidates.end());
        if(sorted_candidates.size()<1)return false;
        *result = sorted_candidates[0].predecessor_bbl_address;
        *covered_miss_count = sorted_candidates[0].covered_miss_ratio_to_predecessor_count;
        sorted_candidates.clear();
        candidate_counts.clear();
        return true;
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
    void push(string line, CFG &dynamic_cfg)
    {
        boost::trim_if(line,boost::is_any_of(","));
        vector<string> parsed;
        boost::split(parsed,line,boost::is_any_of(",\n"),boost::token_compress_on);
        if(parsed.size()<3)panic("ICache Miss LBR Sample contains less than 3 branch records");
        uint64_t missed_cache_line_id = string_to_u64(parsed[0]);
        uint64_t missed_bbl_address = string_to_u64(parsed[1]);
        if(dynamic_cfg.is_self_modifying(missed_bbl_address))return;
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
    void generate_prefetch_locations(Settings *settings, CFG *dynamic_cfg, ofstream &out)
    {
        unordered_map<uint64_t,vector<uint64_t>> bbl_to_prefetch_targets;
        
        vector<pair<uint64_t,uint64_t>> sorted_miss_pcs;
        for(auto it: missed_pc_to_LBR_sample_list)
        {
            sorted_miss_pcs.push_back(make_pair(it.second->size(), it.first));
        }
        
        sort(sorted_miss_pcs.begin(), sorted_miss_pcs.end());
        reverse(sorted_miss_pcs.begin(), sorted_miss_pcs.end());


        uint64_t total_miss_sample_count = total_misses();
        
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
        set<uint64_t> omit_list;
        vector<benefit> candidates;
        for(int i =0; i<sorted_miss_pcs.size(); i++)
        {
            uint64_t missed_pc = sorted_miss_pcs[i].second;
            unordered_map<uint64_t,uint64_t> candidate_counts;
            missed_pc_to_LBR_sample_list[missed_pc]->get_candidate_counts(candidate_counts,settings,omit_list,dynamic_cfg);
            for(auto it : candidate_counts)
            {
                uint64_t predecessor_bbl_address = it.first;
                if(dynamic_cfg->get_fan_out(predecessor_bbl_address,missed_pc) < 0.99)continue;
                uint64_t miss_count = it.second;
                double fan_in = dynamic_cfg->get_fan_out(predecessor_bbl_address, missed_pc);
                candidates.push_back(benefit(missed_pc,predecessor_bbl_address,miss_count,fan_in));
            }
            candidate_counts.clear();
        }
        sort(candidates.begin(), candidates.end());
        reverse(candidates.begin(),candidates.end());
        for(int i = 0; i< candidates.size(); i++)
        {   
            benefit &best_candidate = candidates[i];
            uint64_t prefetch_candidate = best_candidate.predecessor_bbl_address;
            uint64_t missed_pc = best_candidate.missed_bbl_address;
            
            uint64_t inserted_prefetch_count = dynamic_cfg->get_bbl_execution_count(prefetch_candidate);

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
            }
            
            uint64_t covered_miss_count = 0;
            bbl_to_prefetch_targets[prefetch_candidate].push_back(missed_pc);
            covered_miss_count+=missed_pc_to_LBR_sample_list[missed_pc]->remove_covered_misses(prefetch_candidate, settings, dynamic_cfg);
            
            current_running_count += inserted_prefetch_count;
            total_covered_miss_count+=covered_miss_count;

            //if(current_running_count >= permissible_prefetch_count)break;
        }
        candidates.clear();

        vector<multi_bbl_benefit> multi_bbl_candidates;

        for(int i =0; i<sorted_miss_pcs.size(); i++)
        {
            uint64_t missed_pc = sorted_miss_pcs[i].second;
            /*uint64_t still_missing = missed_pc_to_LBR_sample_list[missed_pc]->size();
            out<<missed_pc<<","<<still_missing<<endl;*/
            vector<uint64_t> top_candidate;
            bool candidate_found = missed_pc_to_LBR_sample_list[missed_pc]->get_top_candidate(top_candidate,settings,dynamic_cfg);
            if(candidate_found == false)continue;
            for(int j =0;j<top_candidate.size();j++)
            {    
                uint64_t next_candidate;
                uint64_t covered_miss_count;
                uint64_t dynamic_prefetch_counts;
                bool second_candidate_found = missed_pc_to_LBR_sample_list[missed_pc]->get_top_correlated_candidate(&next_candidate,&covered_miss_count,top_candidate[j],settings,dynamic_cfg);
                if(second_candidate_found == false)continue;
                double fan_in_combined = dynamic_cfg->get_fan_in_for_multiple_prior_bbls(missed_pc,top_candidate[j],next_candidate,&dynamic_prefetch_counts);
                if(fan_in_combined>=0.99)
                {
                    multi_bbl_candidates.push_back(multi_bbl_benefit(missed_pc,top_candidate[j],next_candidate,covered_miss_count,fan_in_combined,dynamic_prefetch_counts));
                }
            }
        }
        sort(multi_bbl_candidates.begin(),multi_bbl_candidates.end());
        reverse(multi_bbl_candidates.begin(),multi_bbl_candidates.end());

        for(int i = 0; i< multi_bbl_candidates.size(); i++)
        {
            multi_bbl_benefit &best_candidate = multi_bbl_candidates[i];
            uint64_t prefetch_candidate_predecessor = best_candidate.predecessor_bbl_address;
            uint64_t prefetch_candidate_predicate = best_candidate.predicate_bbl_address;
            uint64_t missed_pc = best_candidate.missed_bbl_address;
            uint64_t covered_miss_count = best_candidate.covered_miss_counts;
            uint64_t inserted_prefetch_count = best_candidate.dynamic_prefetch_counts;

            current_running_count += inserted_prefetch_count;
            total_covered_miss_count+=covered_miss_count;

            out<<"CS,"<<prefetch_candidate_predicate<<","<<prefetch_candidate_predecessor<<","<<missed_pc<<endl;
        }

        multi_bbl_candidates.clear();


        uint64_t static_prefetch_count = 0;

        uint8_t size_of_prefetch_inst = (7+((settings->max_prefetch_length*1)/8));
        cerr<<total_covered_miss_count<<","<<(100.0*total_covered_miss_count)/total_miss_sample_count<<","<<current_running_count<<","<<(100.0*current_running_count)/settings->total_dyn_ins_count<<endl;
        uint64_t total_count = dynamic_cfg->get_static_count();
        cerr<<"Static counts: "<<static_prefetch_count<<","<<total_count<<","<<(100.0*static_prefetch_count)/total_count<<endl;
        uint64_t total_size = dynamic_cfg->get_static_size();
        double static_prefetch_size = static_prefetch_count*size_of_prefetch_inst;
        cerr<<"Static size: "<<static_prefetch_size<<","<<total_size<<","<<(100.0*static_prefetch_size)/total_size<<endl;
    }
};

class Miss_Profile
{
public:
    Miss_Map profile;
    Miss_Profile(MyConfigWrapper &conf, Settings &settings, CFG &dynamic_cfg)
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
            profile.push(line, dynamic_cfg);
        }
        string output_name = "";
        output_name+=to_string(settings.multiline_mode)+".";
        output_name+=to_string(settings.min_distance_cycle_count)+".";
        output_name+=to_string(settings.max_distance_cycle_count)+".";
        output_name+=to_string(int(settings.dyn_ins_count_inc_ratio*100))+".";
        output_name+=to_string(int(settings.fan_out_ratio*100))+".";
        output_name+="txt";
        ofstream prefetch_out(output_name.c_str());
        profile.generate_prefetch_locations(&settings,&dynamic_cfg, prefetch_out);
    }
};

#endif //LBR_MISS_SAMPLE_H_
