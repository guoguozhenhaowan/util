#ifndef FUZZY_H_
#define FUZZY_H_

#include <iostream>
#include <climits>
#include <string>
#include <vector>
#include <array>

namespace fuzzy{
    class bitap{
        const std::string& text;
        const std::string& pattern;
        int mismatch;
        std::vector<int> results;
        public:
            bitap(const std::string& t, const std::string& p) : text(t), pattern(p){
                mismatch = 0;
            };
            bitap(const std::string& t, const std::string& p, int k) : text(t), pattern(p), mismatch(k){};
            ~bitap(){};
            
            const std::vector<int>& getResults(){
                return results;
            }
            
            int countMatch(){
                return results.size();
            }

            void search(){
                if(mismatch == 0) bitap_bitwise_search();
                else  bitap_fuzzy_bitwise_search();
            }

            void bitap_bitwise_search(){
                if(pattern.length() == 0){
                    return;
                }
                if(pattern.length() > 31){
                    std::cerr << "The pattern is too long" << std::endl;
                    return;
                }
                unsigned long R = 0;
                std::array<unsigned long, CHAR_MAX + 1> pattern_mask;
                pattern_mask.fill(~0);
                for(size_t i = 0; i < pattern.length(); ++i){
                    pattern_mask[pattern[i]] &= ~(1UL << i);
                }
                for(size_t i = 0; i < text.length(); ++i)
                {
                    R |= pattern_mask[text[i]];
                    R <<= 1;
                    if(0 == (R & (1UL << pattern.length()))){
                        results.push_back(i - pattern.length() + 1);
                    }
                }
            }
            
            void bitap_fuzzy_bitwise_search()
            {
                if(pattern.length() == 0){
                    return;
                }
                if(pattern.length() > 31){
                    std::cerr << "The pattern is too long" << std::endl;
                    return;
                }
                std::vector<unsigned long> R;
                R.resize(mismatch + 1);
                std::fill(R.begin(), R.end(), ~1);
                std::array<unsigned long, CHAR_MAX + 1> pattern_mask;
                pattern_mask.fill(~0);
                for(size_t i = 0; i < pattern.length(); ++i){
                    pattern_mask[pattern[i]] &= ~(1UL << i);
                }
                for(size_t i = 0; i < text.length(); ++i)
                {
                    unsigned long old_Rd1 = R[0];
                    R[0] |= pattern_mask[text[i]];
                    R[0] <<= 1;
            
                    for(int d = 1; d <= mismatch; ++d){
                        unsigned long tmp = R[d];
                        R[d] = (old_Rd1 & (R[d] | pattern_mask[text[i]])) << 1;
                        old_Rd1 = tmp;
                    }
                    if(0 == (R[mismatch] & (1UL << pattern.length()))){
                        results.push_back(i - pattern.length() + 1);
                    }
                }
            }
        };
};
#endif
