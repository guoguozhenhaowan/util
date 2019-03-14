#ifndef FUZZY_H_
#define FUZZY_H_

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <bitset>
#include <limits>

template <typename T>
/** class to operate k-mismatch through bitap algorithm */
class bitap{
    const std::string& text;    ///< string to be matched against
    const std::string& pattern; ///< pattern string
    int mismatch;               ///< maximum mismatches allowed
    std::vector<int> results;   ///< matched index results
    public:
        /** Construct a bitap object for exact match 
         * @param t string to be matched against
         * @param p pattern string
         */
        bitap(const std::string& t, const std::string& p) : text(t), pattern(p){
            mismatch = 0;
        };
        
        /** Construct a bitap object for at most k mismatches
         * @param t string to be matched against
         * @param p pattern string
         * @param k maximum mismatches allowed
         */
        bitap(const std::string& t, const std::string& p, int k) : text(t), pattern(p), mismatch(k){};
        ~bitap(){};
        
        /** Get match results
         * @return vector of index match happens
         */
        const std::vector<int>& getResults(){
            return results;
        }
       
        /** Get match numbers
         * @return number of matches 
         */ 
        int countMatch(){
            return results.size();
        }

        /** Search the pattern against the text
         * @param count max allowed mismatches
         */
        void search(int count){
            if(mismatch == 0) bitap_bitwise_search(count);
            else  bitap_fuzzy_bitwise_search(count);
        }

        /** exact bitap match search
         * @param count max allowed mismatches
         */
        void bitap_bitwise_search(int count){
            int get_match = 0;
            if(pattern.length() == 0){
                return;
            }
            if(pattern.length() > std::numeric_limits<T>::digits){
                std::cerr << "The pattern is too long" << std::endl;
                return;
            }
            T R = 0;
            std::array<T, std::numeric_limits<char>::max() + 1>  pattern_mask;
            pattern_mask.fill(~0);
            
            for(size_t i = 0; i < pattern.length(); ++i){
                pattern_mask[pattern[i]] &= ~(1ULL << i);
            }
            
            for(size_t i = 0; i < text.length(); ++i)
            {
                R |= pattern_mask[text[i]];
                R <<= 1;
                if((0 == (R & (1ULL << pattern.length()))) && get_match < count){
                    results.push_back(i - pattern.length() + 1);
                    ++get_match;
                }
            }
        }
        
        /** k-mismatch search
         * @param count max allowed mismatches
         */ 
        void bitap_fuzzy_bitwise_search(int count)
        {
            int get_match = 0;
            if(pattern.length() == 0){
                return;
            }
            if(pattern.length() > std::numeric_limits<T>::digits){
                std::cerr << "The pattern is too long" << std::endl;
                return;
            }
            std::vector<T> R;
            R.resize(mismatch + 1);
            std::fill(R.begin(), R.end(), ~1ULL);
            std::array<T, std::numeric_limits<char>::max() + 1> pattern_mask;
            pattern_mask.fill(~0);
            for(size_t i = 0; i < pattern.length(); ++i){
                pattern_mask[pattern[i]] &= ~(1ULL << i);
            }

            for(size_t i = 0; i < text.length(); ++i)
            {
                T old_Rd1 = R[0];
                R[0] |= pattern_mask[text[i]];
                R[0] <<= 1;
                for(int d = 1; d <= mismatch; ++d){
                    T tmp = R[d];
                    R[d] = (old_Rd1 & (R[d] | pattern_mask[text[i]])) << 1;
                    old_Rd1 = tmp;
                }
                if((0 == (R[mismatch] & (1ULL << pattern.length()))) && get_match < count){
                    results.push_back(i - pattern.length() + 1);
                    ++get_match;
                }
            }
        }
};
#endif
