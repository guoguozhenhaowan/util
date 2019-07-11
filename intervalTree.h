#ifndef INTERVAL_TREE
#define INTERVAL_TREE

#include <vector>
#include <iostream>
#include <algorithm>

/** class to store an interval value: [start, stop] */
template<class T, typename K>
class Interval{
    public:
        K start; ///< interval starting coordinate
        K stop;  ///< interval stopping coordinate
        T value; ///< interval related object(name .etc)

        /* Interval constructor
         * @param s interval starting coordinate
         * @param e interval stopping coordinate
         * @param v interval related object(name .etc)
         */ 
        Interval(K s, K e, const T& v) : start(s), stop(e), value(v){}
        
        /** get interval starting coordinate
         * @param i Interval object
         * @return interval starting coordinate
         */
        static K intervalStart(const Interval<T, K>& i){
            return i.start;
        }

        /** get interval stopping coordinate
         * @param interval starting coordinate
         * @return interval stopping coordinate
         */
        static K intervalStop(const Interval<T, K>& i){
            return i.stop;
        }

        /** output an interval to ostream
         * @param out reference of ostream
         * @param i reference of Interval
         * @return reference of ostream
         */
        friend std::ostream& operator<<(std::ostream& out, const Interval<T, K>& i){
            out << "Interval [" << i.start << ", " << i.stop << "] : " << i.value;
            return out;
        }
};

/** class to sort an Interval object by starting coordinate */
template<class T, typename K>
class IntervalStartSorter{
    public:
        bool operator()(const Interval<T, K>& a, const Interval<T, K>& b){
            return a.start < b.start;
        }
};

/** IntervalTree class */
template<class T, typename K>
class IntervalTree{
    public:
        std::vector<Interval<T, K>> intervals; ///< intervals container with intervals stop >= center && start <= center
        IntervalTree<T, K>* left;              ///< left IntervalTree with intervals stop < center
        IntervalTree<T, K>* right;             ///< right IntervalTree with intervals start > center
        K center;                              ///< center value coordinate

        /** IntervalTree constructor */
        IntervalTree() : left(NULL), right(NULL), center(0) {}
        
        /** IntervalTree destructor */
        ~IntervalTree(){
            if(left){
                delete left;
            }
            if(right){
                delete right;
            }
        }

    private:
        /** copy a tree
         * @param orig reference of IntervalTree object
         * @return pointer to a new IntervalTree copied
         */
        IntervalTree* copyTree(const IntervalTree& orig){
            return new IntervalTree(orig);
        }
        
    public:
        /** copy constructor of IntervalTree
         * @param other reference of IntervalTree object
         */
        IntervalTree(const IntervalTree& other){
            intervals = other.intervals;
            left = (other.left ? copyTree(*other.left) : NULL);
            right = (other.right ? copyTree(*other.right) : NULL);
            center = other.center;
        }

        /** assignment operator of IntervalTree
         * @param other reference of IntervalTree object
         * @return reference of copied IntervalTree object
         */
        IntervalTree& operator=(const IntervalTree& other){
            intervals = other.intervals;
            left = (other.left ? copyTree(*other.left) : NULL);
            right = (other.right ? copyTree(*other.right) : NULL);
            center = other.center;
            return *this;
        }

        /** construct a IntervalTree with predefined parameters
         * @param ivals intervals used to construct IntervalTree
         * @param depth max depth of the IntervalTree
         * @param leftextent smallest value of starting coordinates of all ivals
         * @param rightextent greatest value of ending coordinates of all ivals
         * @minbucket minimum min number of intervals an node allowed to store
         */
        IntervalTree(std::vector<Interval<T, K>>& ivals, int32_t depth = 16, K leftextent = 0, K rightextent = 0, int32_t minbucket = 64){
            left = NULL;
            right = NULL;
            --depth;
            IntervalStartSorter<T, K> intervalStartSorter;
            if(depth == 0 || (ivals.size() < minbucket)){
                std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
                intervals = ivals;
            }else{
                if(leftextent == 0 && rightextent == 0){
                    std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
                }
                K leftp = 0;
                K rightp = 0;
                K centerp = 0;
                if(leftextent || rightextent){
                    leftp = leftextent;
                    rightp = rightextent;
                }else{
                    leftp = ivals.front().start;
                    std::vector<K> stops;
                    stops.resize(ivals.size());
                    std::transform(ivals.begin(), ivals.end(), stops.begin(), Interval<T, K>::intervalStop);
                    rightp = *std::max_element(stops.begin(), stops.end());
                }
                centerp = ivals.at(ivals.size() / 2).start;
                center = centerp;

                std::vector<Interval<T, K>> lefts, rights;
                for(auto i = ivals.begin(); i != ivals.end(); ++i){
                    if(i->stop < center){
                        lefts.push_back(*i);
                    }else if(i->start > center){
                        rights.push_back(*i);
                    }else{
                        intervals.push_back(*i);
                    }
                }
                if(!lefts.empty()){
                    left = new IntervalTree(lefts, depth, leftp, centerp, minbucket);
                }
                if(!rights.empty()){
                    right = new IntervalTree(rights, depth, centerp, rightp, minbucket);
                }
            }
        }

        /** find all intervals in the IntervalTree overlapping with an interval
         * @param start starting coordinate of interval
         * @param stop stopping coordinate of interval
         * @return all intervals in the IntervalTree overlapping with interval [start, stop]
         */
        std::vector<Interval<T, K>> findOverlapping(K start, K stop) const {
            std::vector<Interval<T, K>> ov;
            this->findOverlapping(start, stop, ov);
            return ov;
        }

        /** find all intervals in the IntervalTree overlapping with an interval
         * @param start starting coordinate of interval
         * @param stop stopping coordinate of interval
         * @param overlapping vector to store all intervals in the IntervalTree overlapping with interval [start, stop]
         */
        void findOverlapping(K start, K stop, std::vector<Interval<T, K>>& overlapping) const {
            if(!intervals.empty() && !(stop < intervals.front().start)){
                for(auto i = intervals.begin(); i != intervals.end(); ++i){
                    if(i->stop >= start && i->start <= stop){
                        overlapping.push_back(*i);
                    }
                }
            }
            if(left && start < center){
                left->findOverlapping(start, stop, overlapping);
            }
            if(right && stop > center){
                right->findOverlapping(start, stop, overlapping);
            }
        }


        /** find all intervals in the IntervalTree contained by an interval
         * @param start starting coordinate of interval
         * @param stop stopping coordinate of interval
         * @return all intervals in the IntervalTree contained by interval [start, stop]
         */
        std::vector<Interval<T, K>> findContained(K start, K stop) const {
            std::vector<Interval<T, K>> contained;
            this->findContained(start, stop, contained);
            return contained;
        }

        /** find all intervals in the IntervalTree contained by an interval
         * @param start starting coordinate of interval
         * @param stop stopping coordinate of interval
         * @param contained vector to store intervals in the IntervalTree contained by interval [start, stop]
         */
        void findContained(K start, K stop, std::vector<Interval<T, K>>& contained) const {
            if(!intervals.empty() && !(stop < intervals.front().start)){
                for(auto i = intervals.begin(); i != intervals.end(); ++i){
                    if(i->start >= start && i->stop <= stop){
                        contained.push_back(*i);
                    }
                }
            }
            if(left && start < center){
                left->findContained(start, stop, contained);
            }
            if(right && stop > center){
                right->findContained(start, stop, contained);
            }
        }
};

#endif
