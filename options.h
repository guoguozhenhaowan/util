#ifndef OPTIONS_H
#define OPTIONS_H

/** filter options struct */
struct FilterOpt{
    bool filterShortRead;    ///< filter too short read if true
    bool filterLongRead;     ///< filter too long read if true
    bool trimAdapter;        ///< trim adapter sequnce if true
    bool baseCorrection;     ///< correct base if true
    bool complexityFilter;   ///< filter low complexity tread if true
    int minReadLen;          ///< read length lower than this will be treated as too short
    int maxReadLen;          ///< read length greater than this will be treated as too long
    FilterOpt() = default;
    ~FilterOpt() = default;
};

#endif
