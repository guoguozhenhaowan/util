#include <iostream>
#include <unistd.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>

void writeBits(int v){
    std::ostringstream ss;
    while(v > 0){
        ss << (v % 2);
        v /= 2;
    }
    std::string re = ss.str();
    std::reverse(re.begin(), re.end());
    std::cout << std::setw(8) << std::setfill('0') << re  << std::endl;
}

int main(int argc, char** argv){
    if(argc == 1){
        std::cout << std::string(argv[0]) << " <number> " << std::endl;
        exit(0);
    }
    int val = 0;
    if(std::strlen(argv[1]) == 1){
        val = argv[1][0] - '\0';
    }else{
        val = std::atoi(argv[1]);
    }
    writeBits(val);
}
