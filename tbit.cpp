#include <iostream>
#include <unistd.h>
#include <cstdio>

void writeBits(const size_t v, int fd)
{
    if(!v){
        std::putchar ('0'); 
        return; 
    }

    size_t sz = sizeof(v) * CHAR_BIT;
    size_t rem = 0;
    
    while(sz--){
        if((rem = v >> sz)){
            write(fd, (rem & 1) ? "1" : "0", 1);
        }
    }
}

int main(int argc, char** argv){
    if(argc == 1){
        std::cout << std::string(argv[0]) << " <number> " << std::endl;
        exit(0);
    }
    std::string str(argv[1]);
    int val = 0;
    if(str.length() == 1){
        val = str[0] - val;
    }else{
        val = std::atoi(argv[1]);
    }
    writeBits(val, STDOUT_FILENO);
    std::cout << std::endl;
}
