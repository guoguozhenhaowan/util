#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <fstream>
#include <stack>
#include <map>
#include "util.h"

void parsePhy(const std::string& inphy, std::stack<std::vector<std::string>>& mergeStack){
    const float MAX_DISTANCE_TO_MERGE = 0.01;
    const float MIN_LEAVES_TO_MERGE = 5;
    std::stack<std::pair<std::vector<std::string>, float>> nodeStack;
    std::ifstream fr(inphy);
    std::string line;
    std::vector<std::string> strVec;
    std::pair<std::vector<std::string>, float> lastNode, curNode;
    while(std::getline(fr, line)){
        if(util::starts_with(line, "(") || util::starts_with(line, ",")){
            continue;
        }
        if(util::starts_with(line, ";")){
            break;
        }
        if(util::starts_with(line, ")")){
            if(nodeStack.size() <= 1){
                continue;
            }
            curNode = nodeStack.top();
            nodeStack.pop();
            lastNode = nodeStack.top();
            nodeStack.pop();
            util::split(line, strVec, ":");
            int mergeMarker = (lastNode.second <= MAX_DISTANCE_TO_MERGE) + (curNode.second <= MAX_DISTANCE_TO_MERGE) * 2;
            switch(mergeMarker){
                case 3:
                    std::copy(curNode.first.begin(), curNode.first.end(), std::back_inserter(lastNode.first));
                    if(!util::ends_with(line, ")")){
                        lastNode.second = std::atof(strVec[1].c_str());
                    }
                    nodeStack.push(lastNode);
                    break;
                case 2:
                    if(curNode.first.size() >= MIN_LEAVES_TO_MERGE){
                        mergeStack.push(curNode.first);
                    }
                    nodeStack.push(lastNode);
                    break;
                case 1:
                    if(lastNode.first.size() >= MIN_LEAVES_TO_MERGE){
                        mergeStack.push(lastNode.first);
                    }
                    nodeStack.push(curNode);
                    break;
                default:
                    if(curNode.first.size() >= MIN_LEAVES_TO_MERGE){
                        mergeStack.push(curNode.first);
                    }
                    if(lastNode.first.size() >= MIN_LEAVES_TO_MERGE){
                        mergeStack.push(lastNode.first);
                    }
                    break;
            }
            continue;
        }
        util::split(line, strVec, ":");
        std::vector<std::string> branch = {strVec[0]};
        nodeStack.push(std::make_pair(branch, std::atof(strVec[1].c_str())));
    }
    while(!nodeStack.empty()){
        curNode = nodeStack.top();
        if(curNode.first.size() >= MIN_LEAVES_TO_MERGE){
            mergeStack.push(curNode.first);
        }
        nodeStack.pop();
    }
    fr.close();
}

int main(int argc, char** argv){
    if(argc < 2){
        std::cout << argv[0] << " <inphy> " << std::endl;
        return 0;
    }
    char* inphy = argv[1];
    std::stack<std::vector<std::string>> mergeStack;
    parsePhy(inphy, mergeStack);
    std::vector<std::string> group;
    int i = 0;
    while(!mergeStack.empty()){
        group = mergeStack.top();
        std::cout << "group" << (++i) << ":";
        for(size_t m = 0; m < group.size(); ++m){
            std::cout << group[m] << ",";
        }
        std::cout << std::endl;
        mergeStack.pop();
    }
}
