#ifndef UTESTUTILS_HPP
#define UTESTUTILS_HPP

#include <unordered_map>
#include <mutex>
#include <string>
#include <sstream>

class UTestRaceChecker {
    std::unordered_map<const void*,int> counterAccess;
    std::mutex counterAccessMutex;

    static std::string PtrToString(const void* ptr){
        std::stringstream ss;
        ss << ptr;
        return ss.str();
    }

public:
    UTestRaceChecker() = default;

    ~UTestRaceChecker() noexcept(false) {
        for(auto iter : counterAccess){
            if(iter.second != 0){
                throw std::runtime_error("Error in ~UTestRaceChecker for address : " + PtrToString(iter.first));
            }
        }
    }

    void lock(){
        counterAccessMutex.lock();
    }

    void unlock(){
        counterAccessMutex.unlock();
    }

    UTestRaceChecker(const UTestRaceChecker&) = delete;
    UTestRaceChecker(UTestRaceChecker&&) = delete;
    UTestRaceChecker& operator=(const UTestRaceChecker&) = delete;
    UTestRaceChecker& operator=(UTestRaceChecker&&) = delete;

    void addRead(const void* ptr){
        if(counterAccess.find(ptr) != counterAccess.end()
                && counterAccess[ptr] < 0){
            throw std::runtime_error("Error in addRead for address : " + PtrToString(ptr));
        }
        counterAccess[ptr] += 1;
    }

    void addWrite(const void* ptr){
        if(counterAccess.find(ptr) != counterAccess.end()
               && counterAccess[ptr] != 0){
            throw std::runtime_error("Error in addWrite for address : " + PtrToString(ptr));
        }
        counterAccess[ptr] = -1;
    }

    void releaseRead(const void* ptr){
        if(counterAccess.find(ptr) == counterAccess.end()
               || counterAccess[ptr] <= 0){
            throw std::runtime_error("Error in releaseRead for address : " + PtrToString(ptr));
        }
        counterAccess[ptr] -= 1;
    }

    void releaseWrite(const void* ptr){
        if(counterAccess.find(ptr) == counterAccess.end()
               || counterAccess[ptr] != -1){
            throw std::runtime_error("Error in releaseWrite for address : " + PtrToString(ptr));
        }
        counterAccess.erase(ptr);
    }
};


#endif
