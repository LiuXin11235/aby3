#include <chrono>
#include <iostream>
#include <unordered_map>
#include <string>
#include "../aby3-RTR/debug.h"

class Timer{
public:

    // Singleton pattern.
    static Timer& getInstance() {
        static Timer instance; // Guaranteed to be destroyed.
                               // Instantiated on first use.
        return instance;
    }

    Timer(Timer const&) = delete; // you can not copy Timer instance.
    void operator=(Timer const&) = delete; // you can not assign Timer instance.

    void start(const std::string& key) {
        timestamps[key] = std::chrono::high_resolution_clock::now();
    }

    void end(const std::string& key) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto start_time = timestamps[key];
        durations[key] = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    }

    void print(const std::string& key, const std::string& unit = "microseconds") {
        // check if the key exists
        if(durations.find(key) == durations.end()){
            THROW_RUNTIME_ERROR("Key: " + key + " not found in the timer");
        }
        
        double duration = durations[key];
        if (unit == "milliseconds") {
            duration /= 1000;
        } else if (unit == "seconds") {
            duration /= 1000000;
        } else if (unit == "minutes") {
            duration /= 60000000;
        }
        std::cout << "Time taken by " << key << ": " << duration << " " << unit << std::endl;
    }

    void print(const std::string& unit = "microseconds"){
        for(auto& kv : durations){
            print(kv.first, unit);
        }
    }

private:
    Timer() {} // you can not construct Timer instance from outside.
    std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> timestamps;
    std::unordered_map<std::string, double> durations;
};