//
// Created by konrad_guest on 20/05/2024.
//

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <thread>
#include <chrono>
#include <iostream>
class Progressbar {
private:
    typedef std::chrono::steady_clock Clock;
    typedef Clock::time_point time_point;
    typedef Clock::period period;
    typedef std::chrono::duration<float, period> duration;
    std::uint64_t _total_ticks;
    std::uint64_t _ticks_occurred;
    time_point _begin;
public:
    Progressbar();
    Progressbar(uint64_t ticks);

    void tick();
    void totalTicks(uint64_t ticks);
};
#endif //PROGRESSBAR_H
