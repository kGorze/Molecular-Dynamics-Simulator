//BASED ON THE STACKOVERFLOW ANSWER

#include "Headers/utils/progressbar.h"



Progressbar::Progressbar() {
    _total_ticks = 0;
    _ticks_occurred = 0;
    _begin = Clock::now();

}

Progressbar::Progressbar(uint64_t ticks) {
    _total_ticks = ticks;
    _ticks_occurred = 0;
    _begin = Clock::now();
}

void Progressbar::tick() {
    ++_ticks_occurred;

    // Calculate progress percentage
    float percent_done = (float)_ticks_occurred / _total_ticks;

    // Calculate time remaining
    duration time_taken = std::chrono::steady_clock::now() - _begin;
    duration time_left = time_taken * (1 / percent_done - 1);

    std::chrono::minutes minutes_left = std::chrono::duration_cast<std::chrono::minutes>(time_left);
    std::chrono::seconds seconds_left = std::chrono::duration_cast<std::chrono::seconds>(time_left - minutes_left);

    // Output progress information
    std::cout << "\rProgress: [";
    int width = 50;  // Width of the progress bar
    int progress = width * percent_done;
    for (int i = 0; i < width; ++i) {
        if (i < progress) {
            std::cout << "=";
        }
        else {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(percent_done * 100.0) << " %"
        << " | Time left: " << minutes_left.count() << "m " << seconds_left.count() << "s";
    std::cout.flush();

}

void Progressbar::totalTicks(uint64_t ticks) {
    this->_total_ticks = ticks;
}
