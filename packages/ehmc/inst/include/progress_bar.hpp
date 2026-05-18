#ifndef INST_INCLUDE_PROGRESS_BAR_HPP_
#define INST_INCLUDE_PROGRESS_BAR_HPP_

#include <Rcpp.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>  // For rounding in time display
#include <string> // For optional messages

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

std::string format_time(double seconds);

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

class progress_bar {
public:
  // --- Constructor
  progress_bar(size_t n_iter)
    :   bar_width(50), fill("■"), remainder(" "),
        n_iter(n_iter), current_iter(0), ended(false) {
    start_time = std::chrono::steady_clock::now();
  }
  // --- Destructor
  virtual ~progress_bar(){finish();};

  // --- Setters
  void set_bar_width(size_t width) {
    std::unique_lock<std::mutex> lock{mutex_};
    bar_width = width;
  }

  void set_current_cursor(const std::string& chars) {
    std::unique_lock<std::mutex> lock{mutex_};
    current = chars;
  }

  void fill_bar_progress_with(const std::string& chars) {
    std::unique_lock<std::mutex> lock{mutex_};
    fill = chars;
  }

  void fill_bar_remainder_with(const std::string& chars) {
    std::unique_lock<std::mutex> lock{mutex_};
    remainder = chars;
  }

  void updt_current_iter() {
    std::unique_lock<std::mutex> lock{mutex_};  // CTAD (C++17)
    ++current_iter;
  }

  // --- Display methods
  void show_progress(const std::string& message = "",
                     double value = 0.0) {
    this->updt_current_iter();
    this->update(message, value);
  }

  void update(const std::string& message = "",
              double value = 0.0) {
    // --- Update the message printed in the console

    // --- Check for user interrupt by catching the exception
    try {
      Rcpp::checkUserInterrupt();  // --- Check if R session has been interrupted
    } catch (const Rcpp::exception&) {
      this->finish();  // --- Clean up and exit if interrupted
      return;
    }

    // --- Check if the progress bar has already ended
    {
      std::unique_lock<std::mutex> lock{mutex_};
      if (ended) {
        return; // Do not update if already finished
      }
    }

    // --- Calculate percentage completion
    double progress = (n_iter == 0) ? 1.0 : static_cast<double>(current_iter) / n_iter;
    unsigned percentage = static_cast<unsigned>(progress * 100.0);

    // --- Get the current time and calculate elapsed time
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = now - start_time;

    // --- Estimate total time and remaining time
    double estimated_total_time = elapsed.count() / progress;
    double remaining_time = estimated_total_time - elapsed.count();

    // --- Display the progress bar
    this->display(percentage, remaining_time, message, value);
  }

  void display(unsigned percentage,
               double remaining_time,
               const std::string& message,
               double value) {
    // --- Display the progress bar and the estimated remaining time
    std::unique_lock<std::mutex> lock{mutex_};

    std::cout<< " [";
    size_t completed = static_cast<size_t>(static_cast<double>(bar_width) * static_cast<double>(percentage) / 100.0);

    // --- Fill the bar based on progress
    for (size_t i = 0; i < bar_width; i++) {
      if (i < completed) std::cout << fill;
      else if (i == completed) std::cout << current;
      else std::cout << remainder;
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << percentage << "%";

    // --- Display estimated remaining time (in seconds)
    std::cout << " | ETA: " << std::fixed << format_time(remaining_time) ;

    // --- Print the optional message if provided
    if (!message.empty()) {
      std::cout << std::fixed << " (" << message << value << ")";
    }

    std::cout << "\r";  // Move the cursor to the beginning of the line
    std::cout.flush();
  }

  // --- End the progress bar and print a newline
  void finish() {
    std::unique_lock<std::mutex> lock{mutex_};
    if (!ended) {
      std::cout << std::endl;  // Add a new line at the end
      ended = true; // Mark as ended to prevent further calls
    }
  }

private:
  std::mutex mutex_;
  size_t bar_width;
  std::string fill;
  std::string remainder;
  std::string current;
  size_t n_iter;
  size_t current_iter;
  bool ended;
  std::chrono::steady_clock::time_point start_time;
};

#endif /* INST_INCLUDE_PROGRESS_BAR_HPP_ */
