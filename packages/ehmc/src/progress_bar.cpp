#include <progress_bar.hpp>

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

std::string format_time(double seconds) {
  // --- Convert seconds to hours, minutes, seconds format
  double temp = std::max(seconds, 0.0);  // Avoid negative time
  unsigned hours = static_cast<unsigned>(seconds) / 3600;
  temp -= hours * 3600;
  unsigned minutes = static_cast<unsigned>(seconds) / 60;
  temp -= minutes * 60;

  std::ostringstream oss;
  unsigned d = (seconds < 60) ? 1 : 0;
  if (seconds > 3600) oss << hours << "h";
  if (seconds > 60) oss << minutes << "m";
  oss << std::fixed << std::setprecision(d) << temp << "s";
  return oss.str();
}
