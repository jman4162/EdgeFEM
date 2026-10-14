// Test-only check macros that stay active in Release builds.
//
// assert() is compiled out under -DNDEBUG, which CMake's Release
// configuration sets, so assert-based tests passed vacuously in CI.
// CHECK and CHECK_NEAR always evaluate, print file:line with the failing
// expression (and the values for CHECK_NEAR), and record the failure.
// A failed check does not stop the test.
//
// Each test's main() ends with `return edgefem_test::exit_code();`.
// As a backstop, an atexit handler turns any exit with recorded failures
// into exit status 1, so an early `return 0` cannot hide a failed check.
#pragma once

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace edgefem_test {

inline int &failure_count() {
  static int count = 0;
  return count;
}

inline void exit_guard() {
  if (failure_count() > 0) {
    std::fflush(stdout);
    std::fflush(stderr);
    std::_Exit(1);
  }
}

struct ExitGuardRegistration {
  ExitGuardRegistration() { std::atexit(exit_guard); }
};
inline ExitGuardRegistration exit_guard_registration;

inline void record_failure(const char *file, int line,
                           const std::string &message) {
  ++failure_count();
  std::cerr << file << ":" << line << ": CHECK FAILED: " << message
            << std::endl;
}

inline bool check(bool ok, const char *expr, const char *file, int line) {
  if (!ok)
    record_failure(file, line, expr);
  return ok;
}

template <typename A, typename B, typename Tol>
bool check_near(const A &a, const B &b, const Tol &tol, const char *a_expr,
                const char *b_expr, const char *tol_expr, const char *file,
                int line) {
  const double diff = std::abs(a - b);
  const bool ok = diff <= static_cast<double>(tol); // false for NaN
  if (!ok) {
    std::ostringstream msg;
    msg.precision(10);
    msg << "|" << a_expr << " - " << b_expr << "| <= " << tol_expr << " ("
        << a_expr << " = " << a << ", " << b_expr << " = " << b
        << ", |diff| = " << diff << ", tol = " << tol << ")";
    record_failure(file, line, msg.str());
  }
  return ok;
}

/// Exit status for main(): 0 when no check failed, 1 otherwise.
inline int exit_code() {
  if (failure_count() > 0) {
    std::cerr << failure_count() << " check(s) failed" << std::endl;
    return 1;
  }
  return 0;
}

} // namespace edgefem_test

#define CHECK(cond)                                                            \
  ::edgefem_test::check(static_cast<bool>(cond), #cond, __FILE__, __LINE__)

#define CHECK_NEAR(a, b, tol)                                                  \
  ::edgefem_test::check_near((a), (b), (tol), #a, #b, #tol, __FILE__, __LINE__)
