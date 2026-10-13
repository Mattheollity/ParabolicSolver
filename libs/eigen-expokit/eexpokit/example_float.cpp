// eigen-expokit: example, float
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#include <complex>
#include <limits>

#include "example.hpp"

int main()
{
  std::cout << " *******************************************************************\n";
  std::cout << " float tol=" << 1.L/(1ull << (std::numeric_limits<float>::digits / 2)) << "\n";
  demo<float>();
}
