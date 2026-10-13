// eigen-expokit: example, complex<float>
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
  std::cout << " complex float\n";
  demo<std::complex<float>>();
}
