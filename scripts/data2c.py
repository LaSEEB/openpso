#!/usr/bin/env python3
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

## @file
# Convert data matrices in txt format to C 2D arrays.
#
# @author Nuno Fachada
# @date 2017
# @copyright [Mozilla Public License Version 2.0](https://www.mozilla.org/en-US/MPL/2.0/)

import sys

if len(sys.argv) == 1:
	print("Usage: " + sys.argv[0] + " file.txt")
	exit(1)

f = open(sys.argv[1], "r")
lines = f.readlines()
f.close()

sys.stdout.write("{")

for i in range(len(lines)):

	values = lines[i].split()
	sys.stdout.write("{")

	for j in range(len(values)):

		sys.stdout.write(values[j])
		if j != len(values) - 1:
			sys.stdout.write(", ")

	sys.stdout.write("}")
	if i != len(lines) - 1:
		sys.stdout.write(",\\\n")

sys.stdout.write("}")
