#!/usr/bin/env python3
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

## @file
# Generate regular graph topologies in TGF format.
#
# @author Nuno Fachada
# @date 2019
# @copyright [Mozilla Public License Version 2.0](https://www.mozilla.org/en-US/MPL/2.0/)

import sys

# Check if number of arguments is correct
if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + " n k1 [k2 ...]")
    print("\tn - Population size")
    print("\tkX - Number of neighbors")
    exit(1)

# Get population size and initialize list with number of neighbors
n = int(sys.argv[1])
ks = []

# Parse number of neighbors parameters and put them in list
for k in sys.argv[2:]:
    ks.append(int(k))
    # Number of neighbors must be an odd value
    if ks[-1] % 2 == 0:
        print("All values of k must be odd")
        exit(1)

# Cycle through number of neighbors list
for k in ks:

    # Determine file name for current number of neighbors and open it
    # for writing
    filename = "n" + str(n) + "k" + str(k) + ".tgf"
    f = open(filename, "w")

    # Write nodes and their labels
    for node in range(1, n + 1):
        f.write(str(node) + " node_" + str(node) + "\n")

    # Write separator between nodes and edges/connections
    f.write("#\n")

    # Write edges/connections between nodes
    for node in range(1, n + 1):

        # For current node, write its edges/connections
        for v in range(-k // 2 + 1, k // 2 + 1):

            # Currently no need to specify connection with self
            if v == 0:
                continue

            # Get neighboring node
            neighNode = node + v
            if neighNode <= 0:
                neighNode = neighNode + n
            elif neighNode >= n + 1:
                neighNode = neighNode - n

            # Write connection between node and its neighbor
            f.write(str(node) + " " + str(neighNode) + "\n")

    # Close file for current number of neighbors
    f.close()
