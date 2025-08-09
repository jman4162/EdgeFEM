#!/usr/bin/env python3
import csv
import math
import sys
import matplotlib.pyplot as plt


def main(path):
    theta = []
    gain = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            theta.append(math.radians(float(row["theta_deg"])))
            # simple magnitude -> dB
            eth = float(row.get("eth_mag", 0.0))
            eph = float(row.get("eph_mag", 0.0))
            mag = math.sqrt(eth ** 2 + eph ** 2)
            gain.append(20 * math.log10(mag) if mag > 0 else -120)
    ax = plt.subplot(111, projection="polar")
    ax.plot(theta, gain)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title("Far-field pattern")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: plot_pattern.py pattern.csv")
        sys.exit(1)
    main(sys.argv[1])
