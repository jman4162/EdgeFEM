import numpy as np
import matplotlib.pyplot as plt

def plot_sparams(filepath="waveguide_sparams.s2p"):
    """
    Reads a 2-port Touchstone file and plots the magnitude of S11 and S21.
    """
    try:
        data = np.loadtxt(filepath, skiprows=1)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        print("Please run the waveguide_demo first to generate the S-parameter file.")
        return

    freq = data[:, 0]
    s11 = data[:, 1] + 1j * data[:, 2]
    s21 = data[:, 3] + 1j * data[:, 4]

    plt.figure(figsize=(10, 6))
    plt.plot(freq / 1e9, 20 * np.log10(np.abs(s11)), label="|S11| (dB)")
    plt.plot(freq / 1e9, 20 * np.log10(np.abs(s21)), label="|S21| (dB)")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Magnitude (dB)")
    plt.title("S-Parameters vs. Frequency")
    plt.grid(True)
    plt.legend()
    plt.ylim(-50, 5)

    plot_filepath = "waveguide_sparams.png"
    plt.savefig(plot_filepath)
    print(f"Plot saved to {plot_filepath}")

if __name__ == "__main__":
    plot_sparams()
