import numpy as np
import matplotlib.pyplot as plt

def zplane_plot(z, p, k):  # this might be not very efficient for very long FIRs
    # draw unit circle
    Nf = 2**7
    Om = np.arange(Nf) * (2*np.pi/Nf)
    plt.plot(np.cos(Om), np.sin(Om), 'C7')

    zu, zc = np.unique(z, return_counts=True)  # find and count unique zeros
    for zui, zci in zip(zu, zc):  # plot them individually
        plt.plot(np.real(zui), np.imag(zui), ms=7,
                 color='C0', marker='o', fillstyle='none')
        if zci > 1:  # if multiple zeros exist then indicate the count
            plt.text(np.real(zui), np.imag(zui), zci)

    pu, pc = np.unique(p, return_counts=True)  # find and count unique poles
    for pui, pci in zip(pu, pc):  # plot them individually
        plt.plot(np.real(pui), np.imag(pui), ms=7,
                 color='C3', marker='x')
        if pci > 1:  # if multiple poles exist indicate the count
            plt.text(np.real(pui), np.imag(pui), pci)

    plt.text(0.5, 1.5, 'k=%f' % k)
    plt.axis('square')
    plt.axis([-2, 2, -2, 2])
    plt.xlabel(r'$\Re\{z\}$')
    plt.ylabel(r'$\Im\{z\}$')
    plt.grid(True)