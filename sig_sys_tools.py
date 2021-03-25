"""Some Signals & Systems Helping Routines."""
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Circle
import numpy as np
from scipy import signal


def plot_zplane(z, p, k):
    """Plot pole/zero/gain plot of discrete-time, linear-time-invariant system.

    Note that the for-loop handling might be not very efficient
    for very long FIRs

    z...array of zeros in z-plane
    p...array of poles in z-zplane
    k...gain factor

    """
    # draw unit circle
    Nf = 2**7
    Om = np.arange(Nf) * 2*np.pi/Nf
    plt.plot(np.cos(Om), np.sin(Om), 'C7')

    try:  # TBD: check if this pole is compensated by a zero
        circle = Circle((0, 0), radius=np.max(np.abs(p)),
                        color='C7', alpha=0.15)
        plt.gcf().gca().add_artist(circle)
    except ValueError:
        print('no pole at all, ROC is whole z-plane')

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
        if pci > 1:  # if multiple poles exist then indicate the count
            plt.text(np.real(pui), np.imag(pui), pci)

    plt.text(0, +1, 'k=%f' % k)
    plt.text(0, -1, 'ROC for causal: white')
    plt.axis('square')
    # plt.axis([-2, 2, -2, 2])
    plt.xlabel(r'$\Re\{z\}$')
    plt.ylabel(r'$\Im\{z\}$')
    plt.grid(True)


def interpolate_dft2dtft(X, W):
    """DFT to DTFT interpolation.

    This is the reconstruction filter in frequency domain to get
    a finite length time sequence starting from k=0 out of a N-periodic
    sequence

    X...array containing DFT spectrum
    W...array with normalized digital frequencies
    typically W = np.arange(Nint) * 2*np.pi/Nint with desired Nint

    see e.g.
    Rabiner, Gold, 1975, Theory and Application of Digital Signal Processing
    Prentice Hall, eq. (2.142)

    """
    N = np.size(X)  # we estimate the DFT length from the DFT spectrum
    tmp_2piN = 2*np.pi/N
    tmp_N2 = (N-1)/2
    Xint = np.zeros_like(W, dtype='complex')
    for cW, vW in enumerate(W):  # counter, value
        for cX, vX in enumerate(X):
            W_tmp = vW - tmp_2piN * cX
            Xint[cW] += vX * diric(W_tmp, N) * np.exp(-1j*W_tmp*tmp_N2)
    return Xint


def plot_dtlti_analysis(z, p, k, fs=1, Nf=2**10, Nt=2**5):
    """Plot linear, time-invariant, discrete-time system.

    Impulse, step, frequency response (level/phase/group delay) and z-plane
    of transfer function H(z) given as zeros/poles/gain description

    Note that we use fs only for the frequency responses

    """
    # still hard coded, TBD:
    # figure size
    # group_delay discontinuities and automatic ylim, yticks is not
    # useful sometimes

    plt.figure(figsize=(9, 9), tight_layout=True)

    mu = np.arange(Nf)
    df = fs/Nf
    dW = 2*np.pi/Nf
    f = mu*df  # frequency vector [0...fs)
    W = mu*dW  # digital angular frequency [0...2pi)

    sys = signal.dlti(z, p, k, dt=True)
    sys_ba = signal.TransferFunction(sys)
    b = sys_ba.num   # we need coeff b,a for group_delay()
    a = sys_ba.den

    [W, H] = signal.dlti.freqresp(sys, W)
    [W, gd] = signal.group_delay((b, a), w=W)  # gd in samples
    h = signal.dimpulse(sys, n=Nt)
    he = signal.dstep(sys, n=Nt)

    # plot frequency response: level
    ax1 = plt.subplot(3, 2, 1)
    ax1t = ax1.twiny()

    ax1.grid(True, color='lavender', which='major')
    # plotted below grid :-(
    # ax1t.grid(True, color='mistyrose', which='minor')

    ax1.plot(W/(2*np.pi), 20*np.log10(np.abs(H)), color='C0', lw=2)
    ax1.set_xscale('linear')
    ax1.set_xlim([0, 1])
    ax1.tick_params(axis='x', labelcolor='C0')
    ax1.set_xlabel(r'$\frac{\Omega}{2\pi} = \frac{f}{f_s}$', color='C0')
    ax1.set_xticks(np.arange(11)/10)
    ax1.set_ylabel('level in dB   or   20 lg|H| / dB', color='k')
    ax1.set_axisbelow(True)

    ax1t.plot(f, 20*np.log10(np.abs(H)), color='C3', lw=2)
    ax1t.set_xscale('log')
    ax1t.set_xlim([f[1], f[-1]])
    ax1t.tick_params(axis='x', labelcolor='C3')
    ax1t.set_xlabel('frequency in Hz   or   f / Hz for $f_s$ = ' +
                    '{0:5.2f}'.format(fs)+' Hz', color='C3')
    ax1t.set_axisbelow(True)

    # plot impulse response
    ax2 = plt.subplot(3, 2, 2)
    ax2.stem(np.squeeze(h[0]), np.squeeze(h[1]),
             use_line_collection=True,
             linefmt='C0:', markerfmt='C0o', basefmt='C0:')
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.grid(True)
    ax2.set_xlabel(r'$k$')
    ax2.set_ylabel(r'impulse response $h[k]$')

    # plot frequency response: phase
    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(W, np.unwrap(np.angle(H)))
    ax3.set_xlabel(
        r'digital angular frequency in radian   or   $\Omega$ / rad')
    ax3.set_ylabel(r'phase in radian   or $\angle$ H / rad')
    ax3.set_xlim([0, 2*np.pi])
    ax3.set_xticks(np.arange(9)*np.pi/4)
    ax3.set_xticklabels([r'0', r'$\pi/4$', r'$\pi/2$', r'$3/4\pi$',
                         r'$\pi$', r'$5/4\pi$', r'$3/2\pi$', r'$7/4\pi$',
                         r'$2\pi$'])
    # ax3.set_ylim([-180, +180])
    # ax3.set_yticks(np.arange(-180, +180+45, 45))
    ax3.grid(True)

    # plot step response
    ax4 = plt.subplot(3, 2, 4)
    ax4.stem(np.squeeze(he[0]), np.squeeze(he[1]), use_line_collection=True,
             linefmt='C0:', markerfmt='C0o', basefmt='C0:')
    ax4.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax4.grid(True)
    ax4.set_xlabel(r'$k$')
    ax4.set_ylabel(r'step response $h_\epsilon[k]$')

    # plot frequency response: group delay
    ax5 = plt.subplot(3, 2, 5)
    ax5.plot(W/np.pi, gd/fs)
    ax5.set_xscale('linear')
    ax5.set_xlim([0, 2])
    ax5.set_xlabel(r'$\frac{\Omega}{\pi}$')
    ax5.set_ylabel(r'group delay in seconds   or   $\tau_\mathrm{GD}$ / s')
    ax5.grid(True, which='both')

    # zplane
    ax6 = plt.subplot(3, 2, 6)
    plot_zplane(sys.zeros, sys.poles, sys.gain)  # see function above

    # plt.tight_layout(True)
