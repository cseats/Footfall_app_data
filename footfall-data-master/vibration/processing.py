import numpy
import scipy.signal as signal
import scipy.fftpack as fft
import matplotlib.pyplot as plt
import math
import json
import pandas as pd


def processing(data, fs=100, w='b', do_plots=0):
    #  function [vdv,af,dt] = bs6841vdv(acc,fs,w,do_plots)
    #
    # BS6841 filtering and evaluation of vibration dosage values (vdv)
    # filter acc using the band limiting and frequency weighting filters
    # specified in BS6841.
    #
    # fs is the sampling frequency
    #
    # w is a single character string, 'b' through to 'g'
    # specifying weighting to be used.
    #
    # The following frequency weightings are available:
    #
    # Wb, Wc, Wd, We, Wf, Wg
    #                                                   T
    #                                                (  (   4     )1/4
    # vdv is the vibration dosage value evaluated as (  |  a  dt  )
    #                                                (  )         )
    #                                                   0
    #
    # af is the filtered acceleration, which may have a different time base
    # if decimation is used
    #
    # Sampling frequency is checked to ensure that it is adequate
    # for the band chosen. An error is issued if the Nyquist frequency
    # is less than the ~3dB point of the combined filter, and a warning is
    # issued if it is less than approx 20dB down point.
    #
    # Signal length is checked to ensure that there is at least 2 whole
    # period at the lower 20dB down point. For all except the 'b'
    # weighting this also implies that there are 6 whole periods at the 3dB
    # down point. For the 'b' weighting, because of the step in the frequency
    # response, this implies 28 periods at the lower 3dB point.
    #
    # Over-sampled signals are decimated to ensure filter stability, but not
    # too drastically to try and ensure that pre-warping is not necessary
    # for the bilinear transform from the s domain to the z domain. The
    # responses of the analogue filters (as specified in the standard) and
    # their digital realisations (as implimented by Matlab) should be compared.
    # If there appear to be significant discrepancies then either a) try
    # a higher sampling frequency, b) consider modifying these routines to
    # include pre-warp in the bilinear transform
    #
    # do_plots = 1 to get filtered and non-filtered traces
    # do_plots = 2 to get traces and filter characteristics
    # do_plots = 0 for no_plots
    # do_plots = -1 to plot only filter response for given weighting and
    # sampling freq (a is ignored) into CURRENT figure window
    # do_plots = -2 to do all analysis with minimal output to screen
    #
    # PY, 28/05/02
    # Nicholas Simpson 10/11/14 - dt added to output
    #
    # See also BS6841BANDLIMIT, BS6841FREqWEIGHT

    # check input arguments
    # check sampling frequency
    #  set fminsoft, fminhard and t_min
    #  frequency response of combined band limit and frequency weighting
    #  filters are approx 3dB down at fminhard and 20dB down at fminsoft
    #  If sampling frequency is less than twice fminhard, am error occurs,
    #  if it is less than fminsoft, a warning is issued that sampling
    #  frequency may be insufficient to represent all frequencies in that
    #  particular band.
    #  Signal length is checked to ensure that there are at least 2 periods
    #  of t_min

    k = 0.0
    f = []
    q = []

    if w == 'b':
        fminsoft = 100
        fminhard = 10
        t_min = 1 / 0.2
        f.append(0.4)
        f.append(100)
        f.append(16)
        f.append(16)
        f.append(2.5)
        f.append(4)
        q.append(0.71)
        q.append(0.55)
        q.append(0.9)
        q.append(0.95)
        k = 0.4
    elif w == 'c':
        fminsoft = 50
        fminhard = 7
        t_min = 1 / 0.1
        f.append(0.4)
        f.append(100)
        f.append(8)
        f.append(8)
        q.append(0.71)
        q.append(0.63)
        k = 1.0
    elif w == 'd':
        fminsoft = 11
        fminhard = 2
        t_min = 1 / 0.1
        f.append(0.4)
        f.append(100)
        f.append(2)
        f.append(2)
        q.append(0.71)
        q.append(0.63)
        k = 1.0
    elif w == 'e':
        fminsoft = 10
        fminhard = 1
        t_min = 1 / 0.1
        f.append(0.4)
        f.append(100)
        f.append(1)
        f.append(1)
        q.append(0.71)
        q.append(0.63)
        k = 1.0
    elif w == 'f':
        fminsoft = 0.8
        fminhard = 0.2
        t_min = 1 / 0.04
        f.append(0.08)
        f.append(0.63)
        f.append(1e6)
        f.append(0.25)
        f.append(0.0625)
        f.append(0.1)
        q.append(0.71)
        q.append(0.86)
        q.append(0.8)
        q.append(0.8)
        k = 0.4
    elif w == 'g':
        fminsoft = 50
        fminhard = 7
        t_min = 1 / 0.4
        f.append(0.8)
        f.append(100)
        f.append(1.5)
        f.append(5.3)
        q.append(0.71)
        q.append(0.68)
        k = 0.42
    else:
        print('Invalid weighting')

    if fs < 2 * fminhard:
        print('Increase sampling frequency to at least %3.2f Hz', fs, 2 * fminhard)
        # error('Sampling frequency (%3.2f Hz) insufficient to represent vibration for this weighting')
    elif fs < 2 * fminsoft and do_plots != -2:
        print('Warning, sampling frequency (%3.2f Hz) insufficientto represent higher frequencies for this weighting',
              fs)
        print('Also performance of digital filter at higher frequencies may be compromised')
        print('Check by comparing analogue and digital filter characteristics')
        print('Type help bs6841vdv for more info')
        print('Increase sampling frequency to at least %3.2f Hz wil help', 2 * fminsoft)
    results = []
    for recording in data:
        # check length of signal
        acc = []
        for point in recording:
            acc.append(point['acceleration_z'])
        acc = numpy.array(acc)
        np = len(acc)
        dt = 1 / fs
        t = dt * np
        if t < t_min:
            print('A signal length of at least %2.1 secs is required', t_min)
            print('Insufficient signal length to give representative values of rms/rmq etc.')

        # decimate if fs > 20*fminsoft
        fp, fm = oafft(acc, 1 / fs, 1, 0)  # get original spectral content before decimation
        fsold = fs
        while fs > 20 * fminsoft:
            accdec = []
            for i in range(1, len(acc)):
                accdec[:, i] = signal.resample(acc[:, i], 1, 2)
            acc = accdec
            # clear accdec
            fs = fs / 2

        if fsold != fs:
            print('N.B. Decimation used to reduce sampling frequency from %f to %f ', fsold, fs)

        # get transfer function (analogue) frequency weighting
        bfw, afw = bs6841freqweight(w, k, f, q)

        # get transfer function (analogue) band limiting (descing powers of s)
        bbl, abl = bs6841bandlimit(f, q)

        # calculate analogue frequency responses (descing powers of s)
        f = numpy.logspace(-2, 3, 1024)
        hfw = signal.freqs(bfw, afw, f)
        hbl = signal.freqs(bbl, abl, f)
        # combine
        a = numpy.zeros(9)
        b = numpy.zeros(9)

        for jj in range(0, 5):
            a[jj: jj + 4] = a[jj: jj + 4] + numpy.multiply(afw[jj], abl[0:4])
            b[jj: jj + 4] = b[jj: jj + 4] + numpy.multiply(bfw[jj], bbl[0:4])  # convert to z transform

        bz, az = signal.bilinear(b, a, fs)

        # filter characteristics
        if do_plots == 1 or do_plots == -1:
            hz = signal.freqz(bz, az, f, fs)
            hs = signal.freqs(b, a, f)
            if do_plots == -1:
                plt.subplot(2, 1, 1)
                plt.loglog(f / (2 * math.pi), abs(hfw), f / (2 * math.pi), abs(hbl))
                plt.title(['Analogue Frequency Responses: ', w, ' Weighting'])
                plt.legend('Frequency Weighting', 'Band Limiting')
                plt.axis([10 ^ -2, 10 ^ 2, 10 ^ -2, 10])
                plt.ylabel('Filter Modulus')
                plt.subplot(2, 1, 2)
                plt.loglog(f / (2 * math.pi), abs(hs), f, abs(hz))
                plt.legend('Analogue Filter', 'Digital Filter')
                plt.title('Combined Frequency Responses')
                plt.axis([10 ^ -2, 10 ^ 2, 10 ^ -2, 10])
                plt.ylabel('Filter Modulus')
                af = []
                vdv = []
                print(
                    'doplot = -1, processing of acceleration is ignored, returning empty ''vdv'' and ''af'' variables.')
                print('Type help bs6841vdv for more info.')
                return

        # filter
        af = signal.lfilter(bz, az, acc)
        peak = max(af)
        rms = numpy.power((numpy.mean(numpy.power(af, 2))), .5)
        cf = peak / rms

        # calculate vdv
        dt = 1 / fs
        duration = dt * len(af)
        vdv = (sum(dt * (af ** 4))) ** 0.25

        # calculate frequency response of digital and analogue filters:

        # plot weighting and band limiting responses, digital and analogue
        # characteristics
        if do_plots == 2:
            hz = signal.freqz(bz, az, f, fs)
            hs = signal.freqs(b, a, f)
            fpf, fmf = oafft(af, 1 / fs, 1, 0)
            # plt.figure()
            plt.subplot(3, 1, 1)
            plt.loglog(f / (2 * math.pi), abs(hfw), f / (2 * math.pi), abs(hbl))
            plt.grid()
            plt.title(['Analogue Frequency Responses: ', w, ' Weighting'])
            plt.legend('Frequency Weighting', 'Band Limiting')
            plt.axis([10 ^ -2, 10 ^ 2, 10 ^ -2, 10])
            plt.ylabel('Filter Modulus')
            plt.subplot(3, 1, 2)
            plt.loglog(f / (2 * math.pi), abs(hs), f, abs(hz))
            plt.legend('Analogue Filter', 'Digital Filter')
            plt.title('Combined Frequency Responses')
            plt.axis([10 ^ -2, 10 ^ 2, 10 ^ -2, 10])
            plt.ylabel('Filter Modulus')
            plt.subplot(3, 1, 3)
            plt.loglog(fp, fm, fpf, fmf)
            plt.legend('Original Signal', 'Filtered Signal')
            plt.xlabel('Frequency, Hz')
            plt.ylabel('FFT Modulus')
            # %axis([10^-2, 10^2,10^-2,10])

        if do_plots == 1 or do_plots == 2:
            # plt.figure()
            plt.subplot(2, 1, 1)
            tt = numpy.multiply(dt, range(0, len(af)))
            orig = plt.plot(tt, acc, label='Original')
            filt = plt.plot(tt, af, label='Filtered')
            # plt.legend([orig, filt],['Original Signal', 'Filtered Signal'])
            plt.ylabel('Acceleration, m/s^2')
            plt.title('Filtered Signal for VDV calculation')
            plt.subplot(2, 1, 2)
            plt.plot(tt, af)
            # calculate crest factor
            peak = max(af)
            rms = numpy.power((numpy.mean(numpy.power(af, 2))), .5)
            cf = peak / rms
            plt.title('Filtered Signal, Crest Factor = %g ' % cf)
            plt.ylabel('Acceleration, m/s^2')
            plt.xlabel('Time, secs')

        # plt.show()

        results.append({'vdv': vdv, 'peak_acc': peak, 'crest_factor': cf, 'duration': duration})
        # results.append({'vdv': vdv, 'af': af, 'dt': dt, 'peak_acc': peak, 'crest_factor': cf})
    # print(results)

    pd.Series(results).to_json('data\\output.json')
    with open('data\\output.json', 'w') as outfile:
        json.dump(results, outfile)
    vdvs = []
    cfs = []
    peaks = []
    durations = []
    for entry in results:
        vdvs.append(entry['vdv'])
        cfs.append(entry['crest_factor'])
        peaks.append(entry['peak_acc'])
        durations.append(entry['duration'])

    # plt.figure(1)
    # plt.hist(durations, bins=100)
    # plt.xlabel('Duration of Recording')

    # plt.subplot(3, 1, 1)
    # plt.scatter(durations, vdvs)
    #
    # plt.subplot(3, 1, 2)
    # plt.scatter(cfs, vdvs)
    #
    # plt.subplot(3, 1, 3)
    # plt.scatter(cfs, peaks)
    plt.show()


def bs6841freqweight(w: str, k: float, f: list, q: list):
    a = [0, 0, 0, 0, 0]
    b = [0, 0, 0, 0, 0]
    a_ = [0, 0, 0, 0, 0]
    b_ = [0, 0, 0, 0, 0]
    c = [0, 0, 0]
    d = [0, 0, 0]

    # for i in range(0, 5):
    #     a.append(0)
    #     b.append(0)

    if w == "b" or w == "f":
        b_[0] = 8 * (math.pi ** 2) * f[2] * (f[4] ** 2)
        b_[1] = 4 * math.pi * math.pi * (f[4] ** 2) + 4 * (math.pi ** 2) * f[2] * f[4] / q[2]
        b_[2] = 2 * math.pi * f[2] + 2 * math.pi * f[4] / q[2]
        b_[4] = 1

        for i in range(0, len(b_)):
            b_[i] = b_[i] * 2 * math.pi * k * (f[4] ** 2) * (f[5] ** 2)

        c[0] = 4 * (math.pi ** 2) * (f[3] ** 2)
        c[1] = 2 * math.pi * f[3] / q[1]
        c[2] = 1
        d[0] = 4 * (math.pi ** 2) * (f[5] ** 2)
        d[1] = 2 * math.pi * f[5] / q[3]
        d[2] = 1

        a_[0] = a_[0] + c[0] * d[0]
        a_[1] = a_[1] + c[0] * d[1]
        a_[2] = a_[2] + c[0] * d[2]
        a_[1] = a_[1] + c[1] * d[0]
        a_[2] = a_[2] + c[1] * d[1]
        a_[3] = a_[3] + c[1] * d[2]
        a_[2] = a_[2] + c[2] * d[0]
        a_[3] = a_[3] + c[2] * d[1]
        a_[4] = a_[4] + c[2] * d[2]

        for i in range(0, len(a_)):
            a_[i] = a_[i] * f[2] * (f[4] ** 2)
    else:
        b_[1] = 2 * math.pi * k * (f[3] ** 2)
        b_[0] = b_[1] * 2 * math.pi * f[2]
        a_[0] = f[2]
        a_[1] = f[2]
        a_[2] = f[2]
        a_[0] = a_[0] * 4 * math.pi * math.pi * (f[3] ** 2)
        a_[1] = a_[1] * 2 * math.pi * f[3] / q[1]

    # now we need to flip around a and b coefficients
    for i in range(0, len(a_)):
        a[i] = a_[len(a_) - i - 1]

    for i in range(0, len(b_)):
        b[i] = b_[len(b_) - i - 1]

    return a, b


def bs6841bandlimit(f: list, q: list):
    a = [0, 0, 0, 0, 0]
    b = [0, 0, 0, 0, 0]
    a_ = [0, 0, 0, 0, 0]
    b_ = [0, 0, 0, 0, 0]
    c = [0, 0, 0]
    d = [0, 0, 0]

    # for i in range(0, 4):
    #     a.append(0)
    #     b.append(0)

    b_[2] = 4 * (math.pi ** 2) * (f[1] ** 2)
    c[2] = 1
    c[1] = 2 * math.pi * f[0] / q[0]
    c[0] = 4 * (math.pi ** 2) * (f[0] ** 2)
    d[2] = 1
    d[1] = 2 * math.pi * f[1] / q[0]
    d[0] = 4 * (math.pi ** 2) * (f[1] ** 2)

    a_[0] = a_[0] + c[0] * d[0]
    a_[1] = a_[1] + c[0] * d[1]
    a_[2] = a_[2] + c[0] * d[2]
    a_[1] = a_[1] + c[1] * d[0]
    a_[2] = a_[2] + c[1] * d[1]
    a_[3] = a_[3] + c[1] * d[2]
    a_[2] = a_[2] + c[2] * d[0]
    a_[3] = a_[3] + c[2] * d[1]
    a_[4] = a_[4] + c[2] * d[2]

    for i in range(0, len(a_)):
        a[i] = a_[len(a_) - i - 1]
    for i in range(0, len(b_)):
        b[i] = b_[len(b_) - i - 1]

    return a, b


def oafft(sig, dt: float, zp: int, do_plot: int = 0):
    # zero packs incoming signal, 'sig', with zero packing coefficient, 'zp', # and plots a spectra of the signal,
    # for a given sampling period (dt). # set do_plot = 0 to supress plotting # default for zp is 1, i.e no zero packing
    # If sig is a matrix with n columns, n signals assumed to have same # freq and fm will have m columns # signals must
    # be in COLUMNS # # Modified 21 / 4 / 02, PY x axis started from first frequency point # instead of 0.
    # First point should be at 0, not df.Fixed.Previously # peaks have been one frequency resolution too far to the right
    # Modified 03 / 03 / 06 PY.To be consistent with T / HIS, the magnitude of the # fft is now normalized by the number
    #  of points in it. Modified 03 / 03 / 06 PY.If you have only one signal and send do_plots as > 1 you # will get an
    #  fft with 1 / 3rd octave boundaries shaded

    # transpose if a row vector has been sent
    sig = numpy.array(sig)
    # [r, c] = sig.shape
    # if r == 1:
    #     sig = numpy.transpose(sig)
    #     print('oafft: row vector detected, transposing signal into column vector')

    # sf = 1 / dt
    n = len(sig)
    np = zp * n
    if np % 1 != 0 or np <= 0:
        print('please use a ''zp'' value that gives positive integer number of FFT points')

    ff = abs(fft.fft(sig, np))
    npf = math.ceil((np + 1) / 2)  # number of unique FFT points
    df = (1 / dt) / np
    f = df * (npf - 1)
    fm = 2 * numpy.divide(ff,
                          n)  # normalise by length of signal(NOT zeropacked npf!), but also multiply by 2 to account for
    # the symmetric half of FFT that was removed
    fm[0] = fm[0] / 2  # DC component is not repeated in "symmetric" part of FFT, so its magnitude should not be x2 above
    if np % 2 == 0:  # if np is even, the Nyquist frequency component exists, but is again not repeated in the
        # "symmetric" part of FFT and should not be x2
        fm[len(fm) - 1] = fm[len(fm) - 1] / 2

    # if do_plot != 0:
    #     plt.figure()
    #     plt.plot(f, fm)
    #     plt.xlabel('Frequency, Hz')
    #     plt.ylabel('FFT Magnitude')
    #
    # if do_plot > 1: # shade 1 / 3rd octave bands
    #     plt.figure()
    #     flow = f(2)
    #     fhigh = f(npf)
    #     # find 1 / 3 rd octave bands containing lowest and highest
    #     [fc, fl, fu] = oct3bands(flow, fhigh)
    #     nb = len(fc)
    #     for ii in range(1,nb):
    #         if (math.floor(ii / 2) == ii / 2):
    #             # even
    #             c = 0.15 # alternate between two shading colours
    #         else:
    #             # odd
    #             c = 0.1 # alternate between two shading colours
    #         nni = find(f > fl(ii) & f < fu(ii))
    #
    #         xx = [fl(ii) f(nni) fu(ii)]
    #         yy = [0.0 fm(nni) 0.0]
    #         plt.patches(xx, yy, c)

    return f, fm
