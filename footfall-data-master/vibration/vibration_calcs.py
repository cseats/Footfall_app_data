import math
import numpy
import json
import scipy
import scipy.fftpack
import matplotlib.pyplot as plt


def process(data):
    fails = 0
    for record in data:
        l = len(record)
        times, accels = [], []
        for point in record:
            times.append(point['time'])
            accels.append(point['acceleration'])

        dt = times[l - 1] - times[l - 2]

        if dt == 0:
            fails += 1
            continue

        fs = round(1 / dt)

        if fs < 80 or fs > 120:
            fails += 1
            continue

        a_z = 7
        b_z = 7
        a_z_hp = 7
        b_z_hp = 7

        freq_warning = False

        # my_data = filter_this(record, b_z_hp, a_z_hp)
        tol_df = 0.2
        target_df = 0.1
        n_pt = 0
        n_win = 0
        z_pt = 0
        if (1 / dt) < target_df * tol_df:
            fft_type = 2
            n_pt = math.floor((1 / target_df) / dt)
            n_win = math.floor(len(record) / n_pt)
            overlap = 0.5
        elif (1 / dt) > target_df * tol_df:
            fft_type = 0
            z_pt = math.floor(1 / (target_df * dt))
            n_pt = len(record)
        else:
            fft_type = 1
            n_pt = len(record)

            out_pt = math.ceil(z_pt / 2)

        print(fft_type)
        if fft_type == 1:
            transformed = scipy.fftpack.fft(accels)
            scaled = [x / n_pt for x in transformed]
            plt.plot(scaled)
        elif fft_type == 0:
            transformed = scipy.fftpack.fft(accels)
            scaled = [x / n_pt for x in transformed]
            plt.plot(scaled)
        else:
            scipy.fftpack.fft(accels)
            scaled = [x / n_pt for x in transformed]
            plt.plot(scaled)

    print(fails, 'recordings with incorrect sampling rates were ignored.')
    plt.show()


def ns_rms(data, dt, t):
    samples_in_one_t = math.floor(t / dt)
    n = 1
    for d in data:
        d = d ** 2
        d = numpy.cumsum(d) / samples_in_one_t
        d = math.sqrt(d)
        n += 1
    return data


def filter_this(record, b, a):
    z = 1
    for point in record:
        temp = b * record[1] * z

    data_hp = None

    return z


if __name__ == '__main__':
    with open('../clean_data.json') as data_file:
        clean_data = json.load(data_file)
    process(clean_data)
