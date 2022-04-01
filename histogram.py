# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


# %%

def super_histogram(data, bins='auto', cumulative=0):
    N = data.size
    counts, bin_edges = np.histogram(data, bins=bins)

    if cumulative == 0:
        # RETURN PDF
        bin_widht = (bin_edges[1:] - bin_edges[:-1])
        x = (bin_edges[1:] + bin_edges[:-1]) / 2
        y = counts / (N * bin_widht)
        y_err = np.sqrt(counts * (1 - counts / N)) / (N * bin_widht)
    else:

        if cumulative == 1:
            # RETURN POSITIVE CDF
            out_left = np.sum(data < bin_edges[0])
            counts = out_left + np.cumsum(np.insert(counts, 0, 0))

        elif cumulative == -1:
            # RETURN NEGATIVE CDF
            out_right = np.sum(data > bin_edges[-1])
            counts = out_right + np.cumsum(np.append(counts, 0)[::-1])[::-1]

        x = bin_edges
        y = counts / N
        y_err = np.sqrt(counts * (1 - counts / N)) / N
    return x, y, y_err


# %%

# distr = stats.triang
# params = (0.7, 0, 1.0)

# distr = stats.uniform
# params = (0, 1.0)

# distr = stats.gamma
# params = (1.3, 0.0, 1.0)

# distr = stats.gamma
# params = (1.3, 0.0, 0.2)

# distr = stats.binom
# params = (10000, 0.4, 0.0)


# %%
distr = stats.gamma
params = (1.3, 0.0, 0.2)
count1 = []
count2 = []
count3 = []
for i in range(1000):
    data = distr.rvs(*params, size=15000)

    # edg, cdf, cdf_err = super_histogram(data, bins, -1)
    # theoretical = 1 - distr.cdf(edg, *params)

    # edg, cdf, cdf_err = super_histogram(data, bins, 1)
    # theoretical = distr.cdf(edg, *params)

    edg, cdf, cdf_err = super_histogram(data)
    theoretical = distr.pdf(edg, *params)

    count1.append(np.mean(np.abs(theoretical - cdf) < 1 * cdf_err))
    count2.append(np.mean(np.abs(theoretical - cdf) < 2 * cdf_err))
    count3.append(np.mean(np.abs(theoretical - cdf) < 3 * cdf_err))
print(100 * np.mean(count1), 100 * np.mean(count2), 100 * np.mean(count3))


# %%
distr = stats.gamma
params = (1.3, 0.0, 0.5)
data = distr.rvs(*params, size=200)
bins = np.linspace(0.5, 2.5, 30)
x = np.linspace(0, bins[:-1], 300)
fig, ax = plt.subplots(dpi=300)


edg, cdf, cdf_err = super_histogram(data, bins, -1)
ax.plot(x, 1 - distr.cdf(x, *params), 'C0:')
ax.errorbar(edg, cdf, yerr=cdf_err, fmt='C0.')

edg, cdf, cdf_err = super_histogram(data, bins, 1)
ax.plot(x, distr.cdf(x, *params), 'C1:')
ax.errorbar(edg, cdf, yerr=cdf_err, fmt='C1.')

edg, pdf, pdf_err = super_histogram(data, bins)
ax.plot(x, distr.pdf(x, *params), 'C2:')
ax.errorbar(edg, pdf, yerr=pdf_err, fmt='C2.')


# %%

def super_bincount(data, cumulative=0, n_bins=None, log_binning=False, xlim=(None, None)):
    N = data.size

    min = np.min(data)
    max = np.max(data)
    counts = np.bincount(data.flatten() - min)
    max = len(counts) + min - 1
    x = np.arange(min, min + len(counts))

    out_left = 0
    out_right = 0
    if xlim[0]:
        if xlim[0] > min:
            out_left = np.sum(counts[:xlim[0] - min])
            x = x[xlim[0] - min:]
            counts = counts[xlim[0] - min:]
            min = xlim[0]

    if xlim[1]:
        if xlim[1] < max:
            out_right = np.sum(counts[xlim[1] - max:])
            x = x[:xlim[1] - max]
            counts = counts[:xlim[1] - max]
            max = xlim[1]

    original_n_bins = len(counts)

    if not cumulative:  # COMPUTE PMF

        if n_bins:

            if log_binning:
                add_at = (n_bins * np.log(x / min) /
                          np.log(((max + 1) / min))).astype(int)
            else:
                add_at = np.linspace(
                    0, n_bins, original_n_bins, endpoint=False).astype(int)

            normalization = np.zeros(n_bins, dtype=int)
            binned_counts = np.zeros(n_bins, dtype=int)
            new_bincenter = np.zeros(n_bins, dtype=int)
            np.add.at(binned_counts, add_at, counts)
            np.add.at(new_bincenter, add_at, x)
            np.add.at(normalization, add_at, 1)

            x = new_bincenter / normalization
            y = binned_counts / (N * normalization)
            y_err = np.sqrt(
                binned_counts * (1 - binned_counts / N)) / (N * normalization)

        else:
            y = counts / N
            y_err = np.sqrt(counts * (1 - counts / N)) / N

    else:
        if cumulative == 1:
            counts = np.cumsum(counts) + out_left
        elif cumulative == -1:
            x -= 1
            counts = np.cumsum(counts[::-1])[::-1] + out_right

        y = counts / N
        y_err = np.sqrt(counts * (1 - counts / N)) / N

        if n_bins:
            if log_binning:
                bins = np.unique(np.geomspace(min, max,
                                              n_bins).astype(int)) - min
            else:
                bins = np.linspace(0, len(x) - 1, n_bins).astype(int)
            x = x[bins]
            y = y[bins]
            y_err = y_err[bins]

    return x, y, y_err


distr = stats.logser
params = (0.999, 0)

# distr = stats.binom
# params = (1000, 0.5, 0)

data = distr.rvs(*params, size=1000)
x = np.arange(1, np.max(data))
fig, ax = plt.subplots()
ax.set(xscale='log')
ax2 = ax.twinx()

bin_center, pdf, pdf_err = super_bincount(
    data, cumulative=0, n_bins=20, log_binning=True, xlim=(10, 1000))
ax.plot(x, distr.pmf(x, *params), 'C0:')
ax.errorbar(bin_center, pdf, yerr=pdf_err, fmt='C0.')

bin_center, pdf, pdf_err = super_bincount(
    data, cumulative=1, n_bins=20, log_binning=True, xlim=(10, 1000))
ax2.plot(x, distr.cdf(x, *params), 'C1:')
ax2.errorbar(bin_center, pdf, yerr=pdf_err, fmt='C1.')

bin_center, pdf, pdf_err = super_bincount(
    data, cumulative=-1, n_bins=20, log_binning=True, xlim=(10, 1000))
ax2.plot(x, 1 - distr.cdf(x, *params), 'C2:')
ax2.errorbar(bin_center, pdf, yerr=pdf_err, fmt='C2.')

# %%
