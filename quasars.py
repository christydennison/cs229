from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def load_data():
    train = np.genfromtxt('quasar_train.csv', skip_header=True, delimiter=',')
    test = np.genfromtxt('quasar_test.csv', skip_header=True, delimiter=',')
    wavelengths = np.genfromtxt('quasar_train.csv', skip_header=False, delimiter=',')[0]
    return train, test, wavelengths

def add_intercept(X_):
    X = np.insert(X_, [0], np.ones((X_.shape[0], 1)), axis=1)
    return X

def smooth_data(raw, wavelengths, tau):
    X = add_intercept(wavelengths)
    smooth = np.zeros_like(raw)
    m = raw.shape[0]

    Ws = get_Ws(X, tau)

    for i in range(m):
        spectrum = raw[i]
        smoothed_spectrum = smooth_spectrum(spectrum, X, Ws)
        smooth[i] = smoothed_spectrum[:,1]
    
    return smooth

def LWR_smooth(spectrum, wavelengths, tau):
    X = add_intercept(wavelengths)
    Y = spectrum
    
    Ws = get_Ws(X, tau)
    return smooth_spectrum(Y, X, Ws)

def smooth_spectrum(spectrum, wavelengths, Ws):
    X = wavelengths
    Y = spectrum

    smoothed_spectrum = np.zeros_like(X)
    m = X.shape[0]
    for i in range(m):
        W = np.diag(Ws[i])
        Xt = np.transpose(X)
        XtW = np.dot(Xt, W)
        XtWX = np.dot(XtW, X)
        XtWy = np.dot(np.dot(Xt, W), Y)
        theta = np.dot(np.linalg.inv(XtWX), XtWy)
        smoothed_spectrum[i] = np.dot(X[i], theta)
    return smoothed_spectrum

def get_Ws(wavelengths, tau):
    m = wavelengths.shape[0]
    Ws = np.zeros((m,m))
    for i in range(m):
        for j in range(m):
            x = wavelengths[i,1]
            x_i = wavelengths[j,1]
            Ws[i,j] = np.exp(-(x-x_i)**2/(2*tau**2))
    return Ws

def LR_smooth(Y, X_):
    X = add_intercept(X_)
    theta = np.zeros(2)
    
    # normal eq: theta = (X^T X)^-1 * X^T y
    Xt = np.transpose(X)
    Xty = np.dot(Xt, Y)
    theta = np.dot(np.linalg.inv(np.dot(Xt, X)), Xty)
    yhat = np.dot(X, theta)
    return yhat, theta

def plot_b(X, raw_Y, Ys, desc, filename):
    plt.figure()
    plt.ylabel('Flux')
    plt.xlabel('Wavelengths')
    plt.scatter(X, raw_Y)
    plt.plot(X, Ys, c='r', linewidth=2.0)
    #plt.plot(X, np.asarray(Ys).reshape((X.shape[0],1)), c='r', linewidth=2.0)
    plt.suptitle(desc)
    plt.show()
    #plt.savefig(filename)

def plot_c(Yhat, Y, X, filename):
    plt.figure()
    ############


    #############
    plt.savefig(filename)
    return

def split(full, wavelengths):
    left, right = None, None
    ###############

    ###############
    return left, right

def dist(a, b):
    dist = 0

    for i in range(len(a)):
        dist += (a[i] - b[i])**2

    return dist

def calc_dists(spectra_train, spectra_test):
    m = spectra.shape[0]
    dists = np.zeros((m,m))
    for i in range(m):
        a = spectra[i]
        for j in range(m):
            b = spectra[j]
            dists[i,j] = dist(a, b)
    return dists

def k_neighbors(index, dists, k):
    return np.argsort(dists[index])[:k]

def get_h(index, dists):
    return np.max(dists[index])

def ker(t):
    return np.max(t-1, 0)

def fhat_left(smoothed_spectra, test_right):
    rights = smoothed_spectra[:, 150:450]
    lefts = smoothed_spectra[:, 0:50]
    r_dists = calc_dists(rights)
    m = smoothed_spectra.shape[0]
    n = smoothed_spectra.shape[1]
    lefthats = np.zeros((m,50))
    errors = np.zeros(m)
    # go through all training examples
    for i in range(m):
        # for each neighbor
        top_f_left = 0
        bottom_f_left = 0
        h = get_h(i, r_dists)
        for n_index in k_neighbors(i, r_dists, 3):
            nl = lefts[n_index]
            d = r_dists[i, n_index]
            k_score = ker(d/h)
            top_f_left += k_score*nl
            bottom_f_left += k_score
        fhat_i = top_f_left/bottom_f_left
        lefthats[i] = fhat_i
        errors[i] = dist(lefts[i], fhat_i)
    mean_error = np.mean(errors)
    return lefthats, mean_error

def func_reg(left_train, right_train, right_test):
    m, n = left_train.shape
    lefthat = np.zeros(m)
    ###########################


    ###########################
    return lefthat

