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

def dist(wavelength_indices, a, b):
    dist = 0
    
    for w in wavelength_indices:
        dist += (a[i] - b[i])**2

    return dist

def get_r_indices(wavelengths):
    return range(150,450)

def get_l_indices(wavelengths):
    return range(0,150)

def calc_all_dists(spectra, wavelength_indices):
    dists = {}
    m = spectra.shape[0]
    for i in range(m):
        dists[i] = {}
        a = spectra[i]
        for j in range(m):
            b = spectra[j]
            dists[i,j] = dist(wavelength_indices, a, b)
    return dists

def h(index, dists):
    return 

def func_reg(left_train, right_train, right_test):
    m, n = left_train.shape
    lefthat = np.zeros(m)
    ###########################


    ###########################
    return lefthat

