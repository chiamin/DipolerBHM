import matplotlib.pyplot as plt
import sys, glob
from matplotlib.backends.backend_pdf import PdfPages

def plot_line (xs, ys, n, dx=0.05):
    x = xs[0]
    if n == 0:
        plt.plot (xs, ys, ls='--', c='k', zorder=0)
    elif n % 2 == 1:
        plt.plot (xs, ys, ls='-', c='k', zorder=0)
        for i in range(n//2):
            plt.plot ([x-dx*(i+1)]*2, ys, ls='-', c='k', zorder=0)
            plt.plot ([x+dx*(i+1)]*2, ys, ls='-', c='k', zorder=0)
    elif n % 2 == 0:
        for i in range(n//2):
            plt.plot ([x-dx*(i+1)+0.5*dx]*2, ys, ls='-', c='k', zorder=0)
            plt.plot ([x+dx*(i+1)-0.5*dx]*2, ys, ls='-', c='k', zorder=0)

def plot_diagram (fname):
    with open(fname) as f:
        exist = []
        for line in f:
            tmp = line.split()
            if tmp[0] == '.' or tmp[0] == '*' or tmp[0] == 'x':
                x = int(tmp[1])
                t = float(tmp[2])
                while [x,t] in exist:
                    t += 0.002
                if tmp[0] == '.':
                    plt.plot ([x], [t], marker='o', ls='None', c='k', mfc='w', zorder=100)
                elif tmp[0] == '*':
                    plt.plot ([x], [t], marker='o', ls='None', c='k', mfc='r', zorder=300)
                elif tmp[0] == 'x':
                    plt.plot ([x], [t], marker='x', ls='None', c='k', mew=2, zorder=200)
                exist.append ([x,t])
            elif tmp[0] == '-':
                x = int(tmp[1])
                t1 = float(tmp[2])
                t2 = float(tmp[3])
                n = int(tmp[4])
                plot_line ([x,x], [t1,t2], n)
            elif tmp[0] == '->':
                continue
                x1 = int(tmp[1])
                t1 = float(tmp[2])
                x2 = int(tmp[3])
                t2 = float(tmp[4])
                plt.gca().annotate('', xy=(x2,t2), xytext=(x1,t1), arrowprops=dict(color='tab:blue', width=0.1, headwidth=4))
            elif tmp[0] == '<->':
                x1 = int(tmp[1])
                t1 = float(tmp[2])
                x2 = int(tmp[3])
                t2 = float(tmp[4])
                plt.plot ([x1,x2], [t1,t2], ls='-', c='k')

def plot_all ():
    fnames = glob.glob('test*.wld')
    fnames = sorted(fnames,key=lambda i: int(i.split('test')[-1].split('.wld')[0]))
    fnames = fnames[-5:]
    with PdfPages('wld.pdf') as pdf:
        for fname in fnames:
            print (fname)
            plt.figure()
            plot_diagram(fname)
            pdf.savefig()
            #plt.close()
    plt.show()

def plot_one (fname):
    plt.figure()
    plot_diagram(fname)
    plt.show()

fname = sys.argv[1]
plot_one (fname)
