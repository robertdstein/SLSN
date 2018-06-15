import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from merge_csv import Catalogue
import numpy as np
import math


class Plotting(Catalogue):

    def __init__(self):
        Catalogue.__init__(self)
        self.plot_skymap()
        self.plot_histograms()

    def wrap_around_180(self, ra_deg):
        ra = np.deg2rad(ra_deg)
        ra[ra > np.pi] -= 2 * np.pi
        return ra

    def plot_skymap(self):
        """Plots a map of the distribution of SLSN candidates on the sky.
        Uses the redshift as a colour scale.
        """
        ra = np.array(self.data_table["RA"])
        dec = np.array(self.data_table["Dec"])
        rs = np.log(np.array(self.data_table["Redshift"]))

        plt.figure()
        plt.subplot(111, projection="aitoff")

        cm = plt.cm.get_cmap('RdYlGn_r')
        sc = plt.scatter(
            self.wrap_around_180(ra), np.deg2rad(dec), c=rs, s=35,
            cmap=cm)
        cbar = plt.colorbar(sc)
        cbar.set_label('Log(Redshift(z))')

        path = "plots/SLSN skymap.pdf"
        print "Saving to", path
        plt.tight_layout()
        plt.savefig(path)
        plt.close()

    def plot_histograms(self):
        variables = ["Redshift", 'Absolute Magnitude Peak', "Peak Date",
                     "RA", "Dec", "EBV"]

        n_rows = len(variables) / 2 + len(variables) % 2

        fig = plt.figure()

        for i, var in enumerate(variables):
            plt.subplot(n_rows, 2, i + 1)
            data = np.array(self.data_table[var])
            data = data[[not math.isnan(x) for x in data]]
            plt.hist(data, bins=30)
            plt.xlabel(var)

        fig.set_size_inches(n_rows * 4, 7)
        fig.subplots_adjust(hspace=.5)

        path = "plots/SLSN Histogram.pdf"
        print "Saving to", path
        plt.savefig(path)
        plt.close()


Plotting()
