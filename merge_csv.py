import numpy as np
import csv
import string

class catalogue():
    """Class to extract SLSN entries from various data tables.
    """
    def __init__(self):
        self.entries = []
        self.extract_1208_3217()
        self.extract_1604_08207()

    def extract_1208_3217(self):
        """Extracts SLSN entries for arxiv paper 1208.3217
        """
        path = "source_lists/tabula-1208.3217v1.csv"
        arxiv = "1208.3217"
        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            type = np.nan
            for i, row in enumerate(reader):
                if i < 2:
                    pass
                elif row[1] == "":
                    type = row[0]
                else:
                    new = slsn()
                    new.name = row[0]
                    new.type = type
                    new.arxiv = arxiv
                    new.redshift = float(row[1])
                    new.abs_peak = float(row[2])

                    energy = str(row[3]).split(str(" "))
                    factor = energy[-1][2:]
                    if factor != "":
                        factor = 10 ** float(factor)

                        if len(energy) > 2:
                            down = (
                                float("".join([x for x in energy[0]][:-3])) *
                                factor)
                            up = (
                                float("".join([x for x in energy[1]][:-2])) *
                                factor)
                        else:
                            down = (
                                float("".join([x for x in energy[0]][:-2])) *
                                factor)
                            up = down

                    else:
                        up = np.nan
                        down = np.nan

                    new.radiated_energy_upper = up
                    new.radiated_energy_lower = down
                    new.ref = row[4]

                    self.entries.append(new)

    def extract_1604_08207(self):
        path = "source_lists/tabula-1604.08207.csv"
        arxiv = "1604.08207"
        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 1:
                    pass
                else:
                    new = slsn()
                    new.name = "PTF" + row[0]
                    new.ra = row[1]
                    new.dec = row[2]

                    if len(row[3]) > 2:
                        new.type = "SLSN-R"
                    else:
                        new.type = "SLSN-" + row[3]

                    new.redshift = float(row[4])
                    new.arxiv = arxiv
                    new.ref = arxiv

                    peak = [x for x in row[5]][3:]
                    while not peak[0].isdigit():
                        peak.pop(0)
                    peak = -1. * float("".join(peak))
                    self.abs_peak = peak

                    self.ebv = float(row[13])
                    date = [x for x in row[10]]
                    if len(date) > 10:
                        date.pop(0)
                        print date, "Changed!"
                    print [x for x in row[10]]


                    # print vars(new)
                    # print row
                    # print vars(new)



class slsn:
    """Class for one superluminous supernovae
    """
    def __init__(self):
        self.name = np.nan
        self.type = np.nan
        self.redshift = np.nan
        self.abs_peak = np.nan
        self.radiated_energy_upper = np.nan
        self.radiated_energy_lower = np.nan
        self.ref = np.nan
        self.arxiv = np.nan
        self.peak_date = np.nan
        self.ra = np.nan
        self.dec = np.nan
        self.ebv = np.nan

catalogue()