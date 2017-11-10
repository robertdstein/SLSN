import numpy as np
import csv
from astropy.time import Time
from tabulate import tabulate

class catalogue():
    """Class to extract SLSN entries from various data tables."""
    def __init__(self):
        self.entries = []
        self.data_table = np.nan

        self.extract_1208_3217()
        self.extract_1604_08207()
        self.extract_1605_0250()
        self.extract_1612_05978()

        self.make_table()

    def make_table(self):
        dt = np.dtype([
            ('Name', "S50"),
            ("Type", "S10"),
            ("Redshift", np.float),
            ('Absolute Magnitude Peak', np.float),
            ("Radiated Energy (Lower Limit)", np.float),
            ("Radiated Energy (Upper Limit)", np.float),
            ("Reference", "S50"),
            ("Arxiv Link", "S50"),
            ("Peak Date", "S50"),
            ("RA", "S20"),
            ('Dec', "S20"),
            ("EBV", np.float),
            ("Notes", "S50"),
            ("Alias", "S50")

        ])

        table = np.zeros_like(self.entries, dtype=dt)
        print "There are", len(self.entries), "entries in SLSN catalogue."

        for i, sn in enumerate(self.entries):
            table[i] = np.array(
                (sn.name, sn.type, sn.redshift, sn.abs_peak,
                 sn.radiated_energy_lower, sn.radiated_energy_upper, sn.ref,
                 str(sn.arxiv), str(sn.peak_date), sn.ra, sn.dec, sn.ebv,
                 sn.notes, str(sn.alias)),
                dtype=dt
            )

        to_print = ["Name", "Alias", "Redshift", "Type", "RA", "Dec",
                    "Arxiv Link"]
        table = np.sort(table, order=['Redshift'], axis=0).view()
        # print tabulate(table[to_print], to_print)

    def extract_1208_3217(self):
        """Extracts SLSN entries for arxiv paper 1208.3217 (page 31)"""
        path = "source_lists/tabula-1208.3217v1.csv"
        arxiv = "1208.3217 (p.31)"
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
        """Extracts SLSN entries for arxiv paper 1604.08207 (page 3)"""
        path = "source_lists/tabula-1604.08207.csv"
        arxiv = "1604.08207 (p.3)"
        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 1:
                    pass
                else:
                    new = slsn()
                    new.name = "PTF" + row[0]
                    new.ra = row[1]
                    dec = [x for x in row[2]]

                    if len(dec) > 12:
                        dec.pop(0)
                        dec.pop(0)
                        dec[0] = "-"

                    dec = "".join(dec)

                    new.dec = dec

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
                    new.abs_peak = peak

                    date = [x for x in row[10]]
                    if len(date) > 10:
                        date.pop(0)
                    date = "".join(date)
                    date = Time(date, format="iso")
                    date.out_subfmt = "date"
                    new.peak_date = date

                    new.ebv = float(row[13])

                    notes = row[16]
                    if notes != "":
                        if notes[0] == "=":
                            split = notes.split(".")
                            new.alias.append(split[0][2:])
                            if len(split) > 1:
                                notes = split[1]
                            else:
                                notes = ""
                        else:
                            notes = [x for x in row[16]]
                            while not notes[-1].isalpha():
                                notes.pop(-1)
                            notes = "".join(notes)
                    new.notes = notes

                    self.entries.append(new)

    def extract_1605_0250(self):
        """Extracts SLSN entries for arxiv paper 1605.01250 (pages 5 and 7)"""
        arxiv = "1605.0250 (p.5&7)"

        new_entries = []

        path1 = "source_lists/tabula-1605.05250 (dragged) 1.csv"

        with open(path1, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 1:
                    pass
                elif row[0] == '""':
                    pass
                else:
                    new = slsn()
                    new.name = row[0]
                    new.redshift = float(row[1])
                    new.ref = row[3]
                    new.type = "SLSN-I"
                    new.arxiv = arxiv
                    new_entries.append(new)
        self.entries.extend(new_entries)

    def extract_1612_05978(self):
        """Extracts SLSN entries for arxiv paper 1612.05978 (pages 4 and 5)"""
        arxiv = "1612.05978 (p.4&5)"

        path1 = "source_lists/tabula-1612.05978 (dragged)(1).csv"
        path2 = "source_lists/tabula-1612.05978-1 (dragged) 2(1).csv"
        paths = [path1, path2]

        alias_dict = dict()
        alias_path = "source_lists/1612.05978_aliases.txt"

        def add_to_alias_dict(current_list, counter):
            current_list = current_list.split(",")

            for j, name in enumerate(current_list):
                while name[0] == " ":
                    name = name[1:]
                while name[-1] == " ":
                    name = name[:-2]
                name = name.split(" ")
                name = [s for s in name if s is not ""]
                name = " ".join(name)
                current_list[j] = name

            alias_dict[str(counter)] = current_list

        with open(alias_path, 'rb') as f:
            reader = csv.reader(f, delimiter=';', quotechar='|')
            counter = 1
            current_list = ""
            for i, row in enumerate(reader):
                if i < 2:
                    pass
                elif row[0].isdigit():
                    add_to_alias_dict(current_list, counter)
                    current_list = ""
                    counter += 1
                else:
                    current_list += " ".join(row)
                    current_list += " "
            add_to_alias_dict(current_list, counter)

        remainders = []

        for path in paths:
            with open(path, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in reader:
                    new = slsn()

                    name = row[0]
                    if "SN20" in name:
                        name = name[:8]
                    elif "SSS" in name:
                        name = name[:9]
                    elif "CSS" in name:
                        name = name[:9]
                    elif "LSQ" in name:
                        name = name[:8]
                    elif "MLS" in name:
                        name = name[:9]
                    while not name[-1].isalpha() and not name[-1].isdigit():
                        name = name[:-1]
                    new.name = name
                    rest = row[0][len(name):]
                    if rest.isdigit():
                        remainders.append(int(rest))
                        new.alias = alias_dict[rest]

                    new.ra = row[1]

                    dec = row[2]
                    if len(dec) < 11:
                        dec = "-" + dec
                    new.dec = dec

                    rs = [x for x in row[3] if x.isdigit() or x == "."]
                    new.redshift = float("".join(rs))

                    new.type = row[4]
                    new.ebv = row[5]

                    refs = list("".join(row[7:]))
                    print refs

                    new.arxiv = arxiv

                    self.entries.append(new)

        if len(remainders) != 24:
            raise Exception("Warning! There should be 24 footnote-related "
                            "errors. However, there were not 24 "
                            "name-modification corrections.")



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
        self.notes = ""
        self.alias = []

catalogue()