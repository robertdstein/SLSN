import numpy as np
import csv
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from tabulate import tabulate

class Catalogue():
    """Class to extract SLSN entries from various data tables."""
    def __init__(self):
        self.entries = []
        self.data_table = np.nan

        self.extract_1208_3217()
        self.extract_1604_08207()
        self.extract_1605_0250()
        self.extract_1612_05978()
        self.extract_1612_07321()
        self.extract_1705_06047()
        self.extract_1708_01623()
        self.extract_1708_08971()
        self.extract_the_open_supernova_catalog()

        self.unique_entries = []
        self.combine_entries()

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
            ("Peak Date", np.float),
            ("Discover Date", np.float),
            ("RA", np.float),
            ('Dec', np.float),
            ("EBV", np.float),
            ("Notes", "S50"),
            ("Alias", "S50")
        ])

        # Ignore entries if the ONLY source is one of the following (because
        # their selection criteria are not reliable).
        ignore_list = ["1708.08971 (p.11)"]
        # ignore_list = []

        table = np.zeros_like(self.unique_entries, dtype=dt)
        mask = []

        for i, sn in enumerate(self.unique_entries):
            table[i] = np.array(
                (sn.name, sn.type, float(sn.redshift), sn.abs_peak,
                 sn.radiated_energy_lower, sn.radiated_energy_upper, sn.ref,
                 str(sn.arxiv), float(str(sn.peak_date)),
                 float(str(sn.disc_date)),
                 float(sn.ra), float(sn.dec), sn.ebv,
                 sn.notes, str(sn.alias)),
                dtype=dt
            )

            mask.append(sn.arxiv not in ignore_list)

            if sn.notes != "":
                print "\n", sn.name, "\n", sn.notes

        mask = np.array(mask)
        self.unique_entries = np.array(self.unique_entries)

        to_print = ["Name", "Alias", "Redshift", "Type", "RA", "Dec",
                    # "Absolute Magnitude Peak", "Peak Date", "Discover Date",
                    "Arxiv Link"]

        self.data_table = np.sort(
            table[mask], order=['Name'], axis=0).view()
        print tabulate(self.data_table[to_print], to_print)

        print "There are", len(self.entries), "entries in SLSN Catalogue."
        print "There are", len(self.unique_entries), "unique entries in SLSN " \
                                                     "Catalogue."

        print "We reject any SLSN which are listed only in one of the folowing",
        print "references, on the basis that these sources are unreliable: \n"

        for link in ignore_list:
            print " -", link

        print ""
        print "This leaves us with", len(self.unique_entries[mask]),
        print "reliable SLSN entries in the Catalogue."

        print "We require category to be SLSN"

    def combine_entries(self):
        """Combine entries from

        :return:
        """

        def convert_string(name):
            """Convert a name into the most basic form, to eable comparison
            despite differing style conventions.

            :param name: Name to be converted
            :return: Name, in lowercase, with spaces, brackets, colons and
            dashes removed.
            """
            new = name.lower()
            for character in [" ", "-",  "(", ")", ":"]:
                new = new.replace(character, "")
            return new

        for entry in self.entries:
            included = False
            new_names = [entry.name]
            new_names.extend(entry.alias)
            for name in new_names:
                name_to_check = convert_string(name)

                for existing_entry in self.unique_entries:
                    current_names = [existing_entry.name]
                    current_names.extend(existing_entry.alias)
                    for current in current_names:
                        if included:
                            break
                        if name_to_check == convert_string(current):
                            included = True
                            notes = (entry.notes + ", " +
                                     existing_entry.notes + ", ")
                            for key in vars(entry).keys():
                                new_val = getattr(entry, key)
                                old_val = getattr(existing_entry, key)
                                if new_val is not np.nan:
                                    if key in ["arxiv", "ref"]:
                                        to_add = old_val + ", " + new_val

                                        setattr(existing_entry, key, to_add)
                                    elif key == "alias":
                                        all_names = list(new_names)
                                        all_names.extend(current_names)
                                        all_set = set(all_names)
                                        all_names = list(all_set)
                                        all_names.remove(existing_entry.name)

                                        setattr(existing_entry, key,
                                                all_names)

                                    elif key in ["name", "notes", "coord"]:
                                        pass

                                    elif old_val is np.nan:
                                        setattr(existing_entry, key,
                                                new_val)
                                    elif old_val == new_val:
                                        pass
                                    elif (isinstance(old_val, float) and
                                              isinstance(new_val, float)):
                                        n_digits = min(len(str(old_val)),
                                                       len(str(new_val))) - \
                                                   (len(str(int(old_val))) + 1)

                                        if round(old_val, n_digits) == round(
                                                new_val, n_digits):
                                            if len(str(old_val)) > len(
                                                    str(new_val)):
                                                pass
                                            else:
                                                setattr(existing_entry, key,
                                                        new_val)
                                        else:
                                            notes += ("Numerical Discrepancy "
                                                      "in " + key +
                                                      ": " + str(old_val) +
                                                      ", " + str(new_val) +
                                                      ", ")
                                    else:
                                        notes += ("Discrepancy in " + key +
                                                  ": " + str(old_val) +
                                                  ", " + str(new_val) + ", ")
                            setattr(existing_entry, "notes", notes)
                        # break


            if not included:
                self.unique_entries.append(entry)

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
                    new = SLSN()
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
                    new = SLSN()
                    new.name = "PTF" + row[0]
                    ra = row[1]
                    dec = [x for x in row[2]]

                    if len(dec) > 12:
                        dec.pop(0)
                        dec.pop(0)
                        dec[0] = "-"

                    dec = "".join(dec)

                    new.add_coordinates(ra, dec)

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
                    date.format = "mjd"
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
                    new = SLSN()
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

        # The data table is spread over two pages, so each csv is looped over

        path1 = "source_lists/tabula-1612.05978 (dragged)(1).csv"
        path2 = "source_lists/tabula-1612.05978-1 (dragged) 2(1).csv"
        paths = [path1, path2]

        # Additional aliases are given at the bottom of the table as footnotes.
        # The alias file contains these, in txt format

        alias_dict = dict()
        alias_path = "source_lists/1612.05978_aliases.txt"

        # The reference for each entry is given at the bottom of the table.
        # The ref file contains these, in txt format

        ref_dict = dict()
        ref_path = "source_lists/1612.05978_references.txt"

        def add_to_alias_dict(current_list, counter):
            """Adds an entry to the Alias dictionary, which may itself
            contain more than one name within a list.

            :param current_list: Contains the list of all aliases
            :param counter: The index of the footnote
            """
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

        # Creates the alias dictionary

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

        # Creates the references dictionary

        with open(ref_path, 'rb') as f:
            reader = csv.reader(f, delimiter=';', quotechar='|')
            for row in reader:
                for entry in row[1:-1]:
                    entry = entry.split(":")
                    key = entry[0].strip()[1:-1]
                    value = entry[1]
                    ref_dict[str(key)] = value

        # Creates a list to check each footnote has been found in reading

        remainders = []

        for path in paths:
            with open(path, 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter=',', quotechar='|')
                for row in reader:
                    new = SLSN()

                    # Checks for the footnote numbers, which are added to the
                    #  end of the entries in this column. Also removes
                    # special characters from names.

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

                    # Adds footnote to remainder list, to check that all 24
                    # footnotes were found and removed from the entries. Adds
                    #  the corresponding aliases to the entry

                    if rest.isdigit():
                        remainders.append(int(rest))
                        new.alias = alias_dict[rest]

                    ra = row[1]

                    # Checks for typo in SN2213-1745 RA value, with "." rather
                    # than ":" between seconds

                    if len(ra.split(":")) == 2:
                        ra_split = ra.split(".")
                        ra = ".".join([":".join(ra_split[:2]), ra_split[2]])

                    # Assigns dec, and checks to see if a "-" hasbeen ommited
                    #  from negative declinations
                    dec = row[2]
                    if len(dec) < 11:
                        dec = "-" + dec
                    dec = dec
                    new.add_coordinates(ra, dec)

                    # Removes special characters from redshift

                    rs = [x for x in row[3] if x.isdigit() or x == "."]
                    new.redshift = float("".join(rs))

                    new.type = row[4]
                    new.ebv = float(row[5])

                    # Finds the numbers for the references of each entry,
                    # looks up the corresponding references themselves,
                    # and adds these to the entry.

                    refs = ",".join(row[7:]).strip('"')[1:-1].split(",")

                    references = ""

                    for i in refs:
                        references += ref_dict[i.strip()] + ","

                    new.arxiv = arxiv

                    self.entries.append(new)

        if len(remainders) != 24:
            raise Exception("Warning! There should be 24 footnote-related "
                            "errors. However, there were not 24 "
                            "name-modification corrections.")

    def extract_1612_07321(self):
        """Extracts SLSN entries for arxiv paper 1612.07321 (page 4)"""
        arxiv = "1612.07321 (p. 4)"

        new_entries = []

        path = "source_lists/tabula-1612.07321.csv"

        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 1:
                    pass
                elif row[0][0] != '"':
                    new = SLSN()
                    name = row[0].split("/")

                    # Add note for SN2011kl, which is associated with an
                    # ultra-long GRB
                    if len(name) > 1:
                        GRB_name = "".join([x for x in name[1] if x is not "?"])
                        new.notes = "Associated with " + GRB_name
                        new.name = name[0]

                    # Adds note for ASASSN-15lh, which has a debatable
                    # categorisation as either an SLSN or a TDE
                    elif "?" in name[0]:
                        new.notes = "May be a Tidal Disruption Event (TDE)!"
                        new.name = "".join([x for x in name[0] if x is not "?"])
                    else:
                        new.name = name[0]

                    new.redshift = float(row[1])
                    new.ref = row[-1]
                    new.arxiv = arxiv
                    new.type = "SLSN-I"
                    new_entries.append(new)
                else:
                    previous = new_entries[-1]

                    # Checks for additional aliases, given in the row beneath
                    #  the main entry
                    if row[0] != '""':
                        new_alias = []
                        j = 0
                        while row[j] != '':
                            alias = row[j]
                            alias = alias.replace("(or ", "")
                            alias = alias.replace(")", "")
                            alias = alias.replace('"', "")

                            new_alias.append(alias)
                            j += 1
                        previous.alias = new_alias

                    # Checks for additional references, given in the row
                    # beneath the main entry
                    if row[-1] != '':
                        previous.ref += ", " + row[-1]

        self.entries.extend(new_entries)

    def extract_1705_06047(self):
        """Extracts SLSN entries for arxiv paper 1705.06047 (page 7)"""
        arxiv = "1705.06047 (p.7)"

        new_entries = []

        path = "source_lists/tabula-1705.06047.csv"

        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 2:
                    pass
                else:
                    new = SLSN()
                    new.name = row[0]
                    new.redshift = float(row[1])
                    new.ref = row[3]
                    new.arxiv = arxiv
                    new.type = "SLSN-I"
                    lum = float(row[2]) * 10 ** 44
                    lum = lum * u.erg / u.second
                    new.abs_peak = 4.77 - 2.5 * np.log10(lum/u.L_sun.cgs)

                    new_entries.append(new)

        self.entries.extend(new_entries)

    def extract_1708_01623(self):
        """Extracts SLSN entries for arxiv paper 1708.01623 (page 3)"""
        path = "source_lists/tabula-1708.01623 (dragged).csv"
        arxiv = "1708.01623 (p.3)"
        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 2:
                    pass
                else:
                    # print row
                    new = SLSN()
                    new.name = "PTF" + row[0]
                    new.arxiv = arxiv
                    ra = row[1]

                    dec = [x for x in row[2]]

                    if len(dec) > 12:
                        dec.pop(0)
                        dec.pop(0)
                        dec[0] = "-"

                    dec = "".join(dec)

                    new.add_coordinates(ra, dec)

                    new.redshift = float(row[3])
                    new.type = "SLSN-" + row[4]
                    new.ref = arxiv
                    self.entries.append(new)

    def extract_1708_08971(self):
        """Extract SLSN entries for arxiv paper 1708.08971 (page 11)"""
        path = "source_lists/tabula-1708.08971 (dragged).csv"
        arxiv = "1708.08971 (p.11)"

        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for i, row in enumerate(reader):
                if i < 1:
                    pass
                else:
                    new = SLSN()
                    new.name = row[0]
                    new.redshift = float(row[5])
                    new.ref = row[8]
                    new.arxiv = arxiv

                    date_type = row[6][-1]
                    mjd = Time(float(row[7]), format="mjd")
                    mjd.out_subfmt = "date"

                    ra = row[1]
                    dec = row[2]
                    new.add_coordinates(ra, dec)

                    if date_type == "p":
                        new.peak_date = mjd
                    elif date_type == "d":
                        new.disc_date = mjd
                    else:
                        raise Exception("Unknown date type!")

                    self.entries.append(new)

    def extract_the_open_supernova_catalog(self):
        """Extracts SLSN from the Open supernovae catalog with the type SLSN.
        The catalog can be found at https://sne.space/ .
        """
        arxiv = "1605.01054 (https://sne.space/)"

        new_entries = []

        path = "source_lists/The Open Supernova Catalog.csv"

        with open(path, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            for i, row in enumerate(reader):
                if i < 1:
                    pass

                # Check for weird artifact entry in SLSN catalog
                elif row[0] != "PISN_Jerkstrand":
                    new = SLSN()
                    new.name = row[0]
                    new.alias = list(
                        set([x.strip() for x in row[1].split(",")
                            if x != new.name]))

                    for j, attr in enumerate(["disc_date", "peak_date"]):
                        try:
                            date = "-".join(row[j+2].split("/"))
                            date = Time(date, format="iso")
                            date.out_subfmt = "date"
                            date.format = "mjd"
                            setattr(new, attr, date)
                        except ValueError:
                            if row[j + 2] != "":
                                new.notes += "Problem in " + attr + ", only have "
                                new.notes += row[j + 2] + " (not ISO format), "

                    if row[4] != "":
                        new.abs_peak = float(row[4])

                    coords = [np.nan, np.nan]

                    for j, angle in enumerate(coords):
                        if row[5 + j] != "":
                            values = row[5 + j].split(",")
                            coords[j] = values[0]
                            if len(values) > 1:
                                new.notes += "Several values reported for " + \
                                             ["ra", "dec"][j] + " " + \
                                             str(values) + ", "

                    new.add_coordinates(coords[0], coords[1])

                    for j, attr in enumerate(["redshift", "ebv"]):
                        val = row[7 + (2 * j)]
                        if val != "":
                            all_vals = list(set(val.split(",")))
                            setattr(new, attr, float(all_vals[0]))
                            if len(all_vals) > 1:
                                new.notes += "Several values reported for " + \
                                             attr + " " + str(all_vals) + ", "

                    type_val = row[8]
                    if type_val != "":
                        all_vals = list(set(type_val.split(",")))
                        setattr(new, "type", all_vals[0])
                        if len(all_vals) > 1:
                            new.notes += "Several values reported for " + \
                                         attr + " " + str(all_vals) + ", "

                    ref_val = row[10]
                    if ref_val != "":
                        setattr(new, "ref", ref_val)



                    new.arxiv = arxiv
                    new_entries.append(new)

        self.entries.extend(new_entries)

class SLSN:
    """Class for one superluminous supernovae
    """
    def __init__(self):
        self.name = np.nan
        self.type = np.nan
        self.redshift = np.nan
        self.abs_peak = np.nan
        self.radiated_energy_upper = np.nan
        self.radiated_energy_lower = np.nan
        self.ref = ""
        self.arxiv = ""
        self.peak_date = np.nan
        self.disc_date = np.nan
        self.ra = np.nan
        self.dec = np.nan
        self.coord = np.nan
        self.ebv = np.nan
        self.notes = ""
        self.alias = []

    def add_coordinates(self, ra, dec):
        if np.nan not in [ra, dec]:
            split_ra = ra.split(":")
            if len(split_ra) == 1:
                ra_fmt = u.deg
            else:
                ra_fmt = u.hourangle

            split_dec = dec.split(":")
            if len(split_dec) == 1:
                dec_fmt = u.deg
            else:
                dec_fmt = u.deg

            c = SkyCoord(ra + " " + dec, unit=(ra_fmt, dec_fmt))
            self.ra = c.ra.deg
            self.dec = c.dec.deg
            self.coord = c

if __name__ == '__main__':
    Catalogue()

