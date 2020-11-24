#!/usr/bin/env python3
import re


class Paper:
    def __init__(self, record):
        self.pubmed_id = record['PMID']
        self.title = record['TI']
        self.abstract = record['AB']
        self.journal = record['TA']
        self.doi = record['LID'].split()[0]
        self.dp = record['DP']
        self.so = record['SO']
        self.fau = record['FAU']
        self.au = record['AU']
        self.ad = record['AD']

    @property
    def authors(self):
        authors_list = []
        for fau, au, ad in zip(self.fau, self.au, self.ad):
            tokes = fau.split(", ")
            first = tokes[0]
            if len(tokes) == 2:
                last = tokes[1]
            else:
                last = first
            authors_list.append(
                {
                    "name": au,
                    "firstname": first,
                    "lastname": last,
                    "affiliation": ad
                }
            )
        return authors_list

    @property
    def publication_date(self):
        date = re.compile(". ({}[^;]*)[;.]".format(self.dp))

        return date.findall(self.so)[0]


class Record(dict):
    """A dictionary holding information from a Medline record.

    All data are stored under the mnemonic appearing in the Medline
    file. These mnemonics have the following interpretations:

    """
    def is_valuable(self):
        if 'PMID' not in self:
            return False
        if 'AB' not in self:
            return False
        if 'AU' not in self:
            return False
        if not len(self['AU']) == len(self['AD']):
            return False
        return True


def parse(handle):
    textkeys = (
        "ID",
        "PMID",
        "SO",
        "RF",
        "NI",
        "JC",
        "TA",
        "IS",
        "CY",
        "TT",
        "CA",
        "IP",
        "VI",
        "DP",
        "YR",
        "PG",
        "LID",
        "DA",
        "LR",
        "OWN",
        "STAT",
        "DCOM",
        "PUBM",
        "DEP",
        "PL",
        "JID",
        "SB",
        "PMC",
        "EDAT",
        "MHDA",
        "PST",
        "AB",
        "EA",
        "TI",
        "JT",
    )
    handle = iter(handle)

    key = ""
    record = Record()
    history = []
    for line in handle:
        line = line.rstrip()
        if line[:6] == "      ":  # continuation line
            if key in ["MH", "AD"]:
                # Multi-line MESH term, want to append to last entry in list
                record[key][-1] += line[5:]  # including space using line[5:]
            else:
                record[key].append(line[6:])
        elif line:
            key = line[:4].rstrip()
            if key not in record:
                record[key] = []
            if key == 'AD' and history[-1] == 'AD':
                record[key][-1] += " " + line[6:]
            else:
                record[key].append(line[6:])
            if (key == 'FAU' or key == 'LA') and history[-1] == 'AU':
                if 'AD' in record:
                    record['AD'].append('NA')
                else:
                    record['AD'] = ['NA']
            if key in ('PMID', 'FAU', 'AU', 'AD', 'LA'):
                history.append(key)
        elif record:
            # Join each list of strings into one string.
            for key in record:
                if key in textkeys:
                    record[key] = " ".join(record[key])
            if record.is_valuable():
                yield Paper(record)
            record = Record()
            history = []
    if record:  # catch last one
        for key in record:
            if key in textkeys:
                record[key] = " ".join(record[key])
        if record.is_valuable():
            yield Paper(record)
