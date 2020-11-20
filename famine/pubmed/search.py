#!/usr/bin/env python3
import os
import re
from Bio import Entrez
from medline import parse
from topic import TopicPred

Entrez.email = "liqiming1914658215@gmail.com"
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"

DIR = os.path.dirname(os.path.abspath(__file__))

chineseLast = set(line.rstrip() for line in open(os.path.join(DIR, 'ChineseFamily.csv')))
chineseFirst = set(line.rstrip() for line in open(os.path.join(DIR, 'ChineseFirst.csv')))
journal_if = list(line.rstrip().split(',') for line in open(os.path.join(DIR, 'Journal_if.csv')))
tp = TopicPred(os.path.join(DIR, 'data', 'word_dict.pikle'), os.path.join(DIR, 'data', 'train_model.m'))


def get_topic(journal, abstract):
    if 'cell' in journal or 'immunology' in journal:
        topic = "Biology"
    elif 'physics' in journal:
        topic = "Physics"
    elif 'chemistry' in journal:
        topic = "Chemistry"
    elif 'robotics' in journal:
        topic = "Computer science"
    else:
        topic = tp.classify(abstract)

    return topic


def isChinese(first, last):
    first = first.lower().replace("-", "")
    last = last.lower()
    if last in chineseLast and first in chineseFirst:
        return 'Yes'
    else:
        return 'No'


def Chines_authors(au):
    authors = list()
    au_email = list()
    Chinese_num = 0
    flag = 0
    for author in au:
        flag += 1
        first = author["firstname"]
        last = author["lastname"]
        name = author['name']
        chinese = isChinese(first, last)
        if chinese == 'Yes':
            Chinese_num += 1
        if flag == 3 and Chinese_num == 0:
            break
        affiliation = author["affiliation"]
        search_email = re.findall(r"\S+\@\S+", affiliation)
        if len(set(search_email)) > 1:
            affiliation = affiliation.replace(" ".join(search_email), "")
            au_email.append(name)
            try:
                email = search_email[au_email.index(name)]
            except IndexError:
                email = None
        elif len(set(search_email)) == 1:
            email = search_email[0]
            affiliation = affiliation.replace(email, "")
            affiliation = affiliation.replace(" Electronic address:", "")
            affiliation = affiliation.rstrip()
        else:
            email = None
        if email is not None and email[-1] == '.':
            email = email[:-1]
        authors.append([name, first, last, chinese, affiliation, email])

    return authors


class Pubmed:
    def __init__(self, Journal, min_date, max_date, retmax=10000):
        self.journal = Journal
        self.min_date = min_date
        self.max_date = max_date
        self.retmax = retmax

    @property
    def query(self):
        return "{}[ta] AND {}: {}[dp]".format(
            self.journal, self.min_date, self.max_date)

    def search(self):
        handle = Entrez.esearch(db="pubmed", term=self.query, retmax=self.retmax)
        record = Entrez.read(handle)
        return record["IdList"]

    def detail(self):
        idlist = self.search()
        print("Search '{}', get {} records.".format(self.query, len(idlist)))
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
        records = parse(handle)

        for record in records:
            authors = Chines_authors(record.authors)
            if not authors:
                continue
            topic = get_topic(record.journal, record.abstract)

            article = [record.pubmed_id, record.publication_date,
                       record.title, record.abstract,
                       record.doi, topic]
            yield authors, article


if __name__ == "__main__":
    journal_list = ('Nature', 'Cell', 'Science')
    for journal in journal_list:
        for _, article in Pubmed(journal, '2019/01/01', '2020/11/20').detail():
            print(article[0])
