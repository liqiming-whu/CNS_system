import os
import re
import xml.etree.ElementTree as ET
from Bio import Entrez
from pymed.article import PubMedArticle
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


def replace_tags(text, reverse=False):
    tags = ['<sup>', '</sup>', '<sub>', '</sub>', '<i>', '</i>']
    re_tags = ['[^', '^]', '[_', '_]', '[/', '/]']
    for tag, re_tag in zip(tags, re_tags):
        if reverse:
            text = text.replace(re_tag, tag)
        else:
            text = text.replace(tag, re_tag)

    return text


def fa_Chinese(f_pubmed):
    tree = ET.parse(f_pubmed)
    root = tree.getroot()
    for article in root.iter("PubmedArticle"):
        p = FirstAuthorArticle(xml_element=article)

        if len(p.authors) == 0:
            continue

        if p.title.startswith("Author Correction:"):
            continue

        if len(p.title) == 0:
            continue
        p.title = replace_tags(p.title, reverse=True)

        pubtype_set = set(p.get_pubtypes())

        if len(pubtype_set.intersection(set(['Journal Article', 'Letter']))) == 0:
            continue

        if p.abstract is None:
            continue
        p.abstract = replace_tags(p.abstract, reverse=True)
        authors = list()
        au_email = list()

        Chinese_num = 0
        flag = 0
        for author in p.authors:
            flag += 1
            first = author["firstname"]
            last = author["lastname"]
            if first is None or last is None:
                break
            name = first+" "+last
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
                email = search_email[au_email.index(name)]
            elif len(set(search_email)) == 1:
                email = search_email[0]
                affiliation = affiliation.replace(email, "")
                affiliation = affiliation.replace(" Electronic address:", "")
            else:
                email = None
            if email is not None and email[-1] == '.':
                email = email[:-1]

            print(name+'\n'+str(email)+'\n'+affiliation+'\n')

            authors.append([name, first, last, chinese, affiliation, email])

        if Chinese_num > 0:
            topic = get_topic(p.journal, p.abstract)
            articles = [p.pubmed_id[:8], p.publication_date, p.title, p.abstract, p.doi, topic]
            yield authors, articles


class FirstAuthorArticle(PubMedArticle):
    def get_pubtypes(self):
        path = ".//PublicationType"
        return [
            pubtype.text for pubtype in self.xml.findall(path) if pubtype is not None
        ]


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
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="xml")

        return fa_Chinese(str(handle.read()))
