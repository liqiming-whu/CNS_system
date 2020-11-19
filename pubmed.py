#!/usr/bin/env python3
import re
import os
import csv
import pandas as pd
from Bio import Entrez, Medline
from time import perf_counter

Entrez.email = "liqiming1914658215@gmail.com"
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"


class Search_Pubmed:

    def __init__(self, keywords, retmax=10000):
        self.keywords = keywords
        self.retmax = retmax
        self.count = Search_Pubmed.get_count(keywords)
        self.idlist = Search_Pubmed.search(keywords, retmax)
        print(self)

    def __repr__(self):
        return "Search '{}', get {} results.".format(self.keywords, self.count)

    __str__ = __repr__

    @staticmethod
    def get_count(term):
        handle = Entrez.egquery(term=term)
        record = Entrez.read(handle)
        for row in record["eGQueryResult"]:
            if row["DbName"] == "pubmed":
                count = row["Count"]
        return count

    @staticmethod
    def search(term, retmax):
        handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        record = Entrez.read(handle)
        return record["IdList"]

    @staticmethod
    def detail(idlist):
        t1 = perf_counter()
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
        t2 = perf_counter()
        xml = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="xml")
        t3 = perf_counter()
        # print(t2-t1, t3-t2)

        print(handle.read().split('\n\n'))

        # records = Medline.parse(handle)
        # print(str(handle.read()))
        # with open("text.txt", "w") as f:
        #     f.write(handle.read())
        # with open(name, "w", newline="", encoding="utf-8") as f:
        #     writer = csv.writer(f)
        #     writer.writerow(["Pubmed_ID", "Title", "Author", "Source", "Url"])
        #     for record in records:
        #         print(record)
        #         pmid = str(record.get("PMID", "?"))
        #         title = record.get("TI", "?")
        #         author = ", ".join(record.get("AU", "?"))
        #         source = record.get("SO", "?")
        #         url = "https://pubmed.ncbi.nlm.nih.gov/" + pmid
        #         writer.writerow([pmid, title, author, source, url])

    @staticmethod
    def detail_to_excel(idlist, name):
        handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        id_l = []
        title_l = []
        author_l = []
        date_l = []
        journal_l = []
        url_l = []
        for record in records:
            print(record)
            pmid = str(record.get("PMID", "?"))
            title = record.get("TI", "?")
            author = ", ".join(record.get("AU", "?"))
            date = record.get("DP", "?").split()[0]
            journal = record.get("TA", "?")
            url = "https://pubmed.ncbi.nlm.nih.gov/" + pmid
            hplink = '=HYPERLINK("{}")'.format(url)
            id_l.append(pmid)
            title_l.append(title)
            author_l.append(author)
            journal_l.append(journal)
            date_l.append(date)
            url_l.append(hplink)

        data = {
            "Article": title_l,
            "PubmedID": id_l,
            "作者": author_l,
            "期刊": journal_l,
            "年份": date_l,
            "文章pubmed链接": url_l
        }
        df = pd.DataFrame(data)
        df.to_excel(name)

    def get_records(self):
        handle = Entrez.efetch(db="pubmed", id=self.idlist, rettype="medline", retmode="text")
        records = Medline.parse(handle)

        return records

    def save_abstract(self, path):
        if not os.path.exists(path):
            os.mkdir(path)
        records = self.get_records()
        for record in records:
            try:
                pmid = str(record.get("PMID", "?"))
                print("save ", pmid)
                abstract = record.get("AB", "?")
                with open(path+"/"+pmid+".txt", "w", encoding="utf-8") as f:
                    f.write(abstract)
            except Exception:
                continue

    def save_idlist(self, name):
        with open(name, "w", encoding="utf-8") as f:
            f.write("\n".join(self.idlist))     

    def save_detail_to_csv(self, name):
        records = self.get_records()
        with open(name, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Pubmed_ID", "Title", "Author", "Source", "Url"])
            for record in records:
                try:
                    pmid = str(record.get("PMID", "?"))
                    title = record.get("TI", "?")
                    author = ", ".join(record.get("AU", "?"))
                    source = record.get("SO", "?")
                    url = "https://pubmed.ncbi.nlm.nih.gov/" + pmid
                    writer.writerow([pmid, title, author, source, url])
                except Exception:
                    continue

    def save_detail_to_excel(self, name):
        id_l = []
        title_l = []
        author_l = []
        source_l = []
        url_l = []
        records = self.get_records()
        for record in records:
            try:
                pmid = str(record.get("PMID", "?"))
                title = record.get("TI", "?")
                author = ", ".join(record.get("AU", "?"))
                source = record.get("SO", "?")
                url = "https://pubmed.ncbi.nlm.nih.gov/" + pmid
                hplink = '=HYPERLINK("{}")'.format(url)
                id_l.append(pmid)
            except Exception:
                continue
            title_l.append(title)
            author_l.append(author)
            source_l.append(source)
            url_l.append(hplink)

        data = {
            "Pubmed_ID": id_l,
            "Title": title_l,
            "Author": author_l,
            "Source": source_l,
            "link": url_l
        }
        df = pd.DataFrame(data)
        df.to_excel(name)

    def save_info(self, csvfile):
        records = self.get_records()
        with open(csvfile, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["Pubmed_ID", "Title", "Author", "Source", "Abstract"])
            for record in records:
                try:
                    pmid = str(record.get("PMID", "?"))
                    title = record.get("TI", "?")
                    author = ", ".join(record.get("AU", "?"))
                    source = record.get("SO", "?")
                    abstract = record.get("AB", "?")
                    writer.writerow([pmid, title, author, source, abstract])
                except Exception:
                    print("Save failed. record:", record)
                    continue


if __name__ == "__main__":
    idlist = open("pmids.txt").read().split()
    # Search_Pubmed.detail_to_excel(idlist, "articles.xlsx")
    Search_Pubmed.detail(idlist)