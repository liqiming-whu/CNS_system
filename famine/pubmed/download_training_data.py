#!/usr/bin/env python3
import os
from Bio import Entrez, Medline
import csv
import pandas as pd
import wget
import tarfile

Entrez.email = "liqiming1914658215@gmail.com"
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"


def search(keywords):
    handle = Entrez.esearch(db='pubmed', term=keywords, retmax=2400)
    record = Entrez.read(handle)

    return record["Count"], record["IdList"]


def get_abstract(idlist):
    handle = Entrez.efetch(db='pubmed', id=idlist, rettype="medline",
                           retmode="text")
    records = Medline.parse(handle)

    title_list = []
    abstract_list = []
    for record in records:
        try:
            title = record["TI"]
            abstract = record["AB"]
            title_list.append(title)
            abstract_list.append(abstract)
        except Exception:
            continue

    return title_list, abstract_list


def make_train_data(path, title_list, abstract_list, tag):
    with open(os.path.join(path, "train_data.csv"), "w",
              newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        print("make train data...")
        for title, abstract in zip(title_list[:1800],
                                   abstract_list[:1800]):
            writer.writerow([title, abstract, tag])


def make_test_data(path, title_list, abstract_list, tag):
    with open(os.path.join(path, "test_data.csv"), "w",
              newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        print("make test data...")
        for title, abstract in zip(title_list[1800:2100],
                                   abstract_list[1800:2100]):
            writer.writerow([title, abstract, tag])


def main(path, keywords, tag):
    count, idlist = search(keywords)
    print("count:", count)
    title_list, abstract_list = get_abstract(idlist)
    print("list length:", len(title_list))
    make_train_data(path, title_list, abstract_list, tag)
    make_test_data(path, title_list, abstract_list, tag)
    print("done!")


def get_other_data(path, tag):
    if not os.path.exists("dbpedia_csv.tar.gz"):
        dbpedia_url = 'https://github.com/le-scientifique/torchDatasets/raw/master/dbpedia_csv.tar.gz'
        wget.download(dbpedia_url)
    with tarfile.open("dbpedia_csv.tar.gz", "r:gz") as tar:
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar)

    df = pd.read_csv("dbpedia_csv/test.csv", names=["class", "title", "content"])
    df = df.sample(frac=0.03)
    title = df["title"]
    content = df["content"]

    with open(os.path.join(path, "train_data.csv"), "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        for t, c in zip(title[:1800], content[:1800]):
            writer.writerow([t, c, tag])

    with open(os.path.join(path, "test_data.csv"), "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        for t, c in zip(title[1800:], content[1800:]):
            writer.writerow([t, c, tag])


def merge(path):
    for data in ("train_data.csv", "test_data.csv"):
        with open(data, "w", encoding="utf-8") as f:
            for p in path:
                f.write(open(os.path.join(p, data), encoding="utf-8").read())


if __name__ == "__main__":
    path = [os.path.join("data", i) for i in (
            "biology", "chemistry", "physics", "computer_science", "other")]
    keywords = ("Cell[ta]", "Nat Chem[ta] OR J Am Chem Soc[ta]",
                "J Phys Condens Matter[ta] OR Nat phys[ta] OR Annu Rev Fluid Mech[ta]",
                "IEEE Trans Neural Netw Learn Syst[ta]")
    tags = ["1", "2", "3", "4", "5"]
    for i in range(4):
        if not os.path.exists(path[i]):
            os.mkdir(path[i])
            main(path[i], keywords[i], tags[i])
    if not os.path.exists(path[4]):
        os.mkdir(path[4])
        get_other_data(path[4], tags[4])
    merge(path)
