#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
from datetime import date
from Bio import Entrez


Entrez.email = "liqiming1914658215@gmail.com"                                      
Entrez.api_key = "c80ce212c7179f0bbfbd88495a91dd356708"


class QueryNCBI:
    __slots__ = ['keywords', 'mesh_topic', 'journal', 'year', 'from_date', 'to_date', 'retmax', 'db', 'count', 'idlist']
    def __init__(self, keywords=None, mesh_topic=None, journal=None, year=None, from_date=None, to_date=date.today().strftime("%Y/%m/%d"), retmax=1000, db="pubmed"):
        self.keywords = keywords
        self.mesh_topic = mesh_topic
        self.journal = journal
        self.year = year
        self.from_date = from_date
        self.to_date = to_date
        assert any(getattr(self, i) for i in self.__slots__[:5]), "At least one parameter is required."
        self.retmax = retmax
        self.db = db
        self.count = self.get_count()
        self.idlist = self.search()

    def __repr__(self):
        return f"Search '{self.query}', get {self.count} results."

    __str__ = __repr__

    @property
    def query(self):
        query_list = []
        if self.keywords:
            query_list.append(self.keywords)
        if self.mesh_topic:
            query_list.append(f"{self.mesh_topic}[MeSH Major Topic]")
        if self.journal:
            query_list.append(f"{self.journal}[ta]")
        if self.from_date:
            assert re.compile("\d{4}\/\d{2}\/\d{2}").match(self.from_date), "Date error, fromat: YYYY/MM/DD"
            assert re.compile("\d{4}\/\d{2}\/\d{2}").match(self.to_date), "Date error, fromat: YYYY/MM/DD"
            self.year = None
            query_list.append(f"{self.from_date}: {self.to_date}[dp]")
        if self.year:
            query_list.append(f"{self.year}[dp]")
        
        return " AND ".join(query_list)

    def get_count(self):
        handle = Entrez.egquery(term=self.query)
        record = Entrez.read(handle)
        for row in record["eGQueryResult"]:
            if row["DbName"] == self.db:
                count = row["Count"]
        return count
        
    def search(self):                                                           
       handle = Entrez.esearch(db=self.db, term=self.query, retmax=self.retmax)
       record = Entrez.read(handle)                                            
       return record["IdList"]

    def results_iter(self, type="dict"):
        assert type in ["dict", "xml"], f"format {type} not support."
        id_count = len(self.idlist)
        idlists = [self.idlist[i:i+10000] for i in range(0,id_count, 10000)] if id_count > 10000 else [self.idlist]
        for ids in idlists:
            handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="xml")
            if type == "dict":
                yield Entrez.read(handle)
            else:
                yield handle.read()
    