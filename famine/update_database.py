#!/usr/bin/env python3
import os
from datetime import date, timedelta
import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'famine.settings')
django.setup()
from pubmed.search import Pubmed, journal_if
from articles.models import Article, Author, Author_Article, Journal


def update(journal, mindate=None):
    today = date.today().strftime("%Y/%m/%d")
    fromdate = mindate if mindate else (date.today() - timedelta(days=30)).strftime("%Y/%m/%d")
    journal_object = Journal.objects.get(name=journal)
    for authors, article in Pubmed(journal, fromdate, today).detail():
        print(article[0], article[1])
        Article_obj = Article.objects.filter(pmid=article[0])
        if Article_obj:
            continue
        else:
            Article.objects.create(
                pmid=article[0],
                journal=journal_object,
                pubdate=article[1],
                title=article[2],
                abstract=article[3],
                doi=article[4],
                subject=article[5])

        Article_object = Article.objects.get(pmid=article[0])
        rank = 1
        for author in authors:
            Author_obj = Author.objects.filter(
                name=author[0], affiliation=author[4], email=author[5])
            if not Author_obj:
                Author.objects.create(
                    name=author[0],
                    first=author[1],
                    last=author[2],
                    chinese=author[3],
                    affiliation=author[4],
                    email=author[5])
            Author_object = Author.objects.get(
                name=author[0], affiliation=author[4], email=author[5])
            try:
                if rank == 1:
                    Author_Article.objects.create(
                        author=Author_object, article=Article_object, rank=rank, first='Yes')
                else:
                    Author_Article.objects.create(
                        author=Author_object, article=Article_object, rank=rank)
            except Exception:
                continue
            rank += 1


if __name__ == "__main__":
    for jour_if in journal_if:
        print(jour_if[0], jour_if[1])
        Journal.objects.get_or_create(name=jour_if[0], ifactor=jour_if[1])

    journals = ("Nature", "Science", "Cell")
    for journal in journals:
        update(journal)