from rest_framework import serializers
from .models import Author, Article, Journal, Author_Article


class AuthorSerializer(serializers.ModelSerializer):

    class Meta:
        model = Author
        fields = ['id', 'name', 'first', 'last', 'chinese', 'affiliation', 'email']


class JournalSerializer(serializers.ModelSerializer):

    class Meta:
        model = Journal
        fields = ['id', 'name', 'first', 'last', 'chinese', 'affiliation', 'email']


class ArticleSerializer(serializers.ModelSerializer):

    class Meta:
        model = Article
        fields = ['pmid', 'journal', 'pubdate', 'title', 'authors', 'abstract', 'doi', 'subject']


class Author_ArticleSerializer(serializers.ModelSerializer):

    class Meta:
        model = Author_Article
        fields = ['id', 'author', 'article', 'rank', 'first', 'cofirst']