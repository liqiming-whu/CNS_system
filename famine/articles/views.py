from rest_framework import viewsets
from .serializers import AuthorSerializer, ArticleSerializer, Author_ArticleSerializer
from .models import Author, Journal, Article, Author_Article


class ArticleViewSet(viewsets.ModelViewSet):
    queryset = Article.objects.all()
    serializer_class = ArticleSerializer


class AuthorViewSet(viewsets.ModelViewSet):
    queryset = Author.objects.all()
    serializer_class = AuthorSerializer
