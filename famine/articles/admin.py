from django.contrib import admin

from django.contrib import admin
from .models import Article, Author, Author_Article, Journal


admin.site.register(Article)
admin.site.register(Author)
admin.site.register(Author_Article)
admin.site.register(Journal)
