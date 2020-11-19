#!/usr/bin/env python3
import collections
import os
import pickle
import re
import pandas as pd
from nltk.tokenize import word_tokenize
import nltk.data

DIR = os.path.dirname(os.path.abspath(__file__))
NLTK_DATA = os.path.join(DIR, 'data', 'nltk_data')
nltk.data.path.append(NLTK_DATA)

TRAIN_DATA = "train_data.csv"
TEST_DATA = "test_data.csv"


def clean_str(text):
    text = re.sub(r"[^A-Za-z0-9(),!?\'\`\"]", " ", text) 
    text = re.sub(r"\s{2,}", " ", text)
    text = text.strip().lower()

    return text


def build_word_dict():
    if not os.path.exists("word_dict.pickle"):
        train_df = pd.read_csv(TRAIN_DATA, names=["title", "content", "class"])
        contents = train_df["content"]

        words = list()
        for content in contents:
            for word in word_tokenize(clean_str(content)):
                words.append(word)

        word_counter = collections.Counter(words).most_common()
        word_dict = dict()
        word_dict["<pad>"] = 0
        word_dict["<unk>"] = 1
        word_dict["<eos>"] = 2
        for word, _ in word_counter:
            word_dict[word] = len(word_dict)

        with open(os.path.join("data", "word_dict.pikle"), "wb") as f:
            pickle.dump(word_dict, f)

    else:
        with open(os.path.join("data", "word_dict.pikle"), "rb") as f:
            word_dict = pickle.load(f)

    return word_dict


def build_svm_dataset(step, word_dict):
    if step == "train":
        df = pd.read_csv(TRAIN_DATA, names=["title", "content", "class"])
    else:
        df = pd.read_csv(TEST_DATA, names=["title", "content", "class"])

    df = df.sample(frac=1)
    x = list(map(lambda d: word_tokenize(clean_str(d)), df["content"]))
    x = list(map(lambda d: list(
        map(lambda w: word_dict.get(w, word_dict["<unk>"]), d)), x))
    x = list(map(lambda d: d + [word_dict["<eos>"]], x))
    x = list(" ".join(map(str, i)) for i in x)

    y = list(map(lambda d: d - 1, list(df["class"])))

    return x, y


if __name__ == "__main__":
    build_word_dict()