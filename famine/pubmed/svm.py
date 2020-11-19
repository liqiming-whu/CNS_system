#!/usr/bin/env python3
import os
import numpy as np
import joblib
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC
from data_utils import build_word_dict, build_svm_dataset


class Svm(object):
    def __init__(self):
        self.word_dict = build_word_dict()
        self.vocabulary_size = len(self.word_dict)
        self.train_x, self.train_y = build_svm_dataset("train", self.word_dict)
        self.test_x, self.test_y = build_svm_dataset("test", self.word_dict)

    def make_tfidf(self):
        vectorizer = CountVectorizer()
        tfidftransformer = TfidfTransformer()
        tfidf = tfidftransformer.fit_transform(
            vectorizer.fit_transform(self.train_x))

        return tfidf

    def train(self):
        text_clf = Pipeline([("vect", CountVectorizer()),
                             ("tfidf", TfidfTransformer()),
                             ("clf", SVC(C=1, kernel="linear"))])
        text_clf = text_clf.fit(self.train_x, self.train_y)
        print("train finished.")

        return text_clf

    def predict(self, text_clf):
        predicted = text_clf.predict(self.test_x)
        accuracy = np.mean(predicted == self.test_y)

        return accuracy

    def save_model(self, text_clf):
        joblib.dump(text_clf, os.path.join("data", "train_model.m"))

    def main(self):
        tfidf = self.make_tfidf()
        print("tfidf:", tfidf.shape)
        text_clf = self.train()
        self.save_model(text_clf)
        acc = self.predict(text_clf)

        print("Accuracy:", acc)


if __name__ == "__main__":
    Svm().main()
