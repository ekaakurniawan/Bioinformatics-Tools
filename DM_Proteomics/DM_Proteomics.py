# Copyright (C) 2013 by Eka A. Kurniawan
# eka.a.kurniawan(ta)gmail(tod)com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# Data Mining using NLTK (Natural Language Toolkit) for Proteomics

import nltk

sentence = "Here, we present the structure of the MyoX MyTH4-FERM tandem in complex with the cytoplasmic tail P3 domain of the netrin receptor DCC."
# Tokenization
tokens = nltk.word_tokenize(sentence)
print tokens
# Part of sentence tagging
tagged = nltk.pos_tag(tokens)
print tagged
# Grammar for a S (Sentence)
# --------------------------
# Grammar definition of NP (Noun Phrase):
#  - Could start with one or none DT (Determiner)
#  - Followed by as many NNPs (Proper Nouns)
#  - Followed by as many JJs (Adjectives)
#  - Ended by as many NNs or NNPs (Nouns or Proper Nouns)
# Grammar definition of PP (Preposition Phrase):
#  - A preposition followed by a NP (Noun Phrase)
# Grammar definition of VP (Verb Phrase):
#  - One or many verbs followed by NPs or PPs (Noun Phrases or Preposition Phrases)
grammar = r"""
            NP: {<DT>?<NNP>*<JJ>*<NN|NNP>*}
            PP: {<IN><NP>}
            VP: {<VBP>.*<NP|PP>*}
          """
cp = nltk.RegexpParser(grammar)
# Entity detection using predefined grammar
entities = cp.parse(tagged)
print entities
# Draw sentence tree
entities.draw()

