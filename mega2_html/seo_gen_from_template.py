#!/usr/local/bin/python3

from __future__ import print_function

from optparse import OptionParser
import io, sys, os.path
import re
import pdb

PARSED   = "parsed_wik"
DIR      = "conversions"
TEMPLATE = "template.html"
TEMPLAT2 = "template_mega2.html"
#URME     = "file:///Users/rbaron/mega2/bb/mega2_html"
URME     = "https://watson.hgen.pitt.edu/docs/mega2_html"
INPUTS   = ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF', 'IMPUTE2')

ITEMS     = {}
CITS      = {}
INPUThtml = {}

# Classes:
## Template() class constructor initialized to some text (file) will transform it.
##   All keys defined in .map methods are replaced by their value when writing.  If replacement is
##   a list, insert <ul> {<li>text</li>} ... </ul>

## Table() class constructor initialized to file name. follows table html tags

## HtmlGen class() converts CITS, either url or journal to a formated list of html; one entry
##  per citation; uses template for journal entry & fills with available fields (.mkcite())
## 

# Parse "parsed_wik" file:
## Open "input" file as a "utf8" char file (some of wik[i] may be utf8 characters).
## Ignore comment and blank lines and look at first field:
## if == "cit:"  val = rest of line with white space stripped() and initial&trailing "'" removed
##   use val as a key in CITS hash; corresponding entry is an empty hash
##   subsequent lines define key/value to store in CITS hash entry. key has trailing : removed
##   val is rest of line with white space stripped() and initial&trailing "'" removed

## if == name:  val = rest of line with white space stripped() and initial&trailing "'" removed
##   use val as a key in the ITEMS hash; corresponding entry is an empty array
##   subsequent lines, all elements are added to array with initial&trailing "'" removed

# Data:
## INPUTS (static)     === most ITEMS.keys() that do not end with "format"
## analysis (computed) === ITEMS.keys() ending with "format" NOW not case sensitive
## 
## Names coming from INPUTS and ANALYSIS are used as is for html anchor "link text"
##  For all other uses, Names are converted to lower case, have some embeded punctuation 
##  characters ("/" and " ") altered, and finally analysis have trailing " format" text removed.
##  Note: analysis names must end in "format" for them to be recognized.

# Processing:
#   Create frame_inp_<INPUT>.html; file anchor tag is #inp:
## INPUTS calc:  for each input in INPUTS [Template T2(Template2)]
##   Template.write(options.output + "frame_inp_<input>.html); <input> is lower and " " -> "_"
##    map k/v: **path**/<url>; **anchor**/#inp:<input>; <input> is lower and " " -> "_"
##   Create INPUThtml[<input>] as HtmlGen for citations describing <input>

#   Create frame_ext_<ANALYSIS>.html file; anchor tag is #ext:
## analysis Calc: for each entry in analysis
##   analy is entry; an is first word (to whitespace); anr is an lower and "/" -> "" (remove)
##   Template.write(options.output + "frame_ext_<anr>.html)
##    map k/v: **path**/<url>; **anchor**/#ext:<anr>;
##   Create htmla as HtmlGen for citations describing <analy>
## 
#   Create <INPUT>_<ANALYSIS>.html file
## In same analysis loop (above) with analy/an/anr defined
##   and for each input in INPUTS
##  Define file_name as <inpr>_<anr>.html (i.e. name lower case and slightly editted)
##  file_name goes into Table.col() and "Template" T is expanded and written to file_name.



def parse_info(options):
    if not os.path.isfile(options.input):
        print("Not a File: {0}".format(options.input))
    File = io.open(options.input, "r", encoding='utf8')
    rexp = re.compile("'.*?'")
    for line in File:
        line = line.strip()
        if line == "" or line.startswith("#"):
            continue
        fields = line.split(None, 1)
        key = fields[0]
        if key == 'cit:':
            cit = {}
            name = None
            val = fields[1].strip()[1:-1]
            CITS[val] = cit
        elif key == 'name:':
            name = []
            cit = None
            val = fields[1].strip()[1:-1]
            ITEMS[val] = name
        elif cit is not None:
            val = fields[1].strip()[1:-1]
            cit[key[:-1]] = val
        elif name is not None:
            name.extend((x.strip("'") for x in rexp.findall(line)))
    File.close()

def show_input():
    for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF', 'IMPUTE2'):
        print('\n{0}\n\t{1}\n\t{2}\n\t['.format(inp, ITEMS[inp][1], ITEMS[inp][0]))
#       for ref in ITEMS[inp][0]:
#            if ref in CITS:
#               print('1,' , end="")
#           else:
#               print('0,', end="")
#       print(']\n')

class HtmlGen(object):
    match = None
    repl  = None
    all   = []

    def __init__(self, m, r):
        self.match = m
        self.repl  = r
        self.all   = []

    def __iter__(self):
        return self.all.__iter__()

# cite a url with provided anchor text
    def mega2(self, anchr, urme, tag, anly):
        self.all.append(anchr.format(urme, tag.replace(self.match, self.repl), anly) )

# cite a simple url with default anchor text
    def url(self, ustr):
        self.all.append('  <a class="FlexURL" target="_blank" href="{0}">{1}</a>  '.format(CITS[ustr]['url'], ustr[:-4]) )

# These are all the possible pieces of a doc citation in the accepted order.  If a first element
# is found in the CITS hash for the 'citename', its value is printed using the format string that
# is the second element.  (When the first element is empty, just print the second element.)
    elements = (('author','{0}'), ('vauthors','{0}'), \
                ('last','{0} '),   ('first','{0}'), \
                ('author2', ', {0}'), ('author3', ', {0}'), \
                ('author4', ', {0}'), ('author5', ', {0}'), \
                ('coauthors',', {0}'), ('', '. '),\
                ('date','({0}) '), \
                ('title','{0}. '),  \
                ('journal','{0} '), \
                ('volume','{0}'),  ('issue','({0})'), \
                ('pages',':{0}.'), ('url','<span><a class="FlexURL" target="_blank" href="{0}">{0}</a></span>'), \
                ('doi', ' <span><a class="FlexURL" target="_blank" href="http://doi.org/{0}">doi: {0}</a></span>.'), \
                ('pmid',' <span><a class="FlexURL" target="_blank" href="http://www.ncbi.nlm.nih.gov/pubmed/{0}">PMID: {0}</a></span>.'), \
                ('pmc', ' <span><a class="FlexURL" target="_blank" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC{0}">PMC{0}</a></span>.') )

    def mkcite(self, citename):
        cite = CITS[citename]
        if not cite:
            return
        text = []
        for el in self.elements:
            if el[0] == '':
                text.append(el[1])
            elif el[0] in cite:
                text.append(el[1].format(cite[el[0]]))
        self.all.append("".join(text))

# determine if the citation is a url or journal
    def htmltext(self, cit):
        if cit.endswith('url'):
            self.url(cit)
        else:
            self.mkcite(cit)

class Template(object):
    repl = {}
    href = {}
    lines= []
    MT   = []

    def __init__(self, Filename):
        if not os.path.isfile(Filename):
            print("Not a File: {0}".format(Filename))
#xx     File = io.open(Filename, "r", encoding='utf8')
        File = open(Filename)
        self.lines = File.readlines()
        File.close()

    def map(self, m, r):
        self.repl[m] = r

#   rewrite the stored template lines to a file performing the subsititutions of vals for keys
    def write(self, Filename):
        File = io.open(Filename, "w", encoding='utf8')

        ks = self.repl.keys()

        for line in self.lines:
            newline = line
            href = self.MT

            for k in ks:
                if k in newline:
                    val = self.repl[k]
                    if isinstance(val, str):
                        newline = newline.replace(k, val)
                    else:
#                       pdb.set_trace()
                        newline = newline.replace(k, '<ul >')
                        href = val
            print(newline, file=File, end="")

            for xline in href:
                print('<li>{0}</li>'.format(xline), file=File)
            if href:
                print('</ul>', file=File)

        File.close()


# pretty routine class that mimics/generates the html table tags
class Table(object):
    File   = None
    TabCol = 0
    TabRow = 0

    def __init__(self, Filename):
        self.File = open(Filename, "w")
        self.TabCol = 1
        self.TabRow = 1

    def css(self):
        print("""<!--
  <link rel="stylesheet" href="Mega2_Conversions.css" type="text/css" media="all"/>
-->
<style type="text/css">
table {
	border-collapse: collapse;
}
</style>
<style type="text/css">
th {
	font-weight: bold;
        text-align: center;
        border: 1px solid black;
}
th.ColEven, td.ColEven {
        background: LightGrey;
}
td.ColOdd {
        border: 1px solid black;
        text-align: left;
        margin:  0pt;
        padding: 4pt;
}
td.ColEven {
        border: 1px solid black;
        text-align: left;
        margin:  0pt;
        padding: 4pt;
}
tr.RowOdd {
}
tr.RowEven {
}
</style>
""", file=self.File)

    def head(self):
        print('<table>', file=self.File)
        print('<tr class="RowOdd">', file=self.File)
        print('<th class="ColOdd">LINKAGE</th><th class="ColEven">Mega2</th>', file=self.File)
        print('<th class="ColOdd">PLINK</th><th class="ColEven">VCF or BCF</th>', file=self.File)
        print('<th class="ColOdd">IMPUTE2</th>', file=self.File)
        print('</tr>', file=self.File)
        self.TabRow = self.TabRow + 1

    def row(self):
        print('<tr class="{0}">'.format("RowOdd" if (self.TabRow & 1) else "RowEven"), file=self.File)

    def endrow(self):
        print('</tr>', file=self.File)
        self.TabRow = self.TabRow + 1
        self.TabCol = 1

    def col(self, Filename, analy):
        print('<td class="{0}">'.format("ColOdd" if (self.TabCol & 1) else "ColEven"), file=self.File)
        print('<a href="{0}">to {1}'.format(Filename, analy), file=self.File)

    def endcol(self):
        print('</td>', file=self.File)
        self.TabCol = self.TabCol + 1

    def close(self):
        print('</table>', file=self.File)
        self.File.close()


def main():
    MT = []
    parse = OptionParser()
    parse.add_option("-i", "--input", action="store", default=PARSED,
                     help="input file as wikipedia info");
    parse.add_option("-o", "--output", action="store", default=DIR,
                     help="output file many htmls into DIR");
    parse.add_option("-u", "--urlbase", "--url", action="store", default=URME,
                     help="base url directory");

    options, args = parse.parse_args()


    parse_info(options)
    analysis = [ k for k in sorted(ITEMS.keys()) if k.lower().endswith("format") ]

    T  = Template(TEMPLATE)
    T2 = Template(TEMPLAT2)

    TabFile = Table('table.html')
    TabFile.css()
    TabFile.head()

#   show_input()
#   Create frame_inp_<INPUT>.html file
    for inp in INPUTS:
        inpr = inp.replace(' ', '_').lower()
        T2.map('**path**', options.urlbase)
        T2.map('**anchor**', "#inp:{0}".format(inpr))
        file_name = options.output + "/frame_inp_{0}.html".format(inpr)
        T2.write(file_name)

        htmli = HtmlGen(' ', '_')
# First citation will be pointer to Mega2 Input types documentation
        href_file = "frame_inp_{0}.html".format(inpr)
        htmli.mega2('<a class="FlexURL" target="_blank" href="{0}"> Mega2 Input documentation: {2}</a>',
                    href_file, inpr, inp)

# Add the rest of the citations for an input type
        for cit in ITEMS[inp]:
            htmli.htmltext(cit)
        INPUThtml[inp] = htmli

#    pdb.set_trace()
    for analy in analysis:
#   Create frame_ext_<ANALYSIS>.html file
        an  = analy.split()[0]
        anr = an.replace('/', '').lower()
        print('\n{0}, "{1}"\n\t{2}\n\n'.format(an, analy, ITEMS[analy]))

        T2.map('**path**', options.urlbase)
        T2.map('**anchor**', "#ext:{0}".format(anr))
        file_name = options.output + "/frame_ext_{0}.html".format(anr)
        T2.write(file_name)

#   Create <INPUT>_<ANALYSIS>.html file
        htmla = HtmlGen('/', '')   # / disappears

# First citation will be pointer to Mega2 Analysis documentation
        href_file = "frame_ext_{0}.html".format(anr)
        htmla.mega2('<a class="FlexURL" target="_blank" href="{0}"> Mega2 Analysis documentation: {2}</a>',
                    href_file, anr, analy)

# Add the rest of the citations for an analysis
        for cit in ITEMS[analy]:
            htmla.htmltext(cit)

# New row in Table of "analysis" X "inputs" for the analysis
        TabFile.row()

# Loop over inputs for each analysis; generate col in table and file with html
        for inp in INPUTS:
            inpr = inp.replace(' ', '_').lower()  # " " becomes _

            file_name = options.output + "/{0}_{1}.html".format(inpr, anr)

# Add a new column to the analysis row for each input
            TabFile.col(file_name, analy)

# Now take the main template and rewrite it using the given analysis, input and citations
            T.map('**input**', inp)
            T.map('**analysis**', analy)
            T.map('**ihref**', INPUThtml[inp])
            T.map('**ahref**', htmla)
            T.write(file_name)

            TabFile.endcol()

# row is done
        TabFile.endrow()
# table is done
    TabFile.close()


main()
