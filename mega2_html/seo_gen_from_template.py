#!/usr/local/bin/python3

from __future__ import print_function

from optparse import OptionParser
import io, sys, os.path
import re
import pdb

PARSED   = "parsed_wik"
DIR      = "conversionsX"
TEMPLATE = "template.html"
TEMPLAT2 = "template_mega2.html"
URME     = "file:///Users/rbaron/mega2/bb/mega2_html"
#URME     = "https://watson.hgen.pitt.edu/docs/mega2_html"

DB   = {}
DBREF= {}

def parse_info(options):
#   pdb.set_trace()
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
            DBREF[val] = cit
        elif key == 'name:':
            name = []
            cit = None
            val = fields[1].strip()[1:-1]
            DB[val] = name
        elif cit is not None:
            val = fields[1].strip()[1:-1]
            cit[key[:-1]] = val
        elif name is not None:
            name.extend((x.strip("'") for x in rexp.findall(line)))
    File.close()

def show_input():
    for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
        print('\n{0}\n\t{1}\n\t{2}\n\t['.format(inp, DB[inp][1], DB[inp][0]))
#       for ref in DB[inp][0]:
#            if ref in DBREF:
#               print('1,' , end="")
#           else:
#               print('0,', end="")
#       print(']\n')

class Citation(object):
    match = None
    repl  = None
    all   = []

    def __init__(self, m, r):
        self.match = m
        self.repl  = r
        self.all   = []

    def __iter__(self):
        return self.all.__iter__()

    def mega2(self, anchr, urme, tag, anly):
        self.all.append(anchr.format(urme, tag.replace(self.match, self.repl), anly) )

    def url(self, ustr):
        self.all.append('  <a href="{0}">{1}</a>  '.format(DBREF[ustr]['url'], ustr[:-4]) )

    elements = (('title','{0}. '), ('author','{0}. '), \
                ('last','{0} '),   ('first','{0}, '), ('coauthors','{0}. '), \
                ('date','({0}) '), ('journal','{0} '), \
                ('volume','{0}'),  ('issue','({0})'), \
                ('pages',':{0}.'), ('url','<span><a href="{0}">{0}</a></span>'))
    def mkcite(self, citename):
        cite = DBREF[citename]
        text = []
        for el in self.elements:
            if el[0] in cite:
                text.append(el[1].format(cite[el[0]]))
        self.all.append("".join(text))

class Template(object):
    repl = {}
    href = {}
    lines= []
    MT   = []

    def __init__(self, Filename):
        if not os.path.isfile(Filename):
            print("Not a File: {0}".format(Filename))
        File = open(Filename)
        self.lines = File.readlines()
        File.close()

    def map(self, m, r):
        self.repl[m] = r

    def mapcit(self, m, r):
        self.repl[m] = '<ul >'
        self.href[m] = r

    def write(self, Filename):
        File = io.open(Filename, "w", encoding='utf8')

        ks = self.repl.keys()

        for line in self.lines:
            newline = line
            href = self.MT

            for k in ks:
                if k in newline:
                    newline = newline.replace(k, self.repl[k])
                    href = self.href.get(k, href)

            print(newline, file=File, end="")

            for xline in href:
                print('<li>{0}</li>'.format(xline), file=File)
            if href:
                print('</ul>', file=File)
        File.close()


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
        print('</tr>', file=self.File)
        self.TabRow = self.TabRow + 1

    def row(self):
        print('<tr class="{0}">'.format("RowOdd" if (self.TabRow & 1) else "RowEven"), file=self.File)

    def endrow(self):
        print('</tr>', file=self.File)
        self.TabRow = self.TabRow + 1

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

    options, args = parse.parse_args()


    parse_info(options)

    analysis = [ k for k in sorted(DB.keys()) if k.endswith("format") ]

    T  = Template(TEMPLATE)
    T2 = Template(TEMPLAT2)

    TabFile = Table('table.html')
    TabFile.css()
    TabFile.head()

#   show_input()
#   Create frame_inp_<INPUT>.html file
    for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
        inp = inp.replace(' ', '_').lower()
        T2.map('**path**', URME)
        T2.map('**anchor**', "#inp:{0}".format(inp))
        file_name = options.output + "/frame_inp_{0}.html".format(inp)
        T2.write(file_name)

    for analy in analysis:
#   Create frame_ext_<ANALYSIS>.html file
        an  = analy.split()[0]
        anr = an.replace('/', '').lower()

        T2.map('**path**', URME)
        T2.map('**anchor**', "#ext:{0}".format(anr))
        file_name = options.output + "/frame_ext_{0}.html".format(anr)
        T2.write(file_name)

#   Create <INPUT>_<ANALYSIS>.html file
        ahref = Citation('/', '')

#       ahref.mega2('<a href="{0}#ext:{1}"> Mega2 Analysis documentation: {2}</a>',
#                   URME, an.lower(), analy)
        file_name = "frame_ext_{0}.html".format(anr)
        ahref.mega2('<a href="{0}"> Mega2 Analysis documentation: {2}</a>',
                    file_name, anr, analy)

        for cit in DB[analy]:
            if cit.endswith('url'):
                ahref.url(cit)
            else:
                ahref.mkcite(cit)

        print('\n{0}, "{1}"\n\t{2}\n\n'.format(an, analy, DB[analy]))

        TabFile.row()

        for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
            inpr = inp.replace(' ', '_').lower()

            ihref = Citation(' ', '_')
#           ihref.mega2('<a href="{0}#inp:{1}"> Mega2 Input documentation: {2}</a>',
#                       URME, inp.lower(), inp)
            file_name = "frame_inp_{0}.html".format(inpr)
            ihref.mega2('<a href="{0}"> Mega2 Input documentation: {2}</a>',
                        file_name, inpr, inp)

            for cit in DB[inp]:
                if cit.endswith('url'):
                    ihref.url(cit)
                else:
                    ihref.mkcite(cit)

            file_name = options.output + "/{0}_{1}.html".format(inpr, anr)

            TabFile.col(file_name, analy)

            T.map('**input**', inp)
            T.map('**analysis**', analy)
            T.mapcit('**ihref**', ihref)
            T.mapcit('**ahref**', ahref)
            T.write(file_name)

            TabFile.endcol()

        TabFile.endrow()

    TabFile.close()


main()
