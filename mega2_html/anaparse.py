#!/usr/bin/python

from __future__ import print_function

import sys, os.path
import re
import pdb

ANALYSIS = "wik"
TEMPLATE = "template.html"
TEMPLAT2 = "template_mega2.html"
#URME     = "file:///Users/rbaron/mega2/bb/mega2_html"
URME     = "http://watson.hgen.pitt.edu/docs/mega2_html"

HEAD = re.compile("""\|\s*(.*?)<ref""")
HED  = re.compile("""\|\s*(.*?)\s*\|\|""")
REF  = re.compile("""(<ref.*?</ref>)|(<ref.*?/>)""")
CITE = re.compile("""<ref(.*?)/?>{{(.*?)}}""")
REUS = re.compile("""<ref(.*?)/?>(\s*)""")

URL  = re.compile("""\[(.*?)\]""")

DB   = {}
DBREF= {}

def parse(FILE):
#   pdb.set_trace()
    if not os.path.isfile(FILE):
        print("Not a File: {0}".format(FILE))
    File = open(FILE)
    for line in File:
        line = line.strip()
        if line == "": continue
        # important lines begin |
        if line.startswith("||"):    # url
            pass
        elif line.startswith("|"):   # citation
            reflist = []
            urllist = []
            m = HEAD.search(line)    # | <header name> then "<ref ..."
            if m is None:
                m = HED.search(line) # | <header name> then ... then "||"
                if m is None:
                    print("No Name: {0}".format(line))
                    continue
                else:    # no citation
                    head = m.group(1)
                    DB[head] = (reflist, urllist)
                    print("No citation: {0}".format(head))
            else:
                head = m.group(1)
                headcnt = 0
                DB[head] = (reflist, urllist)
                for mi in REF.finditer(line):        # <ref>(full)</ref> vs <ref (reuse) />
                    headcnt = headcnt + 1
                    if mi.group(2) is None:
                        m = CITE.match(mi.group(1))  # define a full citation
                    else:
                        m = REUS.match(mi.group(2))  # reuse/name an old citation
                    if m is None:
                        print("Bad citation: {0}".format(ref))
                        continue
                    refhash = {}
                    refn = m.group(1).split("=")     # define reference name
                    if len(refn) == 2:
                        key = refn[1].strip()
                    elif headcnt == 1:               # use 1st word of head as reference name
                        key = head.split()[0] + 'Ref'
                    else:                            # and now there are multiple unnamed refs
                        key = head.split()[0]+ str(headcnt) + 'Ref'
                    if key not in DBREF:
                        DBREF[key] = refhash
                    reflist.append(key)
                    ref = m.group(2)
                    fields = ref.split("|")          # | splits fields
                    for field in fields:
                        field = field.strip()
                        kv = field.split("=", 1)
                        if len(kv) < 2:
                            print("{0} discarding {1}".format(head, kv))
                            continue
                        refhash[kv[0]] = kv[1]
        # look for urls in [...]
        for m in URL.finditer(line):
            url = m.group(1)
            if url.startswith("["): continue
            fields = url.split(None, 1)
            DB[head][1].append(fields)

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
        self.all.append('  <a href="{0[0]}">{0[1]}</a>  '.format(ustr) )

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
        File = open(Filename, "w")

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

    parse(ANALYSIS)
    analysis = [ k for k in sorted(DB.iterkeys()) if k.endswith("format") ]

    T  = Template(TEMPLATE)
    T2 = Template(TEMPLAT2)

    TabFile = Table('table.html')
    TabFile.css()
    TabFile.head()

#   show_input()
    for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
        inp = inp.replace(' ', '_').lower()
        T2.map('**path**', URME)
        T2.map('**anchor**', "#inp:{0}".format(inp))
        file_name = "conversions/frame_inp_{0}.html".format(inp)
        T2.write(file_name)

    for analy in analysis:
        an  = analy.split()[0]
        anr = an.replace('/', '').lower()

        T2.map('**path**', URME)
        T2.map('**anchor**', "#ext:{0}".format(anr))
        file_name = "conversions/frame_ext_{0}.html".format(anr)
        T2.write(file_name)

        ahref = Citation('/', '')

#       ahref.mega2('<a href="{0}#ext:{1}"> Mega2 Analysis documentation: {2}</a>',
#                   URME, an.lower(), analy)
        file_name = "frame_ext_{0}.html".format(anr)
        ahref.mega2('<a href="{0}"> Mega2 Analysis documentation: {2}</a>',
                    file_name, anr, analy)

        for url in DB[analy][1]:
            ahref.url(url)

        for citename in DB[analy][0]:
            ahref.mkcite(citename)

        print('\n{0}, "{1}"\n\t{2}\n\t{3}\n'.format(an, analy, DB[analy][1], DB[analy][0]))

#       for ref in DB[analy][0]:
#           if ref in DBREF:
#               print('1,' , end="")
#           else:
#               print('0,', end="")
#       print(']\n')

        TabFile.row()

        for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
            inpr = inp.replace(' ', '_').lower()

            ihref = Citation(' ', '_')
#           ihref.mega2('<a href="{0}#inp:{1}"> Mega2 Input documentation: {2}</a>',
#                       URME, inp.lower(), inp)
            file_name = "frame_inp_{0}.html".format(inpr)
            ihref.mega2('<a href="{0}"> Mega2 Input documentation: {2}</a>',
                        file_name, inpr, inp)

            for url in DB[inp][1]:
                ihref.url(url)

            for citename in DB[inp][0]:
                ihref.mkcite(citename)

            file_name = "conversions/{0}_{1}.html".format(inpr, anr)

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
