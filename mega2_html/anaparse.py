#!/usr/bin/python

from __future__ import print_function

import sys, os.path
import re
import pdb

ANALYSIS = "wik"
TEMPLATE = "template.html"
#URME     = "file:///Users/rbaron/mega2/bb/mega2_html/mega2_documentation.html"
URME     = "http://watson.hgen.pitt.edu/docs/mega2_html/mega2_documentation.html"


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

elements = (('title','{0}. '), ('author','{0}. '), \
            ('last','{0} '),   ('first','{0}, '), ('coauthors','{0}. '), \
            ('date','({0}) '), ('journal','{0} '), \
            ('volume','{0}'),  ('issue','({0})'), \
            ('pages',':{0}.'), ('url','<span><a href="{0}">{0}</a></span>'))
def mkcit(href, citname):
    cit = DBREF[citname]
    text = []
    for el in elements:
        if el[0] in cit:
            text.append(el[1].format(cit[el[0]]))
    href.append("".join(text))

def main():
    MT = []

    parse(ANALYSIS)

    def prune(el):
        return el.endswith("format")
    analysis = filter(prune, sorted(DB.iterkeys()))
    
    if not os.path.isfile(TEMPLATE):
        print("Not a File: {0}".format(TEMPLATE))
    File = open(TEMPLATE)
    Lines = File.readlines()
    File.close()

    TabFile = open('table.html', "w")
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
""", file=TabFile)
    print('<table>', file=TabFile)
    print('<tr class="RowOdd">', file=TabFile)
    print('<th class="ColOdd">LINKAGE</th><th class="ColEven">Mega2</th>', file=TabFile)
    print('<th class="ColOdd">PLINK</th><th class="ColEven">VCF or BCF</th>', file=TabFile)
    print('</tr>', file=TabFile)
    TabCol = 1
    TabRow = 2
#   show_input()

    for analy in analysis:
        an = analy.split()[0]

        ahref = []
        ahref.append( '<a href="{0}#ext:{1}"> Mega2 Analysis documentation: {2}</a>'.format(URME, an.lower().replace('/', ''), analy) )
        for url in DB[analy][1]:
            ahref.append('  <a href="{0[0]}">{0[1]}</a>  '.format(url) )
        for citname in DB[analy][0]:
            mkcit(ahref, citname)
        print('\n{0}, "{1}"\n\t{2}\n\t{3}\n'.format(an, analy, DB[analy][1], DB[analy][0]))

#       for ref in DB[analy][0]:
#           if ref in DBREF:
#               print('1,' , end="")
#           else:
#               print('0,', end="")
#       print(']\n')

        print('<tr class="{0}">'.format("RowOdd" if (TabRow & 1) else "RowEven"), file=TabFile)
        TabRow = TabRow + 1
        for inp in ('LINKAGE', 'Mega2', 'PLINK', 'VCF or BCF'):
            ihref = []
            ihref.append( '<a href="{0}#inp:{1}"> Mega2 Input documentation: {2}</a>'.format(URME, inp.lower().replace(' ', '_'), inp) )
            for url in DB[inp][1]:
                ihref.append('  <a href="{0[0]}">{0[1]}</a>  '.format(url) )
            for citname in DB[inp][0]:
                mkcit(ihref, citname)

            file_name = "conversions/{0}_{1}.html".format(inp.replace(' ', '_'), an.replace('/',''))
            file_name = file_name.lower()

            print('<td class="{0}">'.format("ColOdd" if (TabCol & 1) else "ColEven"), file=TabFile)
            TabCol = TabCol + 1
            print('<a href="{0}">to {1}'.format(file_name, analy), file=TabFile)

            File = open(file_name, "w")
            
            for line in Lines:
                newline = line
                href = MT

                if newline.find('**input**') > -1:
                    newline = newline.replace('**input**', inp)

                if newline.find('**analysis**') > -1:
                    newline = newline.replace('**analysis**', analy)

                if newline.find('**ihref**') > -1:
                    newline = newline.replace('**ihref**', '<ul >')
                    href = ihref

                if newline.find('**ahref**') > -1:
                    newline = newline.replace('**ahref**', '<ul >')
                    href = ahref

                print(newline, file=File, end="")

                for xline in href:
                    print('<li>{0}</li>'.format(xline), file=File)
                if href:
                    print('</ul>', file=File)

            print('</td>', file=TabFile)
            File.close()

        print('</tr>', file=TabFile)

    print('</table>', file=TabFile)
    TabFile.close()


main()
