#!/usr/local/bin/python3

from __future__ import print_function

from optparse import OptionParser
import io, sys, os.path
import re
import pdb

ANALYSIS = "wik"
PARSED   = "parsed_wik"

HEAD = re.compile("""\|\s*(.*?)<ref""")
HED  = re.compile("""\|\s*(.*?)\s*\|\|""")
#REF  = re.compile("""(<ref.*?</ref>)|(<ref.*?/>)""")
REF  = re.compile("""(<ref.*?)(</ref>|/>)""")
CITE = re.compile("""<ref(.*?)/?>{{(.*?)}}""")
REUS = re.compile("""<ref(.*?)/?>(\s*)""")

URL  = re.compile("""\[(.*?)\]""")

DB   = {}
DBREF= {}

def parse_wiki(FILE, Parsed):
#   pdb.set_trace()
    if not os.path.isfile(FILE):
        print("Not a File: {0}".format(FILE))
    File = io.open(FILE, "r", encoding='utf8')
    On = False
    for line in File:
        Force = False
        if not On:
            if line.startswith("'''Mega2, the Manipulation Environment for Genetic Analysis'''"):
                Force = True
                head = "Mega2 Prologue"
            elif 'Input file formats' in line  and '==' in line:  # start
                On = True
                continue
            else:
                continue
        line = line.strip()
        if line == "": continue
        # important lines begin |
        if 'Mega2 documentation' in line  and '==' in line:       # end
            break
        if line.startswith("||"):    # url
            pass
        elif line.startswith("|-"):   # ??
            pass
        elif line.startswith("|") or Force:   # citation
            reflist = []
            urllist = []
            if not Force:
                m = HEAD.search(line)    # | <header name> then "<ref ..."
                if m is None:
                    m = HED.search(line) # | <header name> then ... then "||"
                    if m is None:
                        if line != '|}':  # seems to happen occasionally
                            print("No Name: {0}".format(line))
                        continue
                    else:    # no citation
                        head = m.group(1)
                        DB[head] = (reflist, urllist)
                        print("No citation: {0}, viz \"{1}\"".format(head, line))  # tag only has urls
                else:
                    head = m.group(1)
            headcnt = 0
            DB[head] = (reflist, urllist)
            for mi in REF.finditer(line):        # <ref>(full)</ref> vs <ref (reuse) />
                headcnt = headcnt + 1
                if mi.group(2) == '</ref>':
                    m = CITE.match(mi.group(1))  # define a full citation
                    fresh = True
                else:
                    m = REUS.match(mi.group(0))  # reuse/name an old citation
                    fresh = False
                if m is None:
                    print("Bad citation: {0}".format(ref))
                    continue
                refhash = {}
                refn = m.group(1).split("=")     # define reference name
                if len(refn) >= 2:
                    key = refn[-1].strip() + '_cit'
                    if key in DBREF and fresh:
                        headcnt += 1
                        key = head.split()[0]+ str(headcnt) + '_cit'
                elif headcnt == 1:               # use 1st word of head as reference name
                    key = head.split()[0] + '_cit'
                else:                            # and now there are multiple unnamed refs
                    key = head.split()[0]+ str(headcnt) + '_cit'
                if key not in DBREF:
                    DBREF[key] = refhash
                reflist.append(key)
                ref = m.group(2)
                if ref is None or ref == "":
                    continue
                fields = ref.split("|")          # | splits fields
                for field in fields:
                    field = field.strip()
                    kv = field.split("=", 1)
                    if len(kv) < 2:
                        if field != 'cite journal':  # no =; but this is how a citation begins
                            print("{0} discarding cit element \"{1}\"".format(head, field))
                        continue
                    refhash[kv[0]] = kv[1]
        # look for urls in [...]
        for m in URL.finditer(line):
            url = m.group(1)
            if url.startswith("["): continue
            fields = url.split(None, 1)
            fields.reverse()
            fields[0] = fields[0] + '_url'
            DBREF[fields[0]] = fields[1]
            DB[head][1].append(fields[0])

    fd = io.open(Parsed, "w", encoding='utf8')
    for k in sorted(DBREF.keys()):
        print("{0}: {1!r}".format('cit', k), file=fd)
        v = DBREF[k]
        if type(v) == str:
            print("    url: {0!r}\n".format(v), file=fd)
        else:
            for cik in sorted(v.keys()):
                civ = v[cik]
                if cik == 'pages':
                    civ = civ.replace(u'\u2013', '-')
                print("    {0}: {1!r}".format(cik, civ), file=fd)
            print("", file=fd)

    for k in sorted(DB.keys()):
        print("{0}: {1!r}".format('name', k), file=fd)
        v = DB[k][1]
        for e in v:
            print("{0!r} ".format(e), end='', file=fd)
        v = DB[k][0]
        for e in v:
            print("{0!r} ".format(e), end='', file=fd)
        print("\n", file=fd)
    fd.close()


def main():
    parse = OptionParser()
    parse.add_option("-i", "--input", action="store", default=ANALYSIS,
                     help="input file as raw wikipedia page");
    parse.add_option("-o", "--output", action="store", default=PARSED,
                     help="output file as more legible page");

    options, args = parse.parse_args()

    parse_wiki(options.input, options.output)


main()
