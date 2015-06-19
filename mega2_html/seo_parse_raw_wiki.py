#!/usr/local/bin/python3

from __future__ import print_function

from optparse import OptionParser
import io, sys, os.path
import re
import pdb

ANALYSIS = "wik"
PARSED   = "parsed_wik"

TOPIC = re.compile("""\|\s*(.*?)<ref""")
TOPC  = re.compile("""\|\s*(.*?)\s*\|\|""")

REF  = re.compile("""(<ref.*?)(</ref>|/>)""")
CITE = re.compile("""<ref(.*?)/?>{{(.*?)}}""")
REUS = re.compile("""<ref(.*?)/?>(\s*)""")

URL  = re.compile("""\[(.*?)\]""")

DB   = {}
DBREF= {}

def get_citations(topic, line):
    abbrev = topic.split()[0]
    cnt = 0
    for mi in REF.finditer(line):        # <ref>(full)</ref> vs <ref (reuse) />
        cnt = cnt + 1
        if mi.group(2) == '</ref>':
            m = CITE.match(mi.group(1))  # define a full citation
            fresh = True
        else:
            m = REUS.match(mi.group(0))  # reuse/name an old citation
            fresh = False
        if m is None:
            print("Bad citation: {0}".format(ref))
            continue

        refclause = m.group(1).split("=")     # define reference name
        if len(refclause) >= 2:
            refkey = refclause[-1].strip() + '_cit'
            if refkey in DBREF and fresh:
                cnt += 1
                refkey = abbrev + str(cnt) + '_cit'
        elif cnt == 1:               # use 1st word of topic as reference name
            refkey = abbrev + '_cit'
        else:                            # and now there are multiple unnamed refs
            refkey = abbrev + str(cnt) + '_cit'
#           cnt += 1

        refhash = {}
        if refkey not in DBREF:
            DBREF[refkey] = refhash
        DB[topic][0].append(refkey)

        ref = m.group(2)
        if ref is None or ref == "":
            continue

        fields = ref.split("|")          # | splits fields
        for field in fields:
            field = field.strip()
            kv = field.split("=", 1)
            if len(kv) < 2:
                if field != 'cite journal':  # no =; but this is how a citation begins
                    print("{0} discarding cit element \"{1}\"".format(topic, field))
                continue
            refhash[kv[0].rstrip()] = kv[1].lstrip()

def get_urls(topic, line):
    for m in URL.finditer(line):
        url = m.group(1)
        if url.startswith("["): continue
        fields = url.split(None, 1)
        fields.reverse()
        fields[0] = fields[0] + '_url'
        DBREF[fields[0]] = fields[1]
        DB[topic][1].append(fields[0])

def parse_wiki(FILE, Parsed):
    if not os.path.isfile(FILE):
        print("Not a File: {0}".format(FILE))
    File = io.open(FILE, "r", encoding='utf8')
    On = False
    topic = None
    for line in File:
        if not On:
            if line.startswith("'''Mega2, the Manipulation Environment for Genetic Analysis'''"):
                topic = "Mega2 Prologue"
                DB[topic] = ([], [])  # no access to url list
                get_citations(topic, line)
                continue
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

        if line.startswith("|-"):   # ??
            pass

        elif line.startswith("||"):   # url
            get_urls(topic, line)     # topic on previous line (or here)

        elif line.startswith("|"):    # citation

            m = TOPIC.search(line)    # viz. | <topic name> "<ref ..."
            if m:
                topic = m.group(1)
                DB[topic] = ([], [])
                get_citations(topic, line)
                get_urls(topic, line)
            else:
                m = TOPC.search(line)  # viz. | <topic name> ... "||"
                if m:
                    topic = m.group(1)
                    DB[topic] = ([], [])
                    print("No citation: {0}, viz \"{1}\"".format(topic, line))  # tag only has urls
                    get_urls(topic, line)
                else:
                    if line != '|}':   # seems to happen occasionally
                        print("No Name: {0}".format(line))

def write_DBREF(fd):
    for k in sorted(DBREF.keys()):
        print("{0}: {1!r}".format('cit', k), file=fd)
        v = DBREF[k]
        if type(v) == str:
            print("    url: {0!r}\n".format(v), file=fd)
        else:
            for cik in sorted(v.keys()):
                civ = v[cik]
                print("    {0}: {1!r}".format(cik, civ), file=fd)
            print("", file=fd)

def write_DB(fd):
    for k in sorted(DB.keys()):
        print("{0}: {1!r}".format('name', k), file=fd)

        for e in DB[k][1]:  # url
            print("{0!r} ".format(e), end='', file=fd)

        for e in DB[k][0]:  # citation
            print("{0!r} ".format(e), end='', file=fd)

        print("\n", file=fd)

def main():
    parse = OptionParser()
    parse.add_option("-i", "--input", action="store", default=ANALYSIS,
                     help="input file as raw wikipedia page");
    parse.add_option("-o", "--output", action="store", default=PARSED,
                     help="output file as more legible page");

    options, args = parse.parse_args()

    parse_wiki(options.input, options.output)

    fd = io.open(options.output, "w", encoding='utf8')
    write_DBREF(fd)
    write_DB(fd)
    fd.close()

main()
