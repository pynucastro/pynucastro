import argparse
import re
import requests


PR_URL_BASE = r"https://github.com/pynucastro/pynucastro/pull/"
RELEASE_URL_BASE = r"https://github.com/pynucastro/pynucastro/releases/tag/"


pr = re.compile(r"(\#)(\d+)")
tag = re.compile(r"(\#\# )(\d\.\d+\.\d+)")

def doit(clfile):

    with open(clfile) as cl:
        for line in cl:
            if g := tag.match(line):
                url = rf"{RELEASE_URL_BASE}{g.group(2)}"
                response = requests.get(url)
                if response.status_code < 400:
                    new_line = re.sub(tag, rf"## [\g<2>]({RELEASE_URL_BASE}\g<2>)", line)
                    print(new_line)
                    continue

            new_line = re.sub(pr, rf"[\g<0>]({PR_URL_BASE}\g<2>)", line)
            print(new_line, end="")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("changelog", type=str, nargs=1,
                        help="ChangeLog file")

    args = parser.parse_args()

    doit(args.changelog[0])
