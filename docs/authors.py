"""Output the authors.rst for the documentation"""

import subprocess


IGNORE = ["dependabot[bot]",
          "Travis CI User"]


def _surname_key(name):
    """Attempt to sort by surname, this will catch middle initials
    (with a ".")  and skip them, but treat two full names as a single
    surname

    """

    words = name.split()
    start = 1 if len(words) > 1 else 0

    while start < len(words) - 1:
        word = words[start]
        if len(word) != 2 or not word[0].isalpha() or word[1] != ".":
            break
        start += 1

    return " ".join(words[start:]).casefold()


def get_authors(n):
    """Get the active authors (from last n years) and other authors

    Parameters
    ----------
    n : int
        Number of years in the past to now for active authors

    Returns
    -------
    list(str), list(str)
        Active authors and other authors

    """

    assert isinstance(n, int) and n > 0

    # run git shortlog for last n nears

    result = subprocess.run(
        ["git", "shortlog", "-sn", f"--since={n} years ago", "HEAD"],
        check=True,
        capture_output=True,
        text=True)

    output = result.stdout.splitlines()

    active = []
    for line in output:
        if line.strip() == "":
            continue
        _, name = line.strip().split(maxsplit=1)
        if name in IGNORE:
            continue
        active.append(name)

    active = sorted(active, key=_surname_key)

    # run git shortlog for last n nears

    result = subprocess.run(
        ["git", "shortlog", "-sn", "HEAD"],
        check=True,
        capture_output=True,
        text=True)

    output = result.stdout.splitlines()

    # now of all time

    remaining = []
    for line in output:
        _, name = line.strip().split(maxsplit=1)
        if name in IGNORE or name in active:
            continue
        remaining.append(name)

    remaining = sorted(remaining, key=_surname_key)

    return active, remaining


def write_markdown(n, active, remaining):

    with open("source/authors.rst", "w") as of:

        heading = "Authors"
        of.write(f"{heading}\n")
        of.write("=" * len(heading) + "\n\n")

        heading = f"Current authors (active in the last {n} years)"
        of.write(f"{heading}\n")
        of.write("-" * len(heading) + "\n\n")
        for name in active:
            of.write(f"* {name}\n")
        of.write("\n")

        heading = "Previous authors"
        of.write(f"{heading}\n")
        of.write("-" * len(heading) + "\n\n")
        for name in remaining:
            of.write(f"* {name}\n")
        of.write("\n")


if __name__ == "__main__":

    n = 3
    active, remaining = get_authors(n)
    write_markdown(n, active, remaining)

