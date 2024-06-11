import json
import sys
import re


def auto_links(content, words=["m6A", "FIRE"]):
    """Replace words that are not already markdown links with markdown links to the glossary.md page of the word."""
    for word in words:
        # replace the word with a link to the word if the word is not already within a link
        content = re.sub(r"\.\." + word, f"[{word}](glossary.md#{word.lower()})", content)
    return content


if __name__ == "__main__":
    if len(sys.argv) > 1:  # we check if we received any argument
        if sys.argv[1] == "supports":
            # then we are good to return an exit status code of 0, since the other argument will just be the renderer's name
            sys.exit(0)

    # load both the context and the book representations from stdin
    context, book = json.load(sys.stdin)
    # and now we modify the chapters one by one
    for section in book["sections"]:
        if "Chapter" in section:
            section["Chapter"]["content"] = auto_links(section["Chapter"]["content"])

    # we are done with the book's modification, we can just print it to stdout,
    print(json.dumps(book))
