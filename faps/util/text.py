#!/usr/bin/env python

"""
Small text processing utilities.

"""

import re

def unique(in_list, key=None):
    """Unique values in list ordered by first occurance"""
    uniq = []
    if key is not None:
        keys = []
        for item in in_list:
            item_key = key(item)
            if item_key not in keys:
                uniq.append(item)
                keys.append(item_key)
    else:
        for item in in_list:
            if item not in uniq:
                uniq.append(item)
    return uniq

def ufloat(text):
    """Convert string to float, ignoring the uncertainty part."""
    return float(re.sub('\(.*\)', '', text))


def gfloat(text):
    """Parse a gulp output, where float could also be a rational fraction."""
    # This is much faster than fractions.Fraction
    if "/" in text:
        text = text.split('/')
        return float(text[0])/float(text[1])
    else:
        return float(text)


def try_int(text, default=0):
    """Try to parse an integer but return a default if it fails."""
    try:
        return int(text)
    except ValueError:
        return default


def strip_blanks(lines):
    """Strip lines and remove blank lines."""
    return [line.strip() for line in lines if line.strip() != '']
