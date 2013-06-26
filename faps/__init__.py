"""

Faps is a general purpose atomisitic simulation toolbox with an emphasis on
sorption analysis.

"""

# Turn on keyword expansion to get revision numbers in version strings
# in .hg/hgrc put
# [extensions]
# keyword =
#
# [keyword]
# faps.py =
#
# [keywordmaps]
# Revision = {rev}

# If there is no keyword expansion then just fall back to 0 for final digit.
try:
    __version_info__ = (2, 0, 0, int("$Revision$".strip("$Revision: ")))
except ValueError:
    __version_info__ = (2, 0, 0, 0)

# Keep this standard so that faps.__version__ will give the current version!
__version__ = "%i.%i.%i.%i" % __version_info__
